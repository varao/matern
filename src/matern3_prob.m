function post = matern3_prob(G, side, num_iter)
% MCMC sampler for a stationary Matern type-III repulsive process with probabilistic thinning
% If hardcore = 1
%    Assumes a common radius of repulsion for all points with a non-informative 
%    prior on it
% else
%    each point has its own radius
%
% Within each radius, the thinning probability is a Beta-distributed p_thin

  hardcore = 1;

  % Gamma prior on scale (i.e. Poisson rate)
  sc_alp  = 5;
  sc_beta = 1;

  scale = 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p_thin = .1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dim     = size(G,1);
  num_obs = size(G,2);

  if(length(side) == 1)
    side = side .* ones(dim,1);
  end;

  vXEdge = linspace(0,side(1),5);
  vYEdge = linspace(0,side(2),5);
  mHist2d = hist2d(G',vYEdge,vXEdge) ./ ((vXEdge(2) - vXEdge(1)) * (vYEdge(2) - vYEdge(1)));


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do_fano = 1
  if(do_fano == 1)    % Fano
    for(ll=2:20)
      fano_data(ll-1) = fano(G',side,ll);
    end;
  else             % Ripley
    fano_data = [];
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%

  dist_G                  = sqrt(sq_dist(G,G)); %Since thinning is probabilitic, R's can be arbitrary
  dist_vec                = unique(dist_G);     % Sorted from low to high
  dist_vec                = 1 - exp(-dist_vec); % Prob an exponential falls within [0,a_i) 
  dist_prior              = log([dist_vec(2:end);1]' - dist_vec'); % Log-prob an exponential falls within [a_i, a_i+1)
  dist_hist               = (1:length(dist_vec))-1;    % Number of blocks (-1 because we have both 0 and infty)
                                           % Actually since the distances are unique, this is also the number of observations that weren't thinned.
  dist_G(1:num_obs+1:end) = inf;
  t_G                     = rand(1,num_obs);

  p_vec                   = dist_prior;   % Temp
  num_R                   = 1;

  % Initialize R's randomly
  R = ones(1,num_R) .* min(min(dist_G))/5;

  %%%%%%%%%%%%%%%
     
  vol        = prod(side);
  poiss_mean = scale * vol;

  %%%%%%%%%%%%%%%
  % Debugging
  %%%%%%%%%%%%%%
  count = 0;
  %%%%%%%%%%%%%%%

  for iter = 1:num_iter

    if(rem(iter,500) == 0)
      iter
      min_r
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample thinned events
    num_pts  = poissrnd(poiss_mean);

    if num_pts > 0
      F      = rand(dim,num_pts);
      for(dd = 1:dim)
        F(dd,:) = F(dd,:) * side(dd);
      end;
      t_F    = rand(1,num_pts);
      dist   = sqrt(sq_dist(G,F));
      tim    = t_G(ones(num_pts,1),:)' < t_F(ones(num_obs,1),:);

      shadow = (dist < R) & tim;

      % Some parts of the shadow have higher thinning probability than others
      thin_logp   = sum(shadow) * log(1 - p_thin);  % One minus Shadow, actually (prob_thin = 1 - (1-p_thin)^num)
      thin        = thin_logp < log(rand(1, num_pts));

      F       = F(:,thin);
      t_F     = t_F(:,thin);
      dist    = dist(:,thin);
      tim     = tim (:,thin);
      shadow  = shadow(:,thin);
      num_pts = length(t_F);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      F      = ones(dim,0);
      t_F    = [];
    end;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample times of observations

    for i = 1:num_obs
      if num_pts > 0
        % How many observations (bar i) support each of the thinned events?
        shadow(i,:) = 0;
        shadow_proj = sum(shadow);

        % Focus on thinned events within R of i
        nghbr_F = dist(i,:)   < R;
        t_nF    = t_F(nghbr_F);
        cnt     = shadow_proj(nghbr_F);

        if length(t_nF) == 0
          t_G(i) = rand;
        else
          % Sort their time in increasing order
          [t_nF_sort, t_nF_ord] = sort(t_nF);
          cnt_ord     = cnt(t_nF_ord);
          
          p_base = log(1 - (1 - p_thin) .^ cnt_ord) ;     % Prob of thinned without influence of i
          p_updt = log(1 - (1 - p_thin) .^ (cnt_ord+1));  % Prob of thinned with influence of i
          p_sel  = zeros(1,length(cnt_ord)+1);             % E.g. [0 t1 t2 1]

          p_sel(1) = sum(p_updt);
          for tt = 1:(length(cnt_ord)-1)
            p_sel(tt+1) = sum(p_base(1:tt)) + sum(p_updt(tt+1:end));
          end
          p_sel(end) = sum(p_base);

          p_sel = p_sel - max(p_sel);
          p_sel = exp(p_sel); 
          p_sel = p_sel ./ sum(p_sel);
          bin = sampleDiscrete(p_sel,1);

          if(bin == 1)
            bin_st = 0;
            bin_end = t_nF_sort(bin);
          elseif(bin == length(p_sel))
            bin_st  = t_nF_sort(bin-1);
            bin_end = 1;
          else
            bin_st  = t_nF_sort(bin-1);
            bin_end = t_nF_sort(bin);
          end;
          t_G(i)      = (bin_end - bin_st) * rand + bin_st;
        end;
      else
        t_G(i)      = rand;
      end;
      tim(i,:)    = t_G(i) < t_F;
      shadow(i,:) = (dist(i,:) < R) & tim(i,:);
    end;

%       % How many observations support each of the thinned events?
%       nghbr_G     = dist_G(i,:) < R;
%       nghbr_F     = dist(i,:)   < R;

%       t_nG = t_G(nghbr_G);
%       t_nF = t_F(nghbr_F);

%       t_n  = [t_nG t_nF];
%       msk  = [ones(1,length(t_nG)) zeros(1,length(t_nF))];

%       [t_n_sort, t_n_ord] = sort(t_n);

    min_r_vec = zeros(1, num_obs);
    for i = 1:num_obs
      if num_pts == 0
        min_r = 0;
      else
        shadow_proj = sum(shadow);
        supp        = (shadow_proj == 1);
        min_r = max(supp .* (dist(i,:) .* shadow(i,:))); % Which thinned events are supported only by i (and what's the distance?)
      end;
      min_r_vec(i) = min_r;
    end;

    min_r = max(min_r_vec);
%   min_r = 1 - exp(-min_r);

%   indx = find(dist_vec > min_r,1);
%   tmp  = dist_vec(indx) - min_r;   % Probability you lie in (min_r,d_{index}+1)

%   p_vec(1:indx-1) = -inf;
%   p_vec(indx-1)   = dist_hist(indx-1) * log(1-p_thin) + log(tmp);
%   p_vec(indx:end) = dist_hist(indx:end) .* log(1-p_thin) + dist_prior(indx:end);

%   p_vec = p_vec - max(p_vec);
%   p_vec = exp(p_vec);
%   p_vec = p_vec ./ sum(p_vec);
%   bin = sampleDiscrete(p_vec,1);

%   if bin == length(dist_vec)
%     max_r = 1;
%   else
%     max_r = dist_vec(bin+1);
%   end;
%  if bin > 1
%    count = count + 1;
%  end;
%  if bin > indx-1
%    min_r = dist_vec(bin);
%    count = count + 1;
%  end;

   %R    = ((max_r - min_r) .* rand + min_r)';
   %R    = -log(1-R);

    num_thin  = sum(shadow);
    num_surv  = sum(sum(dist_G < R))/2; % We don't have to worry about times, since for each pair of overlaps there's 1 real opportunity.
    num_mh = 1;                         % Need to recalculate min_r if each mh iteration otherwise
    for(m = 1:num_mh)
      R_new   = R + 0.1*randn;
      if(R_new > min_r)
        shadow_new    = (dist < R_new) & tim;
        num_thin_new  = sum(shadow_new);
        num_surv_new  = sum(sum(dist_G < R_new))/2;
        acc = sum( log(1 - (1-p_thin).^num_thin_new) - log(1 - (1-p_thin).^num_thin)) + (num_surv_new - num_surv) * log(1-p_thin);
        if( log(rand) < acc)
          R = R_new;
          num_thin = num_thin_new;
          num_surv = num_surv_new;
          shadow   = shadow_new;
        end;
      end;
    end;

    max_r = 0;


   %%%%%%%%%%%%
   % Resample p_thin
   %%%%%%%%%%%%
   % indx = find(dist_vec > R,1)-1;
   % if(length(indx) ==0)
   %   indx = 0;
   % end;

   % p_thin = betarnd(1+num_pts, 1+indx,1,1 );

   num_mh = 5;
   for(m = 1:num_mh)
     p_thin_new = rand;
     acc = sum( log(1 - (1-p_thin_new).^num_thin) - log(1 - (1-p_thin).^num_thin)) + (num_surv) * (log(1-p_thin_new) - log(1-p_thin));
     if( log(rand) < acc)
       p_thin = p_thin_new;
     end;
   end;

   scale = gamrnd(num_obs + num_pts + sc_alp, 1/(vol + 1/sc_beta));
   poiss_mean = scale * vol;


    if(hardcore == 1)
      pred = prior_matern3(side, dim, hardcore, 1, R, scale, p_thin);
    else
      pred = prior_matern3(side, dim, hardcore, 1, R_lims, scale, p_thin);
    end;

    if(do_fano == 1)
%    mHist2d = hist2d(pred.G',vYEdge,vXEdge) ./ ((vXEdge(2) - vXEdge(1)) * (vYEdge(2) - vYEdge(1)));;
      for(ll=2:20)
        fn(ll-1) = fano(pred.G',side,ll);
      end;
    else
      fn = [];
    end;


    post(iter).G = G;
    post(iter).side = side;
    post(iter).scale = scale;
    post(iter).fano = fn; %var(mHist2d(:))/mean(mHist2d(:));
    post(iter).pred_mn = mean(mHist2d(:));
    post(iter).pred_vr = var(mHist2d(:));
    post(iter).fano_data = fano_data;
    post(iter).pred = pred.G;
    post(iter).F = F;
    post(iter).R = R;
    post(iter).t_G = t_G;
    post(iter).p_thin = p_thin;
    post(iter).num_pts = num_pts;
  end;
  count
