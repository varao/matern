function post = matern3_prob_np(G, side, num_iter)
% MCMC sampler for a nonstationary Matern type-III repulsive process with probabilistic thinning
% If hardcore = 1
%    Assumes a common radius of repulsion for all points with a non-informative 
%    prior on it
% else
%    each point has its own radius
%
% Within each radius, the thinning probability is a Beta-distributed p_thin

  hardcore = 1;

  % Gamma prior on scale (i.e. Poisson rate)
  sc_alp  = 1;
  sc_beta = 1;

  scale = 1;


  %%%%%%%%%%%%
  % Hyperpriors assume data is scaled so that order 1 point lies per unit area
  %%%%%%%%%%%%

  % Gamma prior on scale
  sc_alp  = 10;
  sc_beta = 10;

  scale = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Log-normal prior on GP hyperparameters
  gphyp        = log([.5;.75]);
  gphyp_cvchol = eye(2);
  length_sc_pr = @(theta) gaussian(theta, gphyp, gphyp_cvchol); 

  aux          = 1;
  jit          = 1e-5;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p_thin = .8;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dim     = size(G,1);
  num_obs = size(G,2);

  if(length(side) == 1)
    side = side .* ones(dim,1);
  end;

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initialize RVs
  E     = zeros(dim,0);
  F     = zeros(dim,0);
  l_E   = [];
  l_F   = [];
  l_G   = zeros(1,num_obs)';
  % GP covariance evaluated at current GF values
  K_EFG          = covSEiso_jit(gphyp, [E,F,G]', jit);
  K_EFG_sqrt     = chol(K_EFG);
  inv_K_EFG_sqrt = inv(K_EFG_sqrt);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  R= 0.0005;

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample GP-thinned events
    num_pts  = poissrnd(poiss_mean);
    E_new   = rand(dim,num_pts);
    for(dd = 1:dim)
      E_new(dd,:) = E_new(dd,:) * side(dd);
    end;

    if num_pts > 0

      % Instantiate GP on E points
      ww       = [l_E; l_F; l_G]'*inv_K_EFG_sqrt;    % Whiten (Zero mean)
      Ks       = covSEiso(gphyp, [E,F,G]', E_new');
      kss      = covSEiso(gphyp, E_new');
      pred     = Ks'*inv_K_EFG_sqrt;
      mn_pred  = pred*ww';
      K_pred   = kss - pred*pred';
      K_pred(1:size(K_pred)+1:end) = K_pred(1:size(K_pred)+1:end) + jit;
      [K_pred_sqrt, err] = chol(K_pred);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Interpolated GP values
      w_E         = randn(size(K_pred_sqrt,1),1);
      l_E         = K_pred_sqrt' * w_E + mn_pred;
      l_E_sig     = sigmoid(l_E);

      thin_gp     = rand(size(l_E_sig)) > l_E_sig;    %Indexes points thinned by GP
      E           = E_new(:,thin_gp);
      l_E         = l_E(thin_gp);

    else
      thin_gp      = [];
      E            = ones(dim,0);
      l_E          = [];
    end;

    K_EFG          = (covSEiso_jit(gphyp, [E,F,G]', jit));
    K_EFG_sqrt     = chol(K_EFG);
    inv_K_EFG_sqrt = inv(K_EFG_sqrt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample Matern thinned events
    num_prim_pts  = poissrnd(poiss_mean);
    F_new         = rand(dim,num_prim_pts);
    for(dd = 1:dim)
      F_new(dd,:) = F_new(dd,:) * side(dd);
    end;

    if num_prim_pts > 0

      t_F          = rand(1,num_prim_pts);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      % F can include only those in the shadow
      dist   = sqrt(sq_dist(G,F_new));
      tim    = t_G(ones(num_prim_pts,1),:)' < t_F(ones(num_obs,1),:);
      shadow = (dist < R) & tim;

      % Some parts of the shadow have higher thinning probability than others
      thin_logp   = sum(shadow) * log(1 - p_thin);  % One minus Shadow, actually (prob_thin = 1 - (1-p_thin)^num)
      thin        = thin_logp < log(rand(1, num_prim_pts));

      F_new  = F_new(:,thin);
      t_F    = t_F(thin);
      dist   = dist(:,thin);
      tim    = tim (:,thin);
      shadow = shadow(:,thin);
      num_prim_pts = length(t_F);
    end;

    if num_prim_pts > 0
      % Instantiate GP on F points (exchange order with shadowing for efficiency)
      ww       = [l_E;l_F; l_G]'*inv_K_EFG_sqrt;    % Whiten (Zero mean)
      Ks       = covSEiso(gphyp, [E,F,G]', F_new');
      kss      = covSEiso(gphyp, F_new');
      pred     = Ks'*inv_K_EFG_sqrt;
      mn_pred  = pred*ww';
      K_pred   = kss - pred*pred';
      K_pred(1:size(K_pred)+1:end) = K_pred(1:size(K_pred)+1:end) + jit;
      [K_pred_sqrt, err] = chol(K_pred);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Interpolated GP values
      w_F         = randn(size(K_pred_sqrt,1),1);
      l_F         = K_pred_sqrt' * w_F + mn_pred;
      l_F_sig     = sigmoid(l_F);

      thin        = rand(size(l_F_sig)) < l_F_sig;    %Indexes points thinned by matern

      F           = F_new(:,thin);
      t_F    = t_F(thin);
      l_F    = l_F(thin);
      dist   = dist(:,thin);
      tim    = tim (:,thin);
      shadow = shadow(:,thin);

    else
      thin         = [];
      F            = ones(dim,0);
      t_F          = [];
      l_F          = [];
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_pts      = length(l_E);
    num_prim_pts = length(l_F);

    K_EFG          = (covSEiso_jit(gphyp, [E,F,G]', jit));
    K_EFG_sqrt     = chol(K_EFG);
    inv_K_EFG_sqrt = inv(K_EFG_sqrt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Now resample the entire GP

    l_pr = randn(1, num_pts + num_prim_pts + num_obs) * K_EFG_sqrt;

    curr_lik_func = @(x) matern_lik(x, num_pts, num_prim_pts);

    curr_log_lik   = curr_lik_func([l_E;l_F; l_G]);

    [xx, phi, curr_log_lik] = elliptical_slice([l_E; l_F; l_G], l_pr, @(x) curr_lik_func(x), curr_log_lik, 0);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GP hyperparameters
    % Assuming 0-mean

     [gphyp, xx, aux, K_EFG_sqrt] = update_theta_aux_surr(gphyp, xx, @(x) curr_lik_func(x), ...
                                                                 @(x) covSEiso_jit(x, [E, F, G]',jit), aux, length_sc_pr, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    l_E   = xx(1:num_pts);
    l_F   = xx(num_pts+1:num_pts+num_prim_pts);
    l_G   = xx(num_pts+num_prim_pts+1:end);

    inv_K_EFG_sqrt = inv(K_EFG_sqrt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample times of observations

    for i = 1:num_obs
      if num_prim_pts > 0
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

    min_r_vec = zeros(1, num_obs);
    for i = 1:num_obs
      if num_prim_pts == 0
        min_r = 0;
      else
        shadow_proj = sum(shadow);
        supp        = (shadow_proj == 1);
        min_r = max(supp .* (dist(i,:) .* shadow(i,:))); % Which thinned events are supported only by i (and what's the distance?)
      end;
      min_r_vec(i) = min_r;
    end;

    min_r = max(min_r_vec);

    num_thin  = sum(shadow);
    num_surv  = sum(sum(dist_G < R))/2; % We don't have to worry about times, since for each pair of overlaps there's 1 real opportunity.
    num_mh = 1;                         % Need to recalculate min_r if each mh iteration otherwise
    for(m = 1:num_mh)
      R_new   = R + 0.02*randn;
%      R_new= 0.0005;
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

   num_mh = 1;
   for(m = 1:num_mh)
     p_thin_new = rand;
     acc = sum( log(1 - (1-p_thin_new).^num_thin) - log(1 - (1-p_thin).^num_thin)) + (num_surv) * (log(1-p_thin_new) - log(1-p_thin));
     if( log(rand) < acc)
       p_thin = p_thin_new;
     end;
   end;
   p_thin = 0.75;
   scale = gamrnd(num_obs + num_prim_pts + sc_alp, 1/(vol + 1/sc_beta));
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


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot grid
  gr = 1;
  grd_sz = 25;
  if gr == 1 & dim == 2
    t_tmp1 = [-.5:side(1)/grd_sz:(side(1)+1.0)];
%    t_tmp1 = [-1.75:side(1)/grd_sz:(side(1)+.75)];  NC
    t_tmp2 = [-0.5:side(2)/grd_sz:(side(2)+.5)];
    [t1,t2] = meshgrid(t_tmp1, t_tmp2);
    t1 = t1(:)'; t2 = t2(:)';
    Kt  = covSEiso(gphyp, [E,F,G]', [t1;t2]');
    ktt = covSEiso(gphyp, [t1;t2]');

    pred_t     = Kt'*inv_K_EFG_sqrt;
    mn_pred_t  = pred_t*inv_K_EFG_sqrt'*[l_E;l_F;l_G];
    K_pred_t   = ktt - pred_t*pred_t';
    K_pred_t(1:size(K_pred_t)+1:end) = K_pred_t(1:size(K_pred_t)+1:end) + jit;

    K_pred_t_sqrt = chol(K_pred_t);
    l_t = mn_pred_t + K_pred_t_sqrt' * randn(size(K_pred_t_sqrt,1),1);
    l_t_sig = scale .* sigmoid(l_t);
    grid.l = l_t_sig;
    grid.x = t1;
    grid.y = t2;

    post(iter).grid = grid;
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    post(iter).G = G;
    post(iter).side = side;
    post(iter).scale = scale;
    post(iter).fano = fn; %var(mHist2d(:))/mean(mHist2d(:));
    post(iter).fano_data = fano_data;
    post(iter).pred = pred.G;
    post(iter).F = F;
    post(iter).R = R;
    post(iter).t_G = t_G;
    post(iter).p_thin = p_thin;
    post(iter).num_prim_pts = num_prim_pts;
    post(iter).l_E_sig = scale .* sigmoid(l_E);
    post(iter).l_F_sig = scale .* sigmoid(l_F);
    post(iter).l_G_sig = scale .* sigmoid(l_G);
    post(iter).gphyp = gphyp;
  end;
  count
