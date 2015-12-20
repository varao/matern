function post = matern3_hc(G, side, num_iter)
% MCMC sampler for a stationary Matern type-III repulsive process
% If hardcore = 1
%    Assumes a common radius of repulsion for all points with a non-informative 
%    prior on it
%
% else
%    each point has its own radius

  %%%%%%%%%%%%%%%%%%

  hardcore = 1;

  % Gamma prior on scale (i.e. Poisson rate)
  sc_alp  = 10;
  sc_beta = 10;

  scale = 1;
  %Hyperpriors on distribution of Matern radii. Both are non-informative

  R_lower_sc         = -1;  % Truncated power on lower limit of R interval
  R_length_sc        = -1;   % Pareto prior on length of R interval
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dim     = size(G,1);
  num_obs = size(G,2);

  if(length(side) == 1)
    side = side .* ones(dim,1);
  end;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vXEdge = linspace(0,side(1),5);
  vYEdge = linspace(0,side(2),5);
  mHist2d = hist2d(G',vYEdge,vXEdge) ./ ((vXEdge(2) - vXEdge(1)) * (vYEdge(2) - vYEdge(1)));


  do_fano = 1
  if(do_fano == 1)    % Fano
    for(ll=2:20)
      fano_data(ll-1) = fano(G',side,ll);
    end;
  else             % Ripley
    fano_data = [];
  end;

  %%%%%%%%%%%%%%%

  dist_G   = sqrt(sq_dist(G,G));
  dist_G(1:num_obs+1:end) = inf;
  t_G      = rand(1,num_obs);

  if hardcore == 1
    num_R = 1;
  else
    num_R = num_obs;
    rad_hyp1 = 2; rad_hyp2 = 1;
  end;

  % Initialize R's randomly
  R = ones(1,num_R) .* min(min(dist_G))/2;

  R_lower = min(R)/2
  R_upper = max(R)*2;
  %%%%%%%%%%%%%%%
     
  vol        = prod(side);
  poiss_mean = scale * vol;
  R_lims = [0,2,0,1];


  for iter = 1:num_iter

    if(rem(iter,500) == 0)
      iter
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample thinned events
    [F, t_F, dist_F, tim, shadow]       = get_thinned(G, t_G, R,  poiss_mean, dim, side, hardcore);
    [t_G, shadow] = get_times(G, t_G, R, shadow, t_F, tim, dist_F, hardcore);
    R  = get_rad(G, t_G, R, shadow, t_G, tim, dist_G, dist_F, R_lims, hardcore);
    scale = gamrnd(num_obs + length(t_F) + sc_alp, 1/(vol + 1/sc_beta));

    mx  = max(R);
    avg = mean(R);
    nm = length(t_F) + length(t_G) ;

    k_inv = 1;
%    R_upper = gprnd(1/(k_inv + nm), mx/(k_inv+nm), mx);

%    R_lower = randpow(.1, min(R));
    if(hardcore == 0)
      R_lower = randpow(.1, min(R));
      R_length = randp(1, 1, num_R + R_length_sc, max(R - R_lower));
      R_upper = R_lower + R_length;
    end;
   
    R_mean = avg; 
    R_std  = .1;

    m0 = 1; n0 = 10; t0 = 1;
    R_var = gamrnd(0.1+nm, 1 + 0.5.*sum((R-avg).*(R-avg)) + nm*(avg-1)*(avg-1)/(2*(nm+n0)));
    R_std = sqrt(R_var);


    if(isnan(R_std)) error 'WTF?'
    end;
    R_mean = n0*m0/(nm+n0) + nm*avg/(nm+n0) + sqrt((n0+nm)*t0)*randn();
    R_lims = [R_lower, R_upper, R_mean, R_std];
%    R_lims = [0, 2, 1, .2];

    poiss_mean = scale * vol;

    if(hardcore == 1)
      pred = prior_matern3(side, dim, hardcore, 1, R, scale, 1);
    else
      pred = prior_matern3(side, dim, hardcore, 1, R_lims, scale, 1);
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
     post(iter).R_mean = R_mean;
     post(iter).R_std = R_std;
    post(iter).t_G = t_G;
    post(iter).num_thin = nm - length(t_G);
  end;
