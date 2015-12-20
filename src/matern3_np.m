function post = matern3_np(G, side, num_iter)
% MCMC sampler for a nonstationary Matern type-III repulsive process with a 
% logistic GP-prior on the intensity
% Assumes a common radius of repulsion for all points with a non-informative 
% prior on it

  vXEdge = linspace(0,side(1),side);
  vYEdge = linspace(0,side(2),side);
  mHist2d = hist2d(G',vYEdge,vXEdge) ./ ((vXEdge(2) - vXEdge(1)) * (vYEdge(2) - vYEdge(1)));


  do_fano = 1
  if(do_fano == 1)    % Fano
    for(ll=2:20)
      fano_data(ll-1) = fano(G',side,ll);
    end;
  else             % Ripley
    fano_data = [];
  end;

  %%%%%%%%%%%%%%%%%%

  hardcore = 1;

  %%%%%%%%%%%%
  % Hyperpriors assume data is scaled so that order 1 point lies per unit area
  %%%%%%%%%%%%

  % Gamma prior on scale
  sc_alp  = 10;
  sc_beta = 10;

  scale = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Log-normal prior on GP hyperparameters
  gphyp        = log([.1;.1]);
  gphyp_cvchol = eye(2);
  length_sc_pr = @(theta) gaussian(theta, gphyp, gphyp_cvchol); 

  aux          = 1;
  jit          = 1e-5;


  %Hyperpriors on distribution of Matern radii. Both are non-informative

  R_lower_sc         = -1;  % Truncated power on lower limit of R interval
  R_length_sc        = -1;  % Pareto prior on length of R interval
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dim     = size(G,1);
  num_obs = size(G,2);

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

  %%%%%%%%%%%%%%%

  dist_G   = sqrt(sq_dist(G,G));
  dist_G(1:num_obs+1:end) = inf;
  t_G      = rand(1,num_obs);

  if hardcore == 1
    num_R = 1;
  else
    num_R = num_obs;
  end;

  % Initialize R's randomly
  R = ones(1,num_R) .* min(min(dist_G))/5;

  R_lower = min(R)/5;
  R_upper = max(R)*2;
  %%%%%%%%%%%%%%%
     
  vol        = prod(side);
  poiss_mean = scale * vol;
  R_lims = [0,2,0,1];

  for iter = 1:num_iter

    if(rem(iter,50) == 0)
      iter
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample Matern-thinned events
    [F_new, t_F, dist_F, tim, shadow]       = get_thinned(G, t_G, R,  poiss_mean, dim, side, hardcore);
    num_prim_pts = length(t_F);

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
      dist_F   = dist_F(:,thin);
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
    num_prim_pts      = length(l_F);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t_G, shadow] = get_times(G, t_G, R, shadow, t_F, tim, dist_F, hardcore);
    R  = get_rad(G, t_G, R, shadow, t_G, tim, dist_G, dist_F, R_lims, hardcore);
    scale = gamrnd(num_obs + num_pts + sc_alp, 1/(vol + 1/sc_beta));
    poiss_mean = scale * vol;

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

    if(hardcore == 1)
      pred = prior_np_matern3(side, dim, hardcore, 0, R, scale, 0, xx, [E,F,G], inv_K_EFG_sqrt, gphyp);
    else
      pred = prior_np_matern3(side, dim, hardcore, 0, R_lims, scale, 0, xx, [E,F,G], inv_K_EFG_sqrt, gphyp);
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
  if gr == 1 & dim == 2
    t_tmp1 = [0:side(1)/10:side(1)];
    t_tmp2 = [0:side(2)/10:side(2)];
    [t1,t2] = meshgrid(t_tmp, t_tmp);
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


    post(iter).scale = scale;
    post(iter).fano = fn; %var(mHist2d(:))/mean(mHist2d(:));
    post(iter).fano_data = fano_data;
    post(iter).pred = pred.G;
    post(iter).pred_l = pred.l_G_sig;
    post(iter).num_pts = num_pts;
    post(iter).E = E;
    post(iter).F = F;
    post(iter).R = R;
    post(iter).t_G = t_G;
    post(iter).l_E_sig = scale .* sigmoid(l_E);
    post(iter).l_F_sig = scale .* sigmoid(l_F);
    post(iter).l_G_sig = scale .* sigmoid(l_G);
    post(iter).gphyp = gphyp;
  end;
