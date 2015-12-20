function post = matern3_np(G, side, num_iter)
% MCMC sampler for a nonstationary Matern type-III repulsive process with a 
% logistic GP-prior on the intensity
% Assumes a common radius of repulsion for all points with a non-informative 
% prior on it

  %%%%%%%%%%%%
  % Hyperpriors assume data is scaled so that order 1 point lies per unit area
  %%%%%%%%%%%%

  % Gamma prior on scale
  sc_alp  = 1;
  sc_beta = 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Log-normal prior on GP hyperparameters
  gphyp        = log([.1;.1]);
  gphyp_cvchol = eye(2);
  length_sc_pr = @(theta) gaussian(theta, gphyp, gphyp_cvchol); 

  aux          = 1;
  jit          = 1e-5;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dim     = size(G,1);
  num_obs = size(G,2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initialize RVs
  scale = 1;
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

  % Initialize R randomly
  R = min(min(dist_G))/5;
  %%%%%%%%%%%%%%%
     
  vol        = side^dim;
  poiss_mean = scale * vol;

  for iter = 1:num_iter

    if(rem(iter,100) == 0)
      iter
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample GP-thinned events
    num_pts  = poissrnd(poiss_mean);
    E_new   = rand(dim,num_pts) .* side;

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
    num_prim_pts  = poissrnd(poiss_mean);
    F_new         = rand(dim,num_prim_pts) .* side;

    if num_prim_pts > 0

      t_F          = rand(1,num_prim_pts);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      % F can include only those in the shadow
      dist   = sqrt(sq_dist(G,F_new));
      tim    = t_G(ones(num_prim_pts,1),:)' < t_F(ones(num_obs,1),:);
      shadow = (dist < R) & tim;

      thin   = any(shadow);

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

    num_prim_pts      = length(l_F);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample times of observations
    for i = 1:num_obs
      if num_prim_pts > 0
        shadow_proj = sum(shadow); % How many observations support each of the thinned events?
        supp        = (shadow(i,shadow_proj == 1));

        if(any(supp))
           max_time = min(tim(supp));
        else
           max_time = 1;
        end;
        t_G(i)      = max_time * rand;
        tim(i,:)    = t_G(i) < t_F;
        shadow(i,:) = (dist(i,:) < R) & tim(i,:);
      else
        t_G(i)      = rand;
      end;
    end;

    scale = gamrnd(num_obs + num_pts + sc_alp, 1/(vol + 1/sc_beta));
    poiss_mean = scale * vol;

    if num_prim_pts == 0
      min_r = 0;
    else
      min_r = max(max(dist(shadow)));
      if length(min_r) == 0
        min_r = 0;
      end;
    end;
    max_r = min(min(dist_G));

    % Non-informative prior on R
    R = (max_r - min_r) * rand + min_r;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot grid
  gr = 1;
  if gr == 1 & dim == 2
    t_tmp = [0:side/10:side];
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
    post(iter).num_pts = num_pts;
    post(iter).E = E;
    post(iter).F = F;
    post(iter).R = R;
    post(iter).t_G = t_G;
    post(iter).max_R = max_r;
    post(iter).min_R = min_r;
    post(iter).l_E_sig = scale .* sigmoid(l_E);
    post(iter).l_F_sig = scale .* sigmoid(l_F);
    post(iter).l_G_sig = scale .* sigmoid(l_G);
    post(iter).gphyp = gphyp;
  end;
