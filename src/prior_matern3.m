function smpl = prior_matern3(side, dim, hardcore, stationary, R_ip, scale, p_thin)

  % p_thin < 1 gives probabilistic thinning

  jit          = 1e-5;
  poiss_mean   = scale * prod(side);

  num_pts  = poissrnd(poiss_mean);
  E        = rand(dim,num_pts);
  for(dd = 1:dim)
    E(dd,:) = E(dd,:) * side(dd);
  end;

  if hardcore == 1
    R = R_ip;
  else
    R_max = R_ip(2);
    R_min = R_ip(1);
    R_mean = R_ip(3);
    R_std = R_ip(4);
  end;

  if stationary == 1
    F    = E;
    l_E  = [];
  else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log-normal prior on GP hyperparameters
    gphyp        = log([.1;.1]);
    gphyp_cvchol = eye(2);
    length_sc_pr = @(theta) gaussian(theta, gphyp, gphyp_cvchol); 

    % GP covariance evaluated at current GF values
    K_E          = covSEiso_jit(gphyp, E', jit);
    K_E_sqrt     = chol(K_E);

    w_E         = randn(size(K_E_sqrt,1),1);
    l_E         = K_E_sqrt' * w_E;
    l_E_sig     = sigmoid(l_E);

    thin_gp     = rand(size(l_E_sig)) > l_E_sig;    %Indexes points thinned by GP
    F           = E(:,thin_gp);
  end;

  num_prim_pts = size(F,2);
  t_F          = sort(rand(1,num_prim_pts));

  if hardcore == 0
    R = rand(1,num_prim_pts) .* (R_max - R_min) + R_min;
 %  R = R_mean + R_std.*randn(1,num_prim_pts);
  end;

  if(num_prim_pts > 0)
    thin   = ones(1,num_prim_pts);
    dist   = sqrt(sq_dist(F,F));
    dist(1:num_prim_pts+1:end) = inf;

    if hardcore == 1
      shadow = (dist < R);
    else
      shadow = (dist < R(ones(num_prim_pts,1),:)');
    end;

    for(ii = 2:num_prim_pts)
      if(p_thin > 0.999999)
        if(any(shadow(1:ii-1,ii)))
          shadow(ii,:) = 0;
          thin(ii) = 0;
        end  
      else
        surv = sum(shadow(1:ii-1,ii));
        if(log(rand) > surv*log(1-p_thin))   % prob(surv) = (1-p_thin)^num
            shadow(ii,:) = 0;
            thin(ii) = 0;
        end  
      end;
    end;
    thin = logical(thin);

    G      = F(:,thin);
    t_G    = t_F(thin);
  else
    G      = F;
    t_G    = t_F;
  end;

  smpl.num_pts = size(G,2);
  smpl.E = E;
  smpl.F = F;
  smpl.G = G;
  smpl.R = R;
  smpl.t_F = t_F;
  smpl.t_G = t_G;
  if(stationary == 0)
    smpl.l_E_sig = scale .* sigmoid(l_E);
    smpl.gphyp = gphyp;
  end;
