function smpl = prior_np_matern3(side, dim, hardcore, stationary, R_ip, scale, p_thin, old_l, old_pts, inv_K_old_sqrt, gphyp)

  % p_thin < 1 gives probabilistic thinning

  jit          = 1e-5;
  poiss_mean   = scale * prod(side);

  if hardcore == 1
    R = R_ip;
  else
    R_max = R_ip(2);
    R_min = R_ip(1);
    R_mean = R_ip(3);
    R_std = R_ip(4);
  end;

  num_pts  = poissrnd(poiss_mean);
  E        = rand(dim,num_pts);
  for(dd = 1:dim)
    E(dd,:) = E(dd,:) * side(dd);
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if num_pts > 0
    % Instantiate GP on F points 
    ww       = old_l'*inv_K_old_sqrt;    % Whiten (Zero mean)
    Ks       = covSEiso(gphyp, old_pts', E');
    kss      = covSEiso(gphyp, E');
    pred     = Ks'*inv_K_old_sqrt;
    mn_pred  = pred*ww';
    K_pred   = kss - pred*pred';
    K_pred(1:size(K_pred)+1:end) = K_pred(1:size(K_pred)+1:end) + jit;
    [K_pred_sqrt, err] = chol(K_pred);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interpolated GP values
    w_E         = randn(size(K_pred_sqrt,1),1);
    l_E         = K_pred_sqrt' * w_E + mn_pred;
    l_E_sig     = sigmoid(l_E);

    thin_gp     = rand(size(l_E_sig)) < l_E_sig;    %Indexes points thinned by matern

    F           = E(:,thin_gp);
    num_prim_pts = size(F,2);
    t_F          = sort(rand(1,num_prim_pts));
    l_F    = l_E(thin_gp);

  else
    thin_gp         = [];
    F            = ones(dim,0);
    t_F          = [];
    l_F          = [];
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    l_G    = l_F(thin);
  else
    G      = F;
    t_G    = t_F;
    l_G    = l_F;
  end;

  smpl.num_pts = size(G,2);
  smpl.E = E;
  smpl.F = F;
  smpl.G = G;
  smpl.R = R;
  smpl.t_F = t_F;
  smpl.t_G = t_G;
  if(stationary == 0)
    smpl.l_G_sig = scale .* sigmoid(l_G);
    smpl.gphyp = gphyp;
  end;
