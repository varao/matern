function [F, t_F, dist, tim, shadow] = get_thinned(G, t_G, R, poiss_mean, dim, side, hardcore)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample thinned events
    num_obs  = size(G,2);
    num_pts  = poissrnd(poiss_mean);

    if num_pts > 0
      F      = rand(dim,num_pts);
      for(dd = 1:dim)
        F(dd,:) = F(dd,:) * side(dd);
      end;
      t_F    = rand(1,num_pts);
      dist   = sqrt(sq_dist(G,F));
      tim    = t_G(ones(num_pts,1),:)' < t_F(ones(num_obs,1),:);  % Which G's support which F's (time only)

      if hardcore == 1
        shadow = (dist < R) & tim;
      else
        shadow = (dist < R(ones(num_pts,1),:)') & tim;  % Which G's support which F's
      end;

      thin   = any(shadow);  % Which F's are supported by any G's

      F      = F(:,thin);
      t_F    = t_F(:,thin);
      dist   = dist(:,thin);
      tim    = tim (:,thin);
      shadow = shadow(:,thin);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      F      = ones(dim,0);
      t_F    = [];
      dist   = zeros(num_pts,0);
      tim    = zeros (num_pts,0);
      shadow = zeros(num_pts,0);
    end;
