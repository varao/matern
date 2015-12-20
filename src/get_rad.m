function R = get_rad(G, t_G, R, shadow, t_surv, tim, dist_G, dist_F, R_lims, hardcore)

  [num_obs, num_pts] = size(shadow);
  R_lower = R_lims(1); R_upper = R_lims(2);
  R_mean  = R_lims(3); R_std = R_lims(4);

  if num_pts > 0
    % How many observations support each of the thinned events?
    shadow_proj = sum(shadow);
    % Which thinned events are supported by only one observation?
    supp        = (shadow_proj == 1);
  end;

  if(hardcore == 1)
    if num_pts == 0
      min_r = 0;
    else
      min_r = max(min(dist_F .* tim));  % The max of the nearest support point of all thinned events.
    end;
    max_r   = min(min(dist_G, [], 2));

    min_r   = max(min_r, R_lower);
    max_r   = min(max_r, R_upper);

    R       = ((max_r - min_r) .* rand + min_r)';
    %R = R_mean + R_std.*randn(1,num_prim_pts);
  else
    for i = 1:num_obs
      shad_curr = shadow(i,:);
      if num_pts == 0
        min_r = 0;
      else
        supp_curr = supp & shad_curr;    % Which thinned events are only supported by i?
%        min_r = max(supp .* (dist_F(i,:) .* shadow(i,:)));  
        if(any(supp_curr))
           min_r    = max(dist_F(i,supp_curr));  
        else
           min_r = 0;
        end;
      end;

      min_r   = max(min_r, R_lower);

      max_r   = min(dist_G(i,t_G(i) < t_surv));
      max_r   = min([max_r, R_upper]);

      
      R(i)       = ((max_r - min_r) .* rand + min_r)';
%     R(i)       =  trunc_norm(R_mean,R_std,min_r,max_r);

      shadow(i,:) = (dist_F(i,:) < R(i)) & tim(i,:);
      shadow_proj = sum(shadow);
      supp        = (shadow_proj == 1);
    end;
  end;
