function [t_G, shadow] = get_times(G, t_G, R, shadow, t_F, tim, dist, hardcore)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Sample times of observations

  [num_obs, num_pts] = size(shadow);

  if num_pts > 0
    % How many observations support each of the thinned events?
    shadow_proj = sum(shadow);
    % Which thinned events are supported only by one observation?
    supp        = (shadow_proj == 1);
  end;

  for i = 1:num_obs
    shad_curr = shadow(i,:);
    if num_pts > 0
      supp_curr = supp & shad_curr; % Which thinned events are only supported by i?
      if(any(supp_curr))
         max_time = min(tim(i,supp_curr));
      else
         max_time = 1;
      end;
      t_G(i)      = max_time * rand;
      tim(i,:)    = t_G(i) < t_F;
    else
      t_G(i)      = rand;
    end;

    if hardcore == 1
      shadow(i,:) = (dist(i,:) < R) & tim(i,:);
    else
      shadow(i,:) = (dist(i,:) < R(i)) & tim(i,:);
    end;

  end;

