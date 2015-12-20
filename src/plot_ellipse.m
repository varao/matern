function plot_ellipse(pts, R, hardcore, col)

  plot(pts(1,:),pts(2,:),'b.')
  hold on;
  num_pts = size(pts,2);

  if hardcore == 1
    for(i = 1:num_pts)
      ellipse(median(R), median(R), 0, pts(1,i), pts(2,i), col);
      ellipse(max(R), max(R), 0, pts(1,i), pts(2,i), col);
    end;
  else
    num_iter = length(R)/num_pts;
    R = reshape(R, num_pts, num_iter);
    for(i = 1:num_pts)
      R_c = R(i,:);
      ellipse(median(R_c), median(R_c), 0, pts(1,i), pts(2,i), col);
%      ellipse(quantile(R_c,.75), quantile(R_c,.75), 0, pts(1,i), pts(2,i), 'b');
%     ellipse(quantile(R_c,.05), quantile(R_c,.05), 0, pts(1,i), pts(2,i), 'm');
    end;
  end;
  axis square
