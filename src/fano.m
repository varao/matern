function fn = fano(data, side, num)

  vXEdge = linspace(0,side(1),num);
  vYEdge = linspace(0,side(2),num);
  mHist2d = hist2d(data,vYEdge,vXEdge);
%  [var(mHist2d(:)) mean(mHist2d(:))]
  fn = var(mHist2d(:))/mean(mHist2d(:));
