side = 1.6;
dim =  2;
poiss_mean = 100;
hardcore = 0;
R_lims = [0,1,.06,.02];

num_iter = 5000;
fano = zeros(num_iter,1);
for(i=1:num_iter)
  pred = prior_matern3(side, dim, hardcore, 1, R_lims, poiss_mean);
  mHist2d = hist2d(pred.G',vYEdge,vXEdge);
  fano(i) = var(mHist2d(:))/mean(mHist2d(:));
end
