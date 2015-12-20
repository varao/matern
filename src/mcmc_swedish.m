load ../data/swedish_trees;

G = swedish_trees';

side = 5;

norm = 1 * max(max(abs(G)));
G = side .* G ./ norm;

%post = matern3_prob(G,side,1000);
%post = matern3_hc(G,side,10000);
post = matern3_np(G,side,100);
