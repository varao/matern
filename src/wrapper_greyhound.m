load GreyhoundSouth.csv
G = GreyhoundSouth';

G(1,:) = (G(1,:) - min(G(1,:))); % ./ (max(G(1,:)) - min(G(1,:)));
G(2,:) = (G(2,:) - min(G(2,:))); % ./ (max(G(2,:)) - min(G(2,:)));
side = max(G,[],2);

%rslt = matern3_prob(G, side, 5000);
%rslt = matern3_hc(G, side, 5000);
%rslt = matern3_prob(G, side, 5000);

rslt = matern3_prob_np(G, side, 5000);
