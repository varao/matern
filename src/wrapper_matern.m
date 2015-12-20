side = 10

%load redwood;
%G = 1.6 .* redwood';

 load swedish_trees;
 swedish_trees(:,1) = swedish_trees(:,1) ./ max(swedish_trees(:,1));
 swedish_trees(:,2) = swedish_trees(:,2) ./ max(swedish_trees(:,2));
G = swedish_trees';
% G = .16 .* swedish_trees';

%load amacrine_off
%load amacrine_on
%amacrine_on(2,:) = amacrine_on(2,:) .* 1.6;
%G = 10 .* amacrine_on;
%load poiss_data
%rslt = matern_retina(amacrine_on, amacrine_off, 1.6,5000);

load nerve_mild1; nerve = nerve_mild1;
load nerve_mod1; nerve = nerve_mod1;
% load nerve_normal; nerve = nerve_normal;
G = nerve';

 G(1,:) = G(1,:) ./ max(G(1,:));
 G(2,:) = G(2,:) ./ max(G(2,:));
 G = G .* side;

%rslt = matern3_prob(G, side,5000);
%rslt = matern3_hc(G, side,5000);
%rslt = matern3_prob(G, [8,5.5],5000);
%rslt = matern3_prob(G, [8,5.5],5000);

rslt = matern3_prob(G, side,500);
