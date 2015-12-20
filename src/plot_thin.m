clf
side = 5;
subplot(1,2,1)
histo2D([post.F]',[0,side],30,[0,side],30,'x','y','t');
title('Repulsion')
hold on; plot(G(1,:),G(2,:),'.');

hold on; plot(G(1,:),G(2,:),'.');
subplot(1,2,2)
histo2D([post.E]',[0,side],30,[0,side],30,'x','y','t');
hold on; plot(G(1,:),G(2,:),'.');
title('Thinning')

set(gcf,'paperposition',[1 1 6 3]);
print -depsc redwood_rep.eps
