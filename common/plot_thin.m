clf
subplot(1,2,1)
%histo2D([state.F]',[0,5],10,[0,5],10,'x','y','t');
%histo2D([state.F]',[0,5],10,[0,5],10,'x','y','t');
histo2D([state.F]',[0,5],10,[0,5],10,'x','y','t');
title('Repulsion')
hold on; plot(G(1,:),G(2,:),'.');
subplot(1,2,2)
histo2D([state.E]',[0,5],10,[0,5],10,'x','y','t');
hold on; plot(G(1,:),G(2,:),'.');
set(gcf,'paperposition',[1 1 6 3]);
print -depsc redwood_thin.eps
