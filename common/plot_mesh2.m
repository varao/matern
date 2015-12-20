function plot_mesh(grd,G)

burnin = 300;

num_samples = length(grd);

gr = grd(1).l;
for i = 1:num_samples
  gr = [gr,grd(i).l];
end;

gr = gr(:,burnin+1:end);
%gr = gr(:,end-400);

gr_mn = (mean(gr,2));
gr_mn = (reshape(gr_mn,11,11));
gr_vr = (std(gr,0,2));
gr_vr = (reshape(gr_vr,11,11));

x     = reshape(grd(1).x,11,11);
y     = reshape(grd(1).y,11,11);

clf;
subplot(1,2,1);
surf(x,y,reshape(gr_mn,11,11));  % Or contour
contour(x,y,reshape(gr_mn,11,11));  % Or contour
title('Mean')
colorbar;
hold on;
plot(G(1,:),G(2,:),'.');
shading flat ; axis square tight ; grid on ; 
set(gca,'XTick',[])
set(gca,'YTick',[])


subplot(1,2,2);
surf(x,y,reshape(gr_vr,11,11));
contour(x,y,reshape(gr_vr,11,11));
title('Std dev.')
colorbar;
hold on;
plot(G(1,:),G(2,:),'.');
min(min(gr_mn))
max(max(gr_mn))
shading flat ; axis square tight ; grid on ; 
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'paperposition',[1 1 6 3]);
print -depsc redwood_mn_vr.eps
