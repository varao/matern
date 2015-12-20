%load redwood_hc_rslt;
%load swedish_hc_rslt;

burnin = 1000;
  jnk = [];
  for i = burnin:length(rslt)
    jnk = [jnk, rslt(i).num_pts];
  end

[rw, x] = hist(jnk,10);
bar(x,rw)
rw = rw ./ sum(rw);
%ix = x < 2.5;
%bar(x(ix),rw(ix))
bar(x,rw)
%xlim([0.3, 0.5]);
%xlim([0.5,0.95]);

set(gcf,'Paperposition',[1 1 3 2]);
print -depsc plot_r.eps
!epstopdf plot_r.eps
