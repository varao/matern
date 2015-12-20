function plot_hist(samples, val)

  clf;
  f = (min(samples));
  c = (max(samples));
  %f = 0; c = .25;
  st = (c-f)/15;
  ed = [f:st:c];

  h = histc(samples,ed);
  h = h ./ sum(h);

  bar(ed, h,'b');
  axis square tight ;
  xlim([f-0.005,c]);
  ym = max(h);
  hold on;
  if(nargin == 2)
    m = floor((val - f)/st);
    plot(m*st + f, h(m+1)+.01, 'or', 'MarkerSize',6, 'MarkerFaceColor', 'red' );
  end
  ylim([0,ym+.05]);
  xlabel('Fano factor');
  ylabel('Posterior predictive distribution');
  xlabel('Thinning radius');
  ylabel('Posterior probability');
  set(gcf,'Paperposition',[1 1 3 2])

  %print -depsc plot_postpred_fano.eps
  %!epstopdf plot_postpred_fano.eps
  print -depsc plot.eps
  !epstopdf plot.eps
