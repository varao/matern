function plot_post()

  clf;
  load LARGE_STAT/stationary_large;
  lu = max([[post_5.scale],[post_10.scale],[post.scale]]);
  ll = min([[post_5.scale],[post_10.scale],[post.scale]]);
  lim = [ll:(lu-ll)/10:lu];

  subplot(1,3,1)
  jnk = histc([post_5.scale],lim);
  jnk = jnk ./ sum(jnk);
  bar(lim,jnk)
  xlim([ll,lu])
  axis square tight ;
  ym = max(jnk);

  jnk = histc([post_10.scale],lim);
  jnk = jnk ./ sum(jnk);
  subplot(1,3,2)
  bar(lim,jnk)
  xlim([ll,lu])
  axis square tight ;
  ym = max([ym,jnk]);

  jnk = histc([post.scale],lim);
  jnk = jnk ./ sum(jnk);
  subplot(1,3,3)
  bar(lim,jnk)
  xlim([ll,lu])
  axis square tight ;
  ym = max([ym,jnk]);

  ym = ym .* 1.1

  subplot(1,3,1)
  ylim([0,ym])
  subplot(1,3,2)
  ylim([0,ym])
  subplot(1,3,3)
  ylim([0,ym])
