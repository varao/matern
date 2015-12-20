%wrapper_matern % Swedish
fn_pr = ([rslt(1000:end).fano]);
tmp = reshape(fn_pr,19,4001)'; 
tmp_u = quantile(tmp,.9); 
tmp_l = quantile(tmp,.1);
tmp_m = quantile(tmp,.5);
fn = ([rslt(1000).fano_data]);
plot_hist(tmp(:,5), fn(5))
