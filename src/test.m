
  dist_vec = sort([0;rand(10,1)]);
  p_vec = log([dist_vec(2:end);1]' - dist_vec');


  p_vec = p_vec - max(p_vec);
  p_vec = exp(p_vec);
  p_vec = p_vec ./ sum(p_vec);
for iter = 1:10000
  bin = sampleDiscrete(p_vec,1);

  if bin == length(p_vec)
    max_r = 1;
  else
    max_r = dist_vec(bin+1);
  end;
  min_r = dist_vec(bin);

  R       = ((max_r - min_r) .* rand + min_r)';
  R       = -log(1-R);
  jnk(iter) = R;

end;
