function t = trunc_exp(t1, t2, A)
  % Sample from a truncated exponential (restricted to the interval [t1,t2] with rate parameter A

  A = -A;
  t = logdiffexp_v(A*t1, log(rand) + A*t2 + log(exp(A*(t1-t2))-1))/A;

  if(not(isreal(t)))
    error 'WTF?';
  end;
