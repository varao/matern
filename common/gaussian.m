function varargout = gaussian(x, mn, cholSigma)

    % Return the log probability of x.
    varargout(1) = {mvnormpdfln(x, mn, cholSigma)};

    % Perhaps also return the derivative.
    if nargout > 1
      varargout(2) = {-solve_triu(cholSigma, ...
                                  solve_tril(cholSigma',x-mn))};

    end

