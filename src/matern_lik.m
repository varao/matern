function lik = matern_lik(x, num_pts, num_prim_pts)

    x_sig    = sigmoid(x);

    l_E_sig   = x_sig(1:num_pts);
    l_F_sig   = x_sig(num_pts+1:num_pts+num_prim_pts);
    l_G_sig   = x_sig(num_pts+num_prim_pts+1:end);

    lik      = sum(log(l_F_sig)) + sum(log(l_G_sig)) + sum(log(1-l_E_sig));
