function K_GH = covSEiso_jit(gphyp, GH,jit)
    K_GH = covSEiso(gphyp, GH);
    K_GH(1:size(K_GH)+1:end) = K_GH(1:size(K_GH)+1:end) + jit;
