function D = DirFun(phi0,p_f,phi)
    D_real = (interp1(phi0,real(p_f),phi,'spline'));
    D_imag = (interp1(phi0,imag(p_f),phi,'spline'));
    D = abs(D_real+D_imag*1i);
end