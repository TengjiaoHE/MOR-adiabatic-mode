function [kh,phizrh,phizsh,phizzsh,rhoh] = SweepMode(zs,zr,z,f,mmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sweeping waterborne modes over different receiver depths %%%%%%%%%%%%%
%% A normal code using the transverse modal projection methed %%%%%%%%%%%
%% Author: Tengjiao He 01/01/2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T. He, X. Liu, R. Nie, W. Guo, et al. Semi-analytical solution for sound 
% propagation from a moving directional source in a shallow-water waveguide
% , J. Sound Vib. 576(35):118259 (2024).


% Any improvements for the code are welcome and bugs can be reported
% through email: hetengjiao@sjtu.edu.cn


% wavenumbers for various values of water depth
kh(:,1) = z.';

% mode eigenfunctions at z = z_s for various values of water depth
phizsh(:,1) = z.';
phizzsh(:,1) = z.';

% mode eigenfunctions at z = z_r for various values of water depth
phizrh(:,1) = z.';

% density at z = z_r for various values of water depth
rhoh(:,1) = z.';

% loop for 
for ii = 1:length(z)

    Param.h = z(ii);

    % waveguide parameters
    
    param;

    % solving eigenvalues and modes

    [kr,phizs,phizzs,phizr,rho] = ModeSolver(f,Param,zs,zr);

    rhoh(ii,2:length(zr)+1) = rho;
    kh(ii,(1:mmax)+1) = kr(1:mmax);
    phizsh(ii,(1:mmax)+1) = phizs(:,1:mmax);
    phizzsh(ii,(1:mmax)+1) = phizzs(:,1:mmax);

    for jj = 1:length(zr)

        phizrh(ii,(1:mmax)+1,jj) = phizr(jj,1:mmax);

    end

end