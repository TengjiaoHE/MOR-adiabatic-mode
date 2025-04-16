function [kr,Phis_m,Phizs_m,Phi_m,rho] = ModeSolver(f,Mode,zs,zr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving waterborne modes in stratified fluid media %%%%%%%%%%%%%%%%%%%
%% A normal code using the transverse modal projection methed %%%%%%%%%%%
%% Author: Tengjiao He 01/01/2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Method from : T. He, X. Liu, R. Nie, W. Guo, et al. Semi-analytical 
% solution for sound propagation from a moving directional source in a 
% shallow-water waveguide, J. Sound Vib. 576(35):118259 (2024).


% Any improvements for the code are welcome and bugs can be reported
% through email: hetengjiao@sjtu.edu.cn



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Layered medium configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 2*pi*f;
D = Mode.z{end}(end); % total depth of the waveguide
H = D;
N = 3*H*f/1500; % num of basis modes

for i = 1:Mode.n

    [ClenCurtis_z,wcc_z] = ClenshawCurtis(length(Mode.z{i}),[Mode.z{i}(1) Mode.z{i}(end)]);  % water column discretion
    ClenCurtis.z{i} = ClenCurtis_z;
    ClenCurtis.wcc{i} = wcc_z;
    ClenCurtis.rho{i} = interp1(Mode.z{i},Mode.rho{i},ClenCurtis.z{i},'spline');
    ClenCurtis.alpha{i} = interp1(Mode.z{i},Mode.alpha{i},ClenCurtis.z{i},'spline');
    ClenCurtis.c{i} = interp1(Mode.z{i},Mode.c{i},ClenCurtis.z{i},'spline')./(1+1i*(ClenCurtis.alpha{i})/54.575);
    ClenCurtis.k{i} = omega./ClenCurtis.c{i};

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multimodal basis function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = (0:N-1).'+1;
n = m.';
for i = 1:Mode.n

    Psi_m{i} = sqrt(2/H)*sin((m)*pi/H.*ClenCurtis.z{i});
    Psi_n{i} = Psi_m{i}.';
    Psi_n_prime{i} = (n).*pi/H.*sqrt(2/H).*cos((n)*pi/H.*ClenCurtis.z{i}.');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coefficient matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% d^2_Psi/dz^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = 0;
for i = 1:Mode.n
    J_it = -0.5*(Mode.z{i}(end)-Mode.z{i}(1))*(1./ClenCurtis.rho{i}.*ClenCurtis.wcc{i}.'.*Psi_n_prime{i}.')*Psi_n_prime{i};
    J = J + J_it;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% lambda^2 Psi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 0;
for i = 1:Mode.n
    L_it = 0.5*(Mode.z{i}(end)-Mode.z{i}(1))*(1./ClenCurtis.rho{i}.*ClenCurtis.wcc{i}.'.*Psi_m{i})*Psi_n{i};
    L = L + L_it;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% w^2/c^2(z) Psi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 0;
for i = 1:Mode.n
    K_it = 0.5*(Mode.z{i}(end)-Mode.z{i}(1))*(1./ClenCurtis.rho{i}.*ClenCurtis.wcc{i}.'.*ClenCurtis.k{i}.^2.*Psi_m{i})*Psi_n{i};
    K = K + K_it;
end

C = L\(J+K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Eigenvalue sort and eigenfunction reconstruction %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a_mn,kr] = eig(C);
cp = omega./real(sqrt(diag(kr)));
[~,id] = sort(cp,"ascend");

% eigenvalues

kr = sqrt(diag(kr));
kr = kr(id);
kr = real(kr)+1i*abs(imag(kr));

% eigenfunction reconstruction

A = a_mn(:,id);
for i = 1:Mode.n
    Phi_n{i} = Psi_n{i}*A;
    Phi_m{i} = Phi_n{i}.';
end

O = 0;
for i = 1:Mode.n
    O_it = 0.5*(Mode.z{i}(end)-Mode.z{i}(1))*(1./ClenCurtis.rho{i}.*ClenCurtis.wcc{i}.'.*Phi_m{i})*Phi_n{i};
    O = O + O_it;
end

% modal normalized coefficient

Nm = diag(sqrtm(O)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Psi_n = sqrt(2/H)*sin(n*pi/H.*zr.');
Phi_m = Psi_n*A./Nm.';

Phis_m = sqrt(2/H)*sin(n*pi/H*zs)*(A)./Nm.';
Phizs_m = n*pi/H*sqrt(2/H).*cos(n*pi/H*zs)*(A)./Nm.';


z = [Mode.z{1}(1:end-1) Mode.z{2}];
rho0 = [Mode.rho{1}(1:end-1) Mode.rho{2}];

rho = interp1(z,rho0,zr,"nearest");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end