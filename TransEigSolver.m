function [varargout] = TransEigSolver(f,MMadm,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving transverse eigenproblems for modal coefficient Rj(x,y) %%%%%%%
%% An HREs solver in Cartesian coordinates with longitudinal invariance %

%% The HREs are solved for the quantity Rj satisfying the 2D Helmholtz equation
%% Rj_{xx} + Rj_{yy} + kj^2 Rj = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perfectly matching layers (PMLs) are used for truncating domain in y at
%% the endpoints MMadm.y(1) and MMadm.y(end) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The transverse modal projection methed is used to solve the transverse
%% eigenproblem for the quantity phij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phij_{yy} + kj^2 phij = kxj^2 phij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% in this function kxj is the eigenvalues sought to be solved %%%%%%%%%%

%% Author: Tengjiao He 01/01/2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Method from :  T. He, J. Liu, S. Ye, X. Qing, and S. Mo, A novel model
% order reduction technique for solving horizontal refraction equations
% in the modeling of three-dimensional underwater acoustic propagation,
% J. Sound Vib. 591: 118617 (2024)

% Any improvements for the code are welcome and bugs can be reported
% through email: hetengjiao@sjtu.edu.cn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Layered medium configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = MMadm.dpml; % thickness of PML
D = MMadm.y(end)+d; % thickness of the domain of interest
H = D+d; % total thickness including two PMLs
varargout{2} = H;
cp = 2*pi*f/MMadm.kj_m;
N = ceil(3*H*f/cp)+1; % num of transverse basis modes

[ClenC_y,wcc_y] = ClenshawCurtis(length(MMadm.y),[MMadm.y(1)+d MMadm.y(end)+d]);  % water column discretion
ClenC_k = interp1(MMadm.y+d,MMadm.kj,ClenC_y,'spline');
[y_pml1,wcc_pml1] = ClenshawCurtis (MMadm.npml,[0 d]);
[y_pml2,wcc_pml2] = ClenshawCurtis (MMadm.npml,[D H]);

% wavenumbers of PMLs
k_pml1 = MMadm.kj(1);
k_pml2 = MMadm.kj(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set-up of PMLs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% positive damping coefficients
p = 4;
q = 4;

% exponential damping function
tao1 = -(d-y_pml1)/d;
s1 = exp(p*tao1)-1i*(exp(q*tao1)-1);
tao2 = -(y_pml2-D)/d;
s2 = exp(p*tao2)-1i*(exp(q*tao2)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multimodal basis psi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% modal projection phi ---> psi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sine basis
m = (1:N).';
n = m.';
psim = sqrt(2/H)*sin((m)*pi/H.*ClenC_y);
psim_pml1 = sqrt(2/H)*sin((m)*pi/H.*y_pml1);
psim_pml1_p = m.*pi/H.*sqrt(2/H).*cos((m)*pi/H.*y_pml1);
psim_pml2 = sqrt(2/H)*sin((m)*pi/H.*y_pml2);
psim_pml2_p = m.*pi/H.*sqrt(2/H).*cos((m)*pi/H.*y_pml2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coefficient matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% d^2_psi/dy^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = m+n; b = m-n; upper = H-d; lower = d;

J_PML1 = -0.5*d*(wcc_pml1.'.*(1./s1).*psim_pml1_p)*psim_pml1_p.';
J_PML2 = -0.5*d*(wcc_pml2.'.*(1./s2).*psim_pml2_p)*psim_pml2_p.';

% sine basis
J1 = upper.*sinc(a/H*upper)-lower.*sinc(a/H*lower);
J2 = upper.*sinc(b/H*upper)-lower.*sinc(b/H*lower);
J  = -pi^2/H^3*(m.*n).*(J1+J2);

J = J+J_PML1+J_PML2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% kxj^2 Psi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_PML1 = 0.5*d*(wcc_pml1.'.*s1.*psim_pml1)*psim_pml1.';
L_PML2 = 0.5*d*(wcc_pml2.'.*s2.*psim_pml2)*psim_pml2.';

% sine basis
L = -1/H*(J1-J2);

L = (L+L_PML1+L_PML2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% krj(x,y)^2 Psi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_PML1 = k_pml1^2*L_PML1;
K_PML2 = k_pml2^2*L_PML2;
K = 0.5*(MMadm.y(end)-MMadm.y(1))*(wcc_y.'.*ClenC_k.^2.*psim)*psim.';
K = K+K_PML1+K_PML2;

C = L\(J+K);
varargout{1} = C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sort eigenvalues and reconstruct eigenfunctions %%%%%%%%%%%%%%%%%%%%%%%
%%% inverse modal projection psi ---> phi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if id == 1

    [a_mn,kx_unsort] = eig(C);
    [~,index] = sort(real(diag(sqrt(kx_unsort))),'descend');
    kx_unsort = sqrt(diag(kx_unsort));

    % transverse eigenvalues kxj
    varargout{3} = real(kx_unsort(index))+1i*abs(imag(kx_unsort(index)));

    A = a_mn(:,index);
    varargout{4} = A;

    % reconstruction of transverse eigenfunction phij
    phin = psim.'*A;
    phim = phin.';

    phin_pml1 = psim_pml1.'*A;
    phim_pml1 = phin_pml1.';
    phin_pml2 = psim_pml2.'*A;
    phim_pml2 = phin_pml2.';

    O = 0.5*(MMadm.y(end)-MMadm.y(1))*(wcc_y.'.*phim)*phin;

    O_PML1 = 0.5*d*(wcc_pml1.'.*s1.*phim_pml1)*phin_pml1;
    O_PML2 = 0.5*d*(wcc_pml2.'.*s2.*phim_pml2)*phin_pml2;
    O = O + O_PML1 + O_PML2;
    varargout{5} = diag(O); % modal normalized coefficient
    varargout{6} = index;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end