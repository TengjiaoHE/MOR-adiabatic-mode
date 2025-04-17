%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An order-reduced solver for HREs in the modeling of 3D underwater %%%%
%% acoustic propagationg in longitudinally invariant environments   %%%%%
%% A code using the transverse modal projection methed %%%%%%%%%%%%%%%%%%
%% Horizontal domain for the vertical mode's coefficients truncated %%%%%
%% by two PMLs                                          %%%%%%%%%%%%%%%%%
%% Author: Tengjiao He 01/01/2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T. He, J. Liu, S. Ye, X. Qing, and S. Mo, A novel model order reduction 
% technique for solving horizontal refraction equations in the modeling of 
% three-dimensional underwater acoustic propagation, J. Sound Vib. 
% 591: 118617 (2024)

% Any improvements for the code are welcome and bugs can be reported
% through email: hetengjiao@sjtu.edu.cn


clear variables;
close all;

clc;

f = 25; % Frequency in [Hz]
omega = 2*pi*f;
c = 1500; % Reference sound speed
k0 = 2*pi*f/c;
lambda = c/f;

% x- and y- grids initialization
ymax = 8000; % Maximal of transverse range
ny = 1024; % Num of discretization points for transverse range
MORparam.y = linspace(0,ymax,ny); % Transverse range grids
MORparam.y0 = 4000; % Source position on the y-axis

xmax = 25000; % Maximal of propagation range
nx = 1024; % Num of discretization points for transverse range
MORparam.x = linspace(0,xmax,nx); % Propagation range grids

% wavenumbers for various values of water depth
kh = dlmread('ASA_wedge_kj/kj_wedge_att.txt');

% mode eigenfunctions at z = z_s = 100 m for various values of water depth
phizsh = dlmread('ASA_wedge_kj/phizs_wedge.txt');

% mode eigenfunctions at z = z_r = 30 m for various values of water depth
phizrh = dlmread('ASA_wedge_kj/phizr_wedge.txt');

mmax = 8; % Num of vertical modes

% geometry of seafloor (a coastal wedge)
h0 = 200;
talph = 200 / 4000;
MORparam.hy = h0 +  (MORparam.y - MORparam.y0) * talph;

% parameters of PMLs
MORparam.dpml = 5*lambda; % Thickness of PML
MORparam.npml = 100; % Num of discretization points for PML
    
p = 0;

h_fig = figure;
set(h_fig,'position',[400 100 1200 300])

% main loop for different vertical modes
for ii = 1:mmax

    disp(['Solving modal coefficients: ' int2str(ii)  ' of ' int2str(mmax)  ' modes'] );
    
    MORparam.kj = interp1(kh(:,1),kh(:,ii+1),MORparam.hy);
    ih1 = find(MORparam.hy>=kh(end,1),1,'first');
    ih2 = find(MORparam.hy<=kh(1,1),1,'last');
    iy0 = find(MORparam.y>=MORparam.y0,1,'first');
    MORparam.kj(1:ih1) = MORparam.kj(ih1);
    MORparam.kj(ih2:end) = MORparam.kj(ih2);
    MORparam.kj_m = max(max(real(MORparam.kj)));
    
    % excitation of vertical modes
    vmode(1:length(MORparam.y),1) = interp1(phizsh(:,1),phizsh(:,ii+1),MORparam.hy);
    vmode_zs = vmode(iy0);
    
    % solving transverse eigenproblems for transverse eigenfunctions Phi
    [~,H,kxj,A,Nm] = TransEigSolver(f,MORparam,1);


    ys = MORparam.y0+MORparam.dpml; % y-axis of the source
    yr = linspace(MORparam.dpml,ymax+MORparam.dpml,ny); % y-axis of the receiver
    hyr = interp1(MORparam.y,MORparam.hy,yr-MORparam.dpml);

    % vertical mode shape
    vmoder(1:length(hyr),1) = interp1(phizrh(:,1),phizrh(:,ii+1),hyr); 
    ih1 = find(hyr>=kh(end,1),1,'first');
    ih2 = find(hyr<=kh(1,1),1,'last');
    vmoder(1:ih1) = vmoder(ih1); 
    vmoder(ih2:end) = vmoder(ih2);
    
    % reconstruction of Rj(x,y) using modal projection phi ---> psi
    n = 1:length(kxj);
    psi_n = sqrt(2/H)*sin((n)*pi/H.*yr.'); % Transverse basis modes psi_n
    phi_m = (psi_n)*(A); % Transverse eigenfunctons phi_m
    Am = sqrt(2/H)*sin((n)*pi/H*ys)*(A)./Nm.';
    Prop = exp(1i*kxj*MORparam.x)./(kxj)*1i/2;
    Rj = 4*pi*(phi_m.*Am*Prop); % Coefficients of vertical modes
    
    % summation of vertical modes
    p = p + vmode_zs * Rj .* repmat(vmoder,1,nx);
    
    % plotting Rj(x,y)
    if ii <= 5

        subplot(1,5,ii)
        pcolor(yr/1000,MORparam.x/1000,mag2db(abs(vmode_zs *Rj).'));
        view(2);
        shading flat;
        caxis( [-70, -30] );
        h_Colorbar = colorbar;
        h_ColorbarTitle = title(h_Colorbar,'dB');
        xlabel('y (km)')
        ylabel( 'x (km)' ) ;
        title( { [ 'Mode' num2str( ii )] } )
        set(gca,'FontSize',12)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off',  ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'ZColor', [.3 .3 .3],'LineWidth', 1)

    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% results plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[all_themes, all_colors] = GetColors();
color1 = all_themes{26};
color2 = all_themes{14};
color3 = all_themes{21};
color = [color1(1,:);color2(1,:);color3(3,:);0.06,0.51,0.04;1,1,0.07;1,0,0;0.7,0,0];
map = GenColormap(color, 512);


tlt = -20*log10(abs( p ));

h_fig = figure;
set(h_fig,'position',[500 250 900 350])
pcolor(MORparam.x/1000,MORparam.y/1000,tlt);
view(2);
set(gca,'ydir','reverse') ;
colormap(flipud(map))
shading flat;
caxis( [ 50, 80 ] );
h_Colorbar = colorbar;
h_ColorbarTitle = title(h_Colorbar,'dB');
xlabel('x (km)')
ylabel( 'y (km)' ) ;
title( { deblank( 'Modal projection for HRE' ); [ 'Freq = ' num2str( f ) ' Hz' ] } )
set(gca,'FontSize',12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off',  ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'ZColor', [.3 .3 .3],'LineWidth', 1)

iy0 = find(yr>=ys,1,'first');
TLtrack = tlt(iy0,:);
figure;
hold all;
xx = dlmread('ASA_wedge_kj/cross_wedge_ASA.TL');
plot(xx(:,1),-xx(:,2),'r','LineWidth',2);
plot(MORparam.x/1000,TLtrack,'k.','MarkerSize',3);
xlabel('x, km');
ylabel('TL, dB');
set(gca,'FontSize',12)
xlim([0 MORparam.x(end)/1000])
ylim([40 100])
set( gca, 'YDir', 'Reverse' )
legend('Ref','Modal projection')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', ...
    'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off',  ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'ZColor', [.3 .3 .3],'LineWidth', 1)



