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

f      = 60; % Frequency in [Hz]
omega  = 2*pi*f;
c      = 1500; % Reference sound speed
lambda = c/f;
nz     = 1;
zr     = 45; % Receiver depth in [m]
zs     = 12; % Source depth in [m]


% number of modes taken into account

mmax = 6;

% solving horizontal wavenumbers and local modes

z = 12:0.5:50; % depth sweep grid

[kh,phizrh,phizsh,phizzsh] = SweepMode(zs,zr,z,f,mmax);

% x- and y- grids initialization

ymax     = 4000; % Maximal of transverse range
ny       = floor(ymax/lambda*5); % Num of discretization points for transverse range
MORparam.y  = linspace(0,ymax,ny); % Transverse range grids
MORparam.y0 = 2000; % Source position on the y-axis

xmax     = 15000; % Maximal of propagation range
nx       = 1024; % Num of discretization points for transverse range
MORparam.x  = linspace(0,xmax,nx); % Propagation range grids

% geometry of internalwaves (an ISW train)

h0 = 15;
cISW = 0.65;
dy = 0*60*cISW:10*cISW:0*60*cISW+60*60*cISW;
MORparam.hy = h0;
for j = 1:12

    MORparam.hy =  MORparam.hy+10 *sech((MORparam.y-1700-25 - (j-1)*400+dy(1))/70).^2;

end

% parameters of PMLs

MORparam.dpml = 5*lambda; % Thickness of PML
MORparam.npml = 100; % Num of discretization points for PML


p = 0; A = zeros(ny,nx,nz); Rj = zeros(ny,nx,2);

h_fig = figure;
set(h_fig,'position',[400 100 1200 300])

% main loop for different vertical modes

for ii = 1:mmax

    disp(['Solving modal coefficients: ' int2str(ii)  ' of ' int2str(mmax)  ' modes'] );

    MORparam.kj = interp1(kh(:,1),kh(:,ii+1),MORparam.hy);
    MORparam.kj_m = max(max(real(MORparam.kj)));

    % indexing the position of the source

    iy0 = find(MORparam.y>=MORparam.y0,1,'first');


    % mode function

    vmode = interp1(phizsh(:,1),phizsh(:,ii+1),MORparam.hy);
    vmodez = interp1(phizzsh(:,1),phizzsh(:,ii+1),MORparam.hy);
    vmode_zs = vmode(iy0);

    % mode amplitude

    [~,H,kxj,U,Nm] = TransEigSolver(f,MORparam,1);

    ys = MORparam.y0+MORparam.dpml; % source depth
    yr = linspace(MORparam.dpml,ymax+MORparam.dpml,ny); % receiver depth
    hyr = interp1(MORparam.y,MORparam.hy,yr-MORparam.dpml);

    n = 1:length(kxj);


    % reconstruction of Rj(x,y) using modal projection phi ---> psi
    Psi_n = sqrt(2/H)*sin(n*pi/H.*yr.');
    Phi_m = Psi_n*U./Nm.';
    Prop = exp(1i*kxj*MORparam.x)./kxj*1i/2;


    Am = sqrt(2/H)*sin(n*pi/H*ys)*U;
    Rj = 4*pi*(Phi_m.*Am*Prop);


    % sweeping the depth

    for jj = 1:length(zr)

        wmode_zr = interp1(phizrh(:,1),squeeze(phizrh(:,ii+1,jj)),hyr);
        ih1 = find(hyr>=kh(1,1),1,'first');
        ih2 = find(hyr<=kh(end,1),1,'last');
        wmode_zr(1:ih1) = wmode_zr(ih1);
        wmode_zr(ih2:end) = wmode_zr(ih2);

        A(:,:,jj) = vmode_zs*Rj .* repmat(wmode_zr,nx,1).';

    end

    p = p + A;

    if ii <= 5

        subplot(1,5,ii)
        pcolor(yr/1000,MORparam.x/1000,mag2db(abs(Rj).'));
        view(2);
        shading flat;
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
tlt = abs( p );
tlt( isnan( tlt ) ) = NaN;           % remove NaNs
tlt( isinf( tlt ) ) = 1e-16;          % remove infinities
icount = find( tlt > 1e-17 );         % for stats, only these values count
tlt( tlt < 1e-17 ) = 1e-17;            % remove zeros
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log
tlmed = median( tlt(icount) );       % median value
tlstd = std( tlt(icount) );          % standard deviation
tlmax = tlmed + 0.75 * tlstd;        % max for colorbar
tlmax = 10 * round( tlmax / 10 );    % make sure the limits are round numbers
tlmin = tlmax - 50;                  % min for colorbar


x_interp = linspace(min(MORparam.x/1000),max(MORparam.x/1000),1024);
y_interp = linspace(min(MORparam.y/1000),max(MORparam.y/1000),512);


h_fig = figure;
set(h_fig,'position',[400 100 1000 400])
izr = find(zr>=zs,1,'first');
[X,Y,Z] = meshgrid(x_interp,y_interp,zr);
P = interp2(MORparam.x/1000,MORparam.y/1000.',squeeze(tlt(:,:,izr)),x_interp,y_interp.');
surf(X(:,:,izr),Y(:,:,izr),Z(:,:,izr),P);
set(gca,'ydir','reverse') ;
view(2);
shading interp
colormap(flipud(Jet))
caxis([ tlmin, tlmax ] )
axis('tight')
hold on
xlabel('x / km')
ylabel('y / km')
box on
set(gca,'FontSize',18)
hold on
