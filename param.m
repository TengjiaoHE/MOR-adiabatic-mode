%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter definition for solving vertical modes %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho     = [1000,1000,1000,1800]; % density of each layer
c       = [1530,1515,1480,1750]; % sound speed of each layer
alpha   = [0,0,0,0.4375]; % wave attenuation in [dB/lambda] for each layer
lam0    = min(c)/f;
h       = [Param.h,Param.h+15,70,500]; % depth of each layer
N_h     = ceil(h/lam0*10); % Num of discretization points for each layer


Param.n = length(h); % Num of layers

for i = 1:Param.n

    if i == 1

        Param.z{i} = linspace(0,h(i),N_h(i));

    else

        Param.z{i} = linspace(h(i-1),h(i),N_h(i));

    end

    Param.c{i} = ones(1,length(Param.z{i}))*c(i);
    Param.rho{i} = ones(1,length(Param.z{i}))*rho(i);
    Param.alpha{i} = ones(1,length(Param.z{i}))*alpha(i);

end

dz = Param.z{end} -Param.z{end}(1);

% attenuation parameters for the ABL
a1 = 1;
a2 = 54.57/(h(end)-h(end-1))-0.0001;
a3 = Param.c{end}(1)/(h(end)-h(end-1))-0.0001;

Param.alpha{end} =  Param.alpha{end}(1)./sqrt(a1-a2/54.57*dz);
% Param.c{end} =  Param.c{end}(1)./sqrt(a1-a3/Param.c{end}(1)*dz);

% Gaussian attenuation
lam2 = Param.c{end}(1)/f;
Param.alpha{end} = fliplr(Param.alpha{end}(1)+100*exp(-(dz).^2/(2*lam2^2)));

% ssp in the thermocline layer where the internal wave occurs
Param.c{2} = -(c(1)-c(3))/(h(2)-h(1))*(Param.z{2}-h(1))+c(1);

% figure
% plot([Param.c{1},Param.c{2},Param.c{3}],[Param.z{1},Param.z{2},Param.z{3}])
% 
% figure
% plot([Param.alpha{end}],[Param.z{end}])