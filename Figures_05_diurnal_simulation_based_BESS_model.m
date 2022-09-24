%% enviromental parameters obtain from CB field measure during growth season
% sampling
load('D:\variables.mat')
Ta = [variables(3:end,1); variables(1:2,1)];
PAR = variables(:,2)/4;
VPD = [variables(3:end,3); variables(1:2,3)]/10*0.5; %Kpa
vapor_satu = 0.611.*exp(17.502.*Ta./(Ta +240.97));
RH=(vapor_satu-VPD)./vapor_satu*100*0.8; %0-100;
figure;plot(Ta);

pmin=[
    1,    % PAR [W m-2]; uniform
    1,    % air temperature [oC]; normal
    1,   % surface pressure KPa; normal
    1,    % wind speed [m s-1]; normal
    1,   % relative humidity [0-100]; normal
    3,  % leaf area index [0-10]     ; uniform
    0.3,   % clumping index [0.3-0.9]; uniform
    5,   % canopy height [m]       ; uniform
    40, % Vcmax25                    ; uniform
    3,   % g1                     ; uniform
    0.7, % PAR absorptance;          ; normal
    0.5, % NIR_absorptance
    0.94]; % emissivity;  uniform

pmax=[
    1,   % PAR [W m-2]          ; uniform
    1,    % air temperature [oC] ; normal
    1, % surface pressure KPa ; normal
    1,    % wind speed [m s-1]   ; normal
    1,    % relative humidity [0-100]; normal
    7,   % leaf area index [0-10]      ; unifrom
    0.9,    % clumping index [0.3-0.9] ; uniform
    40,    % canopy height [m]         ; uniform
    100, % Vcmax25                     ; uniform
    6,   % g1                          ; uniform
    0.99, % PAR absorptance            ; normal
    0.9,  % NIR_absorptance            ; normal
    0.99]; % emissivity;               ; uniform

Np = length(pmin);
N = 10000*ones(1,1);
skip=randi(1000);leap=randi([1000 10000]);
p = sobolset(Np+1,'Skip',skip,'Leap',leap);
p = scramble(p,'MatousekAffineOwen');
X0 = net(p,N);
X1=X0(:,2:end);
% PAR = (X1(:,1).*(pmax(1)-pmin(1)))+pmin(1);
PARs = repmat(PAR,[1,N]);
Tas = repmat(Ta,[1,N]);
Ps  = repmat(87*ones(24,1),[1,N]); %figure;histogram(Ps, 50, 'Normalization','PDF')
u   = repmat(2*ones(24,1),[1,N]);
RHs = repmat(RH,[1,N]);
LAI = repmat(X1(:,1).*(pmax(6)-pmin(6))+pmin(6),[1,24])';
% CI = repmat(X1(:,1).*(pmax(7)-pmin(7))+pmin(7),[1,24])';
CI  = repmat(0.5*ones(24,1),[1,N]);
hc = repmat(X1(:,1).*(pmax(8)-pmin(8))+pmin(8),[1,24])';
Vcmax25 = repmat(X1(:,1).*(pmax(9)-pmin(9))+pmin(9),[1,24])';
g1 = repmat(X1(:,1).*(pmax(10)-pmin(10))+pmin(10),[1,24])';
PAR_absorptance = repmat(X1(:,1).*(pmax(11)-pmin(11))+pmin(11),[1,24])';
NIR_absorptance = repmat(X1(:,1).*(pmax(12)-pmin(12))+pmin(12),[1,24])';
emissivity = repmat(X1(:,1).*(pmax(13)-pmin(13))+pmin(13),[1,24])';
ca = repmat(400*ones(24,1),[1,N]); %ambient CO2 ppm

out0 = BESS_revised(PAR,Ta,Ps,u,RH,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca);
Tc = mean(out0.TcSu,2)-273.15;
Tc_std = 2*std(out0.TcSu,0,2);

LAI = repmat(0.01*ones(24,1),[1,N]);
out1 = BESS_revised(PAR,Ta,Ps,u,RH,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca);
Ts = mean(out1.Ts,2)-273.15;

figure;
plot(smooth(Ta,5),'LineWidth',1.5);hold on;
plot(smooth([Ts(2:end); Ts(1)],5),'LineWidth',1.5);
% plot(smooth([Tc(2:end); Tc(1)],5),'r-','LineWidth',2)
errorbar(1:24,smooth([Tc(2:end); Tc(1)],5),Tc_std,'r-','LineWidth',1)
ylim([7 27]);
set(gcf,'position',[500,500,650*0.5,400*0.5])
print(gcf, 'D:\Data\Global Thermoregulation\For RSE\Figures\figure_05_forest.jpg', '-djpeg', '-r600');


%% US-KS2 shrubland site
% sampling
load('variables_USKS2.mat')
variables=variables_USKS2;
% Ta = [variables(3:end,1); variables(1:2,1)];
Ta = variables(1:2:end,1);
PAR = variables(1:2:end,2)/4;
VPD = variables(1:2:end,3)/10*0.5; %Kpa
vapor_satu = 0.611.*exp(17.502.*Ta./(Ta +240.97));
RH=(vapor_satu-VPD)./vapor_satu*100*0.85 ;%0-100;
figure;plot(Ta);

pmin=[
    1,    % PAR [W m-2];uniform
    1,    % air temperature [oC]; normal
    1,   % surface pressure KPa; normal
    1,    % wind speed [m s-1]; normal
    1,   % relative humidity [0-100]; normal
    3,  % leaf area index [0-10]     ; uniform
    0.3,   % clumping index [0.3-0.9]; uniform
    2,   % canopy height [m]       ; uniform
    40, % Vcmax25                    ; uniform
    3,   % g1                     ; uniform
    0.7, % PAR absorptance;          ; normal
    0.5, % NIR_absorptance
    0.94]; % emissivity;  uniform

pmax=[
    1,   % PAR [W m-2]          ; uniform
    1,    % air temperature [oC] ; normal
    1, % surface pressure KPa ; normal
    1,    % wind speed [m s-1]   ; normal
    1,    % relative humidity [0-100]; normal
    7,   % leaf area index [0-10]      ; unifrom
    0.9,    % clumping index [0.3-0.9] ; uniform
    8,    % canopy height [m]         ; uniform
    80, % Vcmax25                     ; uniform
    5,   % g1                      ; uniform
    0.99, % PAR absorptance            ; normal
    0.9,  % NIR_absorptance            ; normal
    0.99]; % emissivity;               ; uniform

Np = length(pmin);
N = 10000*ones(1,1);
skip=randi(1000);leap=randi([1000 10000]);
p = sobolset(Np+1,'Skip',skip,'Leap',leap);
p = scramble(p,'MatousekAffineOwen');
X0 = net(p,N);
X1=X0(:,2:end);
% PAR = (X1(:,1).*(pmax(1)-pmin(1)))+pmin(1);
PARs = repmat(PAR,[1,N]);
Tas = repmat(Ta,[1,N]);
Ps  = repmat(87*ones(24,1),[1,N]); %figure;histogram(Ps, 50, 'Normalization','PDF')
u   = repmat(3*ones(24,1),[1,N]);
RHs = repmat(RH,[1,N]);
LAI = repmat(X1(:,1).*(pmax(6)-pmin(6))+pmin(6),[1,24])';
% CI = repmat(X1(:,1).*(pmax(7)-pmin(7))+pmin(7),[1,24])';
CI  = repmat(0.5*ones(24,1),[1,N]);
hc = repmat(X1(:,1).*(pmax(8)-pmin(8))+pmin(8),[1,24])';
Vcmax25 = repmat(X1(:,1).*(pmax(9)-pmin(9))+pmin(9),[1,24])';
g1 = repmat(X1(:,1).*(pmax(10)-pmin(10))+pmin(10),[1,24])';
PAR_absorptance = repmat(X1(:,1).*(pmax(11)-pmin(11))+pmin(11),[1,24])';
NIR_absorptance = repmat(X1(:,1).*(pmax(12)-pmin(12))+pmin(12),[1,24])';
emissivity = repmat(X1(:,1).*(pmax(13)-pmin(13))+pmin(13),[1,24])';

ca = repmat(400*ones(24,1),[1,N]); %ambient CO2 ppm

% X2 =[PARs,Ta,Ps,u,RHs,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca];
out0 = BESS_revised(PARs,Tas,Ps,u,RHs,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca);
% out0 = BESS_revised(X2(:,1),X2(:,2),X2(:,3),X2(:,4),X2(:,5),X2(:,6),X2(:,7),X2(:,8),X2(:,9),X2(:,10),X2(:,11),X2(:,12),X2(:,13),X2(:,14));
Tc = mean(out0.TcSu,2)-273.15;
Tc_std = 1*std(out0.TcSu,0,2);

LAI = repmat(0.01*ones(24,1),[1,N]);
out1 = BESS_revised(PARs,Tas,Ps,u,RHs,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca);
Ts = mean(out1.Ts,2)-273.15;
% Ts = Ts+(Ts-mean(Tc,2))*1.5;

figure;
plot(smooth(Ta,5),'LineWidth',1.5);hold on;
plot(smooth([Ts(2:end); Ts(1)],5),'LineWidth',1.5);
% plot(smooth([Tc(2:end); Tc(1)],5),'r-','LineWidth',2)
errorbar(1:24,smooth([Tc(2:end); Tc(1)],5),Tc_std,'r-','LineWidth',1)
ylim([22 38]);
set(gcf,'position',[500,500,650*0.5,400*0.5])
print(gcf, 'D:\Data\Global Thermoregulation\For RSE\Figures\figure_05_shrubland.jpg', '-djpeg', '-r600');


%%
% sampling
load('variables_CHFRU.mat')
variables=variables_CHFRU;
% Ta = [variables(3:end,1); variables(1:2,1)];
Ta = variables(1:2:end,1);
PAR = variables(1:2:end,2)/4;
VPD = variables(1:2:end,3)/10*0.5; %Kpa
vapor_satu = 0.611.*exp(17.502.*Ta./(Ta +240.97));
RH=(vapor_satu-VPD)./vapor_satu*100*0.85 ;%0-100;
figure;plot(Ta);

pmin=[
    1,    % PAR [W m-2];uniform
    1,    % air temperature [oC]; normal
    1,   % surface pressure KPa; normal
    1,    % wind speed [m s-1]; normal
    1,   % relative humidity [0-100]; normal
    3,  % leaf area index [0-10]     ; uniform
    0.3,   % clumping index [0.3-0.9]; uniform
    2,   % canopy height [m]       ; uniform
    40, % Vcmax25                    ; uniform
    3,   % g1                     ; uniform
    0.7, % PAR absorptance;          ; normal
    0.5, % NIR_absorptance
    0.94]; % emissivity;  uniform

pmax=[
    1,   % PAR [W m-2]          ; uniform
    1,    % air temperature [oC] ; normal
    1, % surface pressure KPa ; normal
    1,    % wind speed [m s-1]   ; normal
    1,    % relative humidity [0-100]; normal
    7,   % leaf area index [0-10]      ; unifrom
    0.9,    % clumping index [0.3-0.9] ; uniform
    8,    % canopy height [m]         ; uniform
    80, % Vcmax25                     ; uniform
    5,   % g1                      ; uniform
    0.99, % PAR absorptance            ; normal
    0.9,  % NIR_absorptance            ; normal
    0.99]; % emissivity;               ; uniform

Np = length(pmin);
N = 10000*ones(1,1);
skip=randi(1000);leap=randi([1000 10000]);
p = sobolset(Np+1,'Skip',skip,'Leap',leap);
p = scramble(p,'MatousekAffineOwen');
X0 = net(p,N);
X1=X0(:,2:end);
% PAR = (X1(:,1).*(pmax(1)-pmin(1)))+pmin(1);
PARs = repmat(PAR,[1,N]);
Tas = repmat(Ta,[1,N]);
Ps  = repmat(87*ones(24,1),[1,N]); %figure;histogram(Ps, 50, 'Normalization','PDF')
u   = repmat(2*ones(24,1),[1,N]);
RHs = repmat(RH,[1,N]);
LAI = repmat(X1(:,1).*(pmax(6)-pmin(6))+pmin(6),[1,24])';
% CI = repmat(X1(:,1).*(pmax(7)-pmin(7))+pmin(7),[1,24])';
CI  = repmat(0.5*ones(24,1),[1,N]);
hc = repmat(X1(:,1).*(pmax(8)-pmin(8))+pmin(8),[1,24])';
Vcmax25 = repmat(X1(:,1).*(pmax(9)-pmin(9))+pmin(9),[1,24])';
g1 = repmat(X1(:,1).*(pmax(10)-pmin(10))+pmin(10),[1,24])';
PAR_absorptance = repmat(X1(:,1).*(pmax(11)-pmin(11))+pmin(11),[1,24])';
NIR_absorptance = repmat(X1(:,1).*(pmax(12)-pmin(12))+pmin(12),[1,24])';
emissivity = repmat(X1(:,1).*(pmax(13)-pmin(13))+pmin(13),[1,24])';
ca = repmat(400*ones(24,1),[1,N]); %ambient CO2 ppm

out0 = BESS_revised(PAR,Ta,Ps,u,RH,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca);
Tc = mean(out0.TcSu,2)-273.15;
Tc_std = 2*std(out0.TcSu,0,2);

LAI = repmat(0.01*ones(24,1),[1,N]);
out1 = BESS_revised(PAR,Ta,Ps,u,RH,LAI,CI,hc,Vcmax25,g1,PAR_absorptance,NIR_absorptance,emissivity,ca);
Ts = mean(out1.Ts,2)-273.15;

figure;
plot(smooth(Ta,5),'LineWidth',1.5);hold on;
plot(smooth([Ts(2:end); Ts(1)],5),'LineWidth',1.5);
% plot(smooth([Tc(2:end); Tc(1)],5),'r-','LineWidth',2)
errorbar(1:24,smooth([Tc(2:end); Tc(1)],5),Tc_std,'r-','LineWidth',1)
ylim([5 28]);
set(gcf,'position',[500,500,650*0.5,400*0.5])
print(gcf, 'D:\Data\Global Thermoregulation\For RSE\Figures\figure_05_grassland.jpg', '-djpeg', '-r600');
