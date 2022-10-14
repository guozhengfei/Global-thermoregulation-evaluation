function out = BESS_revised(PAR,Ta,Ps,u,RH,LAI,CI,hc,vcmax25,g1,PAR_absorptance,NIR_absorptance, emissivity_leaf,ca)
%% This function is based on the original code of BESS model (Ryu et al., 2011);link:https://www.environment.snu.ac.kr/bess-flux
% dir add dif PAR [umol m-2 s-1]
PAR = PAR*4.56;
% Diffuse photosynthetically active radiation [umol m-2 s-1]
% Source: BESS radiation product
PARDiff = PAR*0.25;
% Direct photosynthetically active radiation [umol m-2 s-1]
% Source: BESS radiation product
PARDir = PAR - PARDiff;
% dir add dif NIR [umol m-2 s-1]
NIR = PAR; %*0.48/52
% Diffuse near infrared radiation [umol m-2 s-1]
% Source: BESS radiation product
NIRDiff = NIR.*0.2;
% Direct near infrared radiation [umol m-2 s-1]
% Source: BESS radiation product
NIRDir = NIR - NIRDiff;

% Surface pressue [Pa]
% Source: ERA5_land product
Ps = Ps*1000;
% Air temperature [oC] --> [K]
% Source: MOD07_L2 product, with pressue level correction and cross-calibration using ERA Interim Reanalysis 2 metre temperature dataset
Ta = Ta + 273.15;
Ta(Ta>323) = 323;
% Dewpoint temperature [K]
% Source: ERA5 product
% Transfer RH to Td
T = Ta -273.15;
Td = 243.04*(log(RH/100)+((17.625*T)./(243.04+T)))./(17.625-log(RH/100)-((17.625*T)./(243.04+T)))+273.15;
% cos_SZA Solar zenith angle
% Source: MCD19A2_GRANULES
cosdSZA = cos(30*pi/180);
% Land surface temperature [K]
% Source: MOD11A1 product

% Wind speed [m s-1]
% Source: ERA5 U-wind at 10 m and V-wind at 10 m datasets
WS = u;
% Leaf area index [m2 m-2]
% Source: MOD15A2 product, with temporally filtered (Appendix in Jiang and Ryu 2016)
LAI = LAI*1; %LAI should be 0~10
% Black-sky visible albedo [-]
% Source: MCD43B3 product, with gap filled 
ALB_BSA_VIS =(1-PAR_absorptance)*(1/1.95)*2;
% Black-sky near infrared albedo [-] (direct beam albedo)
% Source: MCD43B3 product, with gap filled 
ALB_BSA_NIR = 1-NIR_absorptance;
% White-sky visible albedo [-] (diffusion beam albedo)
% Source: MCD43B3 product, with gap filled 
%ALB_WSA_VIS = ALB_WSA_VIS_glb(~isnan(mask))*0.001;
ALB_WSA_VIS = ALB_BSA_VIS*0.95*2; % calculate according to the global fluxnet2015 site albedo
% White-sky near infrared albedo [-]
% Source: MCD43B3 product, with gap filled 
%ALB_WSA_NIR = ALB_WSA_NIR_glb(~isnan(mask))*0.001;
% Climate region classification
% Source: Köppen climate classification
% Climate = importdata('Input/Climate.FLUXNET2015.mat');
% Climate = repmat(Climate,365*length(years),1);
% Land cover classification
% Source: MCD12Q1 (IGPB), with cropland/Natural vegetation mosaic replaced by secondary type
% Digital elevation model [m]
% Source: SRTM GTOPO 30 Arc Second
% DEM = importdata('Input/DEM.FLUXNET2015.mat');
% DEM = repmat(DEM,365*length(years),1);
% Clumping index [-]
% Source: Jingming Chen's group
CI = CI*1;
% C4 plant fraction [-]
% Source: Still's global distribution of C3 and C4 vegetation dataset, only applied for grasslands and croplands
% fC4 = importdata('Input/fC4.FLUXNET2015.mat');
% fC4 = repmat(fC4,365*length(years),1);
% fC4(Landcover~=10 & Landcover~=12) = 0;
fC4 = zeros(size(Ta)); %
% Canopy height [m]
% Source: GLAS/ICESat global forest canopy height dataset, only applied for forest
hc = hc*1;
% Annual maximum and minimum leaf area index [m2 m-2]
% Source: derived from LAI data
% LAI_Max = lai_95_glb(~isnan(mask));
% LAI_Min = lai_05_glb(~isnan(mask));
% Peak carbonxylation capacity at 25 degree for C3 and C4 leaf [umol m-2 s-1]
% Source: plant functional type dependent look-up table (Appendix in Jiang and Ryu 2016)
% peakVcmax25_C3Leaf = importdata('Input/peakVcmax25_C3Leaf.FLUXNET2015.mat');
% peakVcmax25_C4Leaf = importdata('Input/peakVcmax25_C4Leaf.FLUXNET2015.mat');
Vcmax25_C3Leaf = vcmax25;
g1 = g1*1;
% Ambient carbon dioxide concentration [ppm]
% Source: OCO-2 XCO2 dataset and NOAA CO2 growth rate dataset
% XCO2 = importdata('Input/XCO2.FLUXNET2015.mat');
% Ca = 400 + zeros(size(Ta));
Ca = ca + zeros(size(Ta));
% list = [1,31;32,59;60,90;91,120;121,151;152,181;182,212;213,243;244,273;274,304;305,334;335,365];
% for mon = 1:12
%     temp = repmat(XCO2(mon,:),list(mon,2)-list(mon,1)+1,1);
%     Ca(list(mon,1):list(mon,2),:) = temp;
% end
% Ca = repmat(Ca,length(years),1);  
% dXCO2 = [1.85,2.36,2.28,1.56,2.45,1.74,2.11,1.77,1.67,2.39,1.69,2.40,2.51,1.90,2.73];
% dCa = nan([1,length(years)]);
% for year = years
%     dCa(year-1999) = sum(dXCO2(15-(2015-year)+1:end));
% end
% dCa = repmat(reshape(repmat(dCa,365,1),[],1),1,length(name));
% Ca = Ca - dCa;
% 
% LAT = repmat(lats',365*length(years),1); 
% LON = repmat(lons',365*length(years),1); 
% DOY = repmat([1:365]',length(years),length(name));
% YEAR = repmat(reshape(repmat(years,365,1),365*length(years),1),1,length(name));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The core module of BESS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% m_CanopyShortwaveRadiation
%%            

% Lookup table for leaf scattering coefficient and soil reflectance (Table S1 in Ryu et al 2011)
LUT = [[0.12 0.15 0.12 0.15 0.12 0.175 0.15 0.1 0.1 0.175 0.175 0.175 NaN 0.175 NaN 0.1];
       [0.45 0.7 0.45 0.7 0.55 0.83 0.7 0.62 0.62 0.83 0.83 0.83 NaN 0.83 NaN 0.62].*0.9;    % 0.9 is to consider UV (Youngryel)
       [0.11 0.11 0.11 0.11 0.11 0.11 0.25 0.25 0.25 0.25 0.11 0.1 NaN 0.1 NaN 0.25];
       [0.23 0.23 0.23 0.23 0.23 0.23 0.42 0.42 0.42 0.42 0.23 0.15 NaN 0.15 NaN 0.42].*0.9];
% Leaf scattering coefficient for PAR

sigma_P = 0.1418 + zeros(size(LAI)); % add by me
sigma_N = 0.6146 + zeros(size(LAI)); % add by me
rho_PSoil = 0.1586 + zeros(size(LAI)); % add by me
rho_NSoil = 0.2578 + zeros(size(LAI)); % add by me

% for i = 1:16
%     sigma_P(Landcover==i) = LUT(1,i); % mean of first row in LUT :0.1418
% end 
% % Soil reflectance for PAR
% rho_PSoil = single(nan(size(LAI)));
% for i = 1:16
%     rho_PSoil(Landcover==i) = LUT(3,i);
% end  
% % Leaf scattering coefficient for NIR
% sigma_N = single(nan(size(LAI)));
% for i = 1:16
%     sigma_N(Landcover==i) = LUT(2,i);
% end    
% % Soil reflectance for NIR
% rho_NSoil = single(nan(size(LAI)));
% for i = 1:16
%     rho_NSoil(Landcover==i) = LUT(4,i);
% end  
  
% Beam radiation extinction coefficient of canopy
kb = 0.5 ./ cosdSZA;    % Table A1 in Ryu et al 2011
% Extinction coefficient for beam and scattered beam PAR
kk_Pb = 0.46 ./ cosdSZA;    % Table A1 in Ryu et al 2011
% Extinction coefficient for diffuse and scattered diffuse PAR
kk_Pd = 0.72;    % Table A1 in Ryu et al 2011
% Extinction coefficient for beam and scattered beam NIR
kk_Nb = kb .* sqrt(1-sigma_N);    % Table A1 in Ryu et al 2011
% Extinction coefficient for diffuse and scattered diffuse NIR
kk_Nd = 0.35 * sqrt(1-sigma_N);    % Table A1 in Ryu et al 2011
% Nitrogen extinction coefficient
kn = kk_Pb;    % Table A1 in Ryu et al 2011

%% PAR

% Total absorbed incoming PAR
Q_PDn = (1-ALB_BSA_VIS).*PARDir.*(1-exp(-kk_Pb.*LAI.*CI)) + (1-ALB_WSA_VIS).*PARDiff.*(1-exp(-kk_Pd.*LAI.*CI));    % Eq. (2) in Ryu et al 2011
% Absorbed incoming beam PAR by sunlit leaves
Q_PbSunDn = PARDir.*(1-sigma_P).*(1-exp(-kb.*LAI.*CI));    % Eq. (3) in Ryu et al 2011
% Absorbed incoming diffuse PAR by sunlit leaves
Q_PdSunDn = PARDiff.*(1-ALB_WSA_VIS).*(1-exp(-(kk_Pd+kb).*LAI.*CI)).*kk_Pd./(kk_Pd+kb);    % Eq. (4) in Ryu et al 2011
% Absorbed incoming scattered PAR by sunlit leaves
Q_PsSunDn = PARDir.*((1-ALB_BSA_VIS).*(1-exp(-(kk_Pb+kb).*LAI.*CI)).*kk_Pb./(kk_Pb+kb)-(1-sigma_P).*(1-exp(-2.*kb.*LAI.*CI))/2);    % Eq. (5) in Ryu et al 2011
Q_PsSunDn(Q_PsSunDn<0) = 0;
% Absorbed incoming PAR by sunlit leaves
Q_PSunDn = Q_PbSunDn + Q_PdSunDn + Q_PsSunDn;    % Eq. (6) in Ryu et al 2011
% Absorbed incoming PAR by shade leaves
Q_PShDn = Q_PDn - Q_PSunDn;    % Eq. (7) in Ryu et al 2011
Q_PShDn(Q_PShDn<0) = 0;
% Incoming PAR at soil surface
I_PSoil = (1-ALB_BSA_VIS).*PARDir + (1-ALB_WSA_VIS).*PARDiff - (Q_PSunDn+Q_PShDn);
% Absorbed PAR by soil
Q_PSoil = (1-rho_PSoil) .* I_PSoil;
% Absorbed outgoing PAR by sunlit leaves
Q_PSunUp = I_PSoil .* rho_PSoil .* exp(-kk_Pd.*LAI.*CI);    % Eq. (8) in Ryu et al 2011
% Absorbed outgoing PAR by shade leaves
Q_PShUp = I_PSoil .* rho_PSoil .* (1-exp(-kk_Pd.*LAI.*CI));    % Eq. (9) in Ryu et al 2011
% Total absorbed PAR by sunlit leaves
Q_PSun = Q_PSunDn + Q_PSunUp;    % Eq. (10) in Ryu et al 2011
% Total absorbed PAR by shade leaves
Q_PSh = Q_PShDn + Q_PShUp;    % Eq. (11) in Ryu et al 2011

APAR_Sun = Q_PSun; 
APAR_Sh = Q_PSh;
APAR_Soil = Q_PSoil;

%% NIR

L_CI = LAI.*CI;
exp__kk_Nd_L_CI = exp(-kk_Nd.*L_CI);

% Absorbed incoming NIR by sunlit leaves
Q_NSunDn = NIRDir.*(1-sigma_N).*(1-exp(-kb.*L_CI)) + NIRDiff.*(1-ALB_BSA_NIR).*(1-exp(-(kk_Nd+kb).*L_CI)).*kk_Nd./(kk_Nd+kb) + NIRDir.*((1-ALB_BSA_NIR).*(1-exp(-(kk_Nb+kb).*L_CI)).*kk_Nb./(kk_Nb+kb)-(1-sigma_N).*(1-exp(-2*kb.*L_CI))/2);    % Eq. (14) in Ryu et al 2011
% Absorbed incoming NIR by shade leaves
Q_NShDn = (1-ALB_BSA_NIR).*NIRDir.*(1-exp(-kk_Nb.*L_CI)) + (1-ALB_BSA_NIR).*NIRDiff.*(1-exp__kk_Nd_L_CI) - Q_NSunDn;    % Eq. (15) in Ryu et al 2011
% Incoming NIR at soil surface
I_NSoil = (1-ALB_BSA_NIR).*NIRDir + (1-ALB_BSA_NIR).*NIRDiff - (Q_NSunDn+Q_NShDn);
% Absorbed NIR by soil
Q_NSoil = (1-rho_NSoil) .* I_NSoil;
% Absorbed outgoing NIR by sunlit leaves
Q_NSunUp = I_NSoil .* rho_NSoil .* exp__kk_Nd_L_CI;    % Eq. (16) in Ryu et al 2011
% Absorbed outgoing NIR by shade leaves
Q_NShUp = I_NSoil .* rho_NSoil .* (1-exp__kk_Nd_L_CI);    % Eq. (17) in Ryu et al 2011
% Total absorbed NIR by sunlit leaves
Q_NSun = Q_NSunDn + Q_NSunUp;    % Eq. (18) in Ryu et al 2011
% Total absorbed NIR by shade leaves
Q_NSh = Q_NShDn + Q_NShUp;     % Eq. (19) in Ryu et al 2011

ANIR_Sun = Q_NSun;
ANIR_Sh = Q_NSh;
ANIR_Soil = Q_NSoil;

%%

%%
% m_CanopyMaxCarboxylationCapacity
%%

kb = 0.5 ./ cosdSZA;
kn = 0.713;
kn_kb_Lc = kn + kb.*LAI;
LAI_Vcmax25_C3Leaf = LAI .* Vcmax25_C3Leaf;
Vcmax25_C3Tot =  LAI_Vcmax25_C3Leaf .* (1-exp(-kn)) ./ kn;    % Eq. (32) in Ryu et al 2011
Vcmax25_C3Sun = LAI_Vcmax25_C3Leaf .* (1-exp(-kn_kb_Lc)) ./ kn_kb_Lc;    % Eq. (33) in Ryu et al 2011
Vcmax25_C3Sh = Vcmax25_C3Tot - Vcmax25_C3Sun;    % Eq. (34) in Ryu et al 2011
% LAI_Vcmax25_C4Leaf = LAI .* Vcmax25_C4Leaf;
% Vcmax25_C4Tot =  LAI_Vcmax25_C4Leaf .* (1-exp(-kn)) ./ kn;    % Eq. (32) in Ryu et al 2011
% Vcmax25_C4Sun = LAI_Vcmax25_C4Leaf .* (1-exp(-kn_kb_Lc)) ./ kn_kb_Lc;    % Eq. (33) in Ryu et al 2011
% Vcmax25_C4Sh = Vcmax25_C4Tot - Vcmax25_C4Sun;    % Eq. (34) in Ryu et al 2011

%%

%%
% m_Thermodynamics
%%

% ratio molecular weight of water vapour dry air
mv_ma = 0.622;    % [-] (Wiki)

% latent heat of vaporization
lambda = 1.91846e6 * (Ta./(Ta-33.91)).^2;    % [J kg-1] (Henderson-Sellers, 1984)

% ambient vapour pressure
ea = 2.1718e10 * exp(-4157./(Td-33.91));    % [Pa] (Henderson-Sellers, 1984)
% ea = 0.611*exp(17.502*(Ta+273.15)./(Ta+273.15 +240.97))*1000; % [Pa] (Campbell and Norman., 2012)
% saturated vapour pressure
es = 2.1718e10 * exp(-4157./(Ta-33.91));    % [Pa] (Henderson-Sellers, 1984)

% water vapour deficit
VPD = es - ea;    % [Pa]

% relative humidity
RH = ea ./ es;    % [-]

% 1st derivative of saturated vapour pressure
desTa = 4157 * es .* (Ta-33.91).^(-2);    % [pa K-1]

% 2nd derivative of saturated vapour pressure
ddesTa = 4157 * (desTa.*(Ta-33.91).^(-2)+(-2)*es.*(Ta-33.91).^(-3));    % [pa K-2]

% specific humidity
q = (mv_ma.*ea) ./ (Ps-0.378*ea);    % [-] (Garratt, 1994)

% air density
rhoa = Ps ./ (287.05.*Ta);    % [kg m-3] (Garratt, 1994)

% specific heat of dry air 
Cpd = 1005 + (Ta-250).^2 / 3364;    % [J kg-1 K-1] (Garratt, 1994)

% specific heat of air 
Cp = Cpd .* (1+0.84*q);    % [J kg-1 K-1] (Garratt, 1994)

% ratio_molecular_weight_of_water_vapour_dry_air
a = 0.622;    % [-] (Wiki)

% Psychrometric constant
gamma = Cp.*Ps./(a*lambda);    % [pa K-1]

% Precipitable water
eta = 0.465 * ea ./ Ta;    % [] (Prata, 1996)

% Clear-sky emissivity
epsa = 1 - (1+eta).*exp(-(1.2+3.0*eta).^0.5);    % [-] (Prata, 1996)

% Constraints
VPD(VPD<0)=0;
RH(RH>1)=1;
RH(RH<0)=0;
rhoa(rhoa<0)=0;
Cp(Cp<0)=0;
gamma(gamma<0)=0;
desTa(desTa<0)=0;
ddesTa(ddesTa<0)=0;

%%

%%
% m_CarbonWaterFluxes
%%    

% Extinction coefficient for beam PAR or NIR for black leaves
kb = 0.5 ./ cosdSZA;    % Table A1 in Ryu et al 2011
fSun = CI .* exp(-kb.*LAI.*CI);    % Eq. (1) in Ryu et al 2011
fSum_T = CI .* exp(-kb.*LAI*0.1.*CI);  % add by me, divide canopy to 10 layer, we foucus on the first layer for Tcan calculate
epss = 0.94;    % Table A1 in Ryu et al 2011
epsf = emissivity_leaf;    % Table A1 in Ryu et al 2011
mskC4 = fC4 > 0.01;

% Initialization
Tf_C3Sun = Ta;
Tf_C3Sh = Ta;
Tf_C4Sun = Ta;
Tf_C4Sh = Ta;
Ts = Ta;
Tf = Ta;
Ci_C3Sun = 0.7 * Ca + zeros(size(Ta),'single');
Ci_C3Sh = 0.7 * Ca + zeros(size(Ta),'single');
Ci_C4Sun = 0.4 * Ca + zeros(size(Ta),'single');
Ci_C4Sh = 0.4 * Ca + zeros(size(Ta),'single');

% Initial estimation of ra
k = 0.4;    % von Kármán constant
z0 = hc * 0.1;    % (Norman)
ustar = WS.*k./(log(10./z0));
ra = WS./ustar.^2+2./(k.*ustar);    % Eq. (2-4) in Ryu et al 2008
ra(ra>100) = 100;
% for wet plant
fw = RH.^4;
fw(fw<=0.6^4)=0;
%rhc = 1.0/(0.01*LAI*fwet);
% Iteration    
for iter = 1:3

    %%
    % m_CanopyLongwaveRadiation
    %%
    
    % Extinction coefficient for longwave radiation
    kd = 0.78;    % Table A1 in Ryu et al 2011
    
    % Stefan_Boltzmann_constant
    sigma = 5.670373e-8;    % [W m-2 K-4] (Wiki)

    % Long wave radiation flux densities from air, soil and leaf
    La = epsa * sigma .* Ta.^4;    % Stefan Boltzmann law
    Ls = epss * sigma .* Ts.^4;
    Lf = epsf * sigma .* Tf.^4;

    % Absorbed longwave radiation by sunlit leaves
    ALW_Sun = (Ls-Lf).*kd.*(exp(-kd.*LAI)-exp(-kb.*LAI))./(kd-kb) + kd.*(La-Lf).*(1-exp(-(kb+kd).*LAI))./(kd+kb);    % Eq. (44) in Kowalczyk et al 2006
    % Absorbed longwave radiation by shade leaves
    ALW_Sh = (1-exp(-kd.*LAI)) .* (Ls+La-2*Lf) - ALW_Sun;    % Eq. (45) in Kowalczyk et al 2006

    %%
    
    % C3_Sunlit, C3_Shade, C4_Sunlit, C4_Shade
    for i = 1:2
        if i == 1
            flgC3 = true;
            flgSun = true;
            Tf = Tf_C3Sun;
            Ci = Ci_C3Sun;
            Vcmax25 = Vcmax25_C3Sun;
            APAR = APAR_Sun;
            ANIR = ANIR_Sun;
            ALW = ALW_Sun;
        elseif i == 2
            flgC3 = true;
            flgSun = false;
            Tf = Tf_C3Sh;
            Ci = Ci_C3Sh;
            Vcmax25 = Vcmax25_C3Sh;
            APAR = APAR_Sh;
            ANIR = ANIR_Sh;
            ALW = ALW_Sh;
%         elseif i == 3
%             flgC3 = false;
%             flgSun = true;
%             Tf = Tf_C4Sun;
%             Ci = Ci_C4Sun;
%             Vcmax25 = Vcmax25_C4Sun;
%             APAR = APAR_Sun;
%             ANIR = ANIR_Sun;
%             ALW = ALW_Sun;
%         else
%             flgC3 = false;
%             flgSun = false;
%             Tf = Tf_C4Sh;
%             Ci = Ci_C4Sh;
%             Vcmax25 = Vcmax25_C4Sh;
%             APAR = APAR_Sh;
%             ANIR = ANIR_Sh;
%             ALW = ALW_Sh;
        end
%         if i >= 3
%             Tf(~mskC4) = nan;
%             Ci(~mskC4) = nan;
%             Vcmax25(~mskC4) = nan;
%             APAR(~mskC4) = nan;
%             ANIR(~mskC4) = nan;
%             ALW(~mskC4) = nan;
%         end
        
        if flgC3
        
            %%
            % m_C3Photosynthesis
            %%
            
            % O2 concentration
            O = 209460;    % [umol mol-1] 
            % Gas constant
            R = 8.314e-3;    % [kJ K-1 mol-1] 

            % Unit convertion
            Ci = Ci*1e-6 .* Ps;    % [umol mol-1] -> [Pa]
            O = O*1e-6 .* Ps;    % [umol mol-1] -> [Pa]

            % Kinetic constants and temperature correction (von Caemmerer et al., 2009)
            KC25 = 267*1e-6*1e5;    % [ubar] -> [Pa]
            KCEa = 80.99;    % [kJ mol-1]
            KO25 = 164*1e-3*1e5;    % [mbar] -> [Pa]
            KOEa = 23.72;    % [kJ mol-1]
            GammaS25 = 36.9*1e-6*1e5;    % [ubar] -> [Pa]
            GammaSEa = 24.6;    % [kJ mol-1]
            KC = KC25 .* exp(KCEa*(Tf-298)./(R*Tf*298));
            KO = KO25 .* exp(KOEa*(Tf-298)./(R*Tf*298));
            GammaS = GammaS25 .* exp(GammaSEa*(Tf-298)./(R*Tf*298));
            VcmaxEa = 58.52;    % [kJ mol-1]
            JmaxEa = 37;    % [kJ mol-1]
            RdEa = 66.4;    % [kJ mol-1]
            H = 220;    % [kJ mol-1]
            S = 0.710;    % [kJ K-1 mol-1]
            Vcmax = Vcmax25 .* exp(VcmaxEa*(Tf-298)./(R*Tf*298)) .* (1+exp((S*298-H)/(R*298))) ./ (1+exp((S*Tf-H)./(R*Tf)));
            Jmax = 3.6*Vcmax25 .* exp(JmaxEa*(Tf-298)./(R*Tf*298)) .* (1+exp((S*298-H)/(R*298))) ./ (1+exp((S*Tf-H)./(R*Tf)));    % 3.6 is from (Groenendijk et al., 2011)
            Rd = 0.0089*Vcmax25 .* exp(RdEa*(Tf-298)./(R*Tf*298)) .* (1+exp((S*298-H)/(R*298))) ./ (1+exp((S*Tf-H)./(R*Tf)));    % 0.0089 is from (De Pury and Farquhar, 1997)

            Kf = 0.05;    % Rate constant for fluorescence
            Kd = 0.95;    % Rate constant for thermal deactivation at Fm
            Kp = 4.0;    % Rate constant for photochemisty
            po0 = Kp ./ (Kp+Kd+Kf);    % Dark photochemistry fraction (Genty et al., 1989)
            alf = 0.5 * po0;    % Quantum yield 
            J = (alf.*APAR.*Jmax) ./ (alf.*APAR+2.1*Jmax);    % Electron transport rate (von Caemmerer and Farquhar, 1981)

            % Farquhar's model
            C_GammaS = Ci - GammaS;
            coefJ = C_GammaS./(4*(Ci+2*GammaS));
            Ac = Vcmax .* C_GammaS./(Ci+KC.*(1+O./KO));
            Aj = J .* coefJ;
            a = 0.99;    % Co-limitation effect is not essential for canopy level (De Pury and Farquhar, 1997)
            b = -(Aj+Ac);
            c = Aj.*Ac;
            Ag = (-b + sign(b) .* sqrt(b.^2 - 4.*a.*c))./(2.*a);
            An = Ag - Rd;
            An(An<0) = 0;

            % SIF (from Joe Berry's code)
            Ja = Ag ./ coefJ;
            ps = po0 .* Ja./J;
            x = 1 - ps./po0;
            Kn = (6.2473*x-0.5944) .* x;
            fm = Kf./(Kf+Kd+Kn);
            fs = fm .* (1-ps);    % Fluorescence yield
            SIF = APAR .* fs;
            
            %%

        else
        
            %%
            % m_C4Photosynthesis
            %%

            % Partial pressure of O2 (Collatz et al., 1991)
            O2 = 20.9 * 1e3;    % [Pa] 

            Tref = 298;
            slti = 0.2;
            shti = 0.2;
            Thl = 288;
            Thh = 303;
            Trdm = 328;
            qt = 0.1 * (Tf-Tref);
            TH = 1 + exp((-220E3+703*Tf)./(8.314*Tf));
            TL = 1 + exp(slti.*(Thl-Tf));

            Rd = 0.015 * Vcmax25 .* 1.8.^qt ./(1+exp(1.3*(Tf-Trdm)));
            Vcmax = Vcmax25 .* 2.1.^qt ./(TL.*TH);

            spfy = 2600 * 0.75 .^qt;            % This is, in theory, Vcmax/Vomax.*Ko./Kc, but used as a separate parameter
            gam = 0.5 ./spfy .*O2; %[bar]       compensation point [bar]
            kp = Vcmax25/56*1E6.* 1.8.^qt;

            Kf = 0.05;
            Kd = 0.95;
            Kp = 4.0;
            po0 = Kp./(Kf+Kd+Kp);         % dark photochemistry fraction (Genty et al., 1989)
            Je = 0.5*po0 .* APAR;          % electron transport rate

            effcon = 0.17;                    % Berry and Farquhar (1978): 1/0.167
            Ve = Je.*(Ci-gam)./(Ci+2*gam) .* effcon;
            Vs = kp.*Ci;
            a = 0.8;
            b = -(Vcmax+Ve);
            c = Vcmax.*Ve;
            V = (-b + sign(b) .* sqrt(b.^2 - 4.*a.*c))./(2.*a).*(Ci>gam);
            a = 0.98;
            b = -(V+Vs);
            c = V.*Vs;
            Ag = (-b + sign(b) .* sqrt(b.^2 - 4.*a.*c))./(2.*a);
            An = Ag - Rd;

            % SIF
            Ja = Ag ./((Ci-gam)./(Ci+2*gam))./effcon;        % actual electron transport rate
            ps = po0.*Ja./Je;               % this is the photochemical yield
            x = 1 - ps./po0;
            Kn = (6.2473*x-0.5944) .* x;
            fm = Kf./(Kf+Kd+Kn);
            fs = fm .* (1-ps);
            SIF = APAR .* fs;
            
        end    
        
        %%
        % m_Conductance
        %%
        
        % Ball-Berry constant
        if flgC3
            m = 10;    % [-] Slope of Ball-Berry equation
            b0 = 0.01;    % [mol m-2 s-1] Intercept of Ball-Berry equation
        else
            m = 4;
            b0 = 0.04;
        end
        
        % Convert factor
        cf = 0.446 * (273.15./Tf).*(Ps/101325);

        % Boundary layer H2O conductance
        gb = 1./ra*1e2 .* cf;    % [mol m-2 s-1]

        % Leaf surface CO2 concentration
        Cs = Ca - An./gb;    % [umol./mol]
        
        % Stomatal H2O conductance
%         gs = m*RH.*An./Cs + b0;    % [mol m-2 s-1]
        gs = 1.6*(1+g1./sqrt(VPD*0.001)).*An./Cs + 0.1*b0; % added by me
        % Intercellular CO2 concentration
        Ci = Cs - 1.6*An./gs;    % [umol./mol]
        
        % Stomatal resistance to vapour transfer from cell to leaf surface
        rs = 1./(gs./cf*1e-2);    % [s m-1]
        rs(rs<0) = 0;

        %%
        
        %%
        % m_EnergyBalance
        %%
                            
        % Stefan_Boltzmann_constant
        sigma = 5.670373e-8;    % [W m-2 K-4] (Wiki)

        % Canopy net radiation
        %Rn = APAR/4.5 + ANIR + ALW - 4*0.98*sigma*Ta.^3.*(Tf-Ta);
        Rn = APAR/4.56 + ANIR/4.56 + ALW - 4*emissivity_leaf*sigma.*Ta.^3.*(Tf-Ta);
        % To reduce redundant computation
        ddesTa_ra2 = ddesTa .* ra.^2;
        gamma_ra_rc = gamma .* (ra+rs);
        rhoa_Cp_gamma_ra_rc = rhoa .* Cp .* gamma_ra_rc;

        % Solution
        a = 1/2 .* ddesTa_ra2./rhoa_Cp_gamma_ra_rc;    % Eq. (10b) in Paw and Gao (1988)
        b = -1 - ra.*desTa./gamma_ra_rc - ddesTa_ra2.*Rn./rhoa_Cp_gamma_ra_rc;    % Eq. (10c) in Paw and Gao (1988)
        c = rhoa .* Cp./gamma_ra_rc.*VPD + desTa.*ra./gamma_ra_rc.*Rn + 1/2*ddesTa_ra2./rhoa_Cp_gamma_ra_rc.*Rn.^2;    % Eq. (10d) in Paw and Gao (1988)
        LE = (-b+sign(b).*sqrt(b.^2-4*a.*c))./(2.*a);    % Eq. (10a) in Paw and Gao (1988)
        LE = real(LE);
        LE(LE>Rn) = Rn(LE>Rn);
        LE(Rn<0) = 0;
        LE(LE<0) = 0;
        LE(Ta<273.15) = 0;

        H = Rn - LE;
        dT = ra./(rhoa.*Cp) .* H;    % Eq. (6) in Paw and Gao (1988)
        dT(dT>10) = 10;
%         dT(dT<0) = 0;
        Tf = Ta + dT;
        %%

        if i == 1
            An_C3Sun = An;
            Rn_C3Sun = Rn;
            LE_C3Sun = LE;
            H_C3Sun = H;
            Tf_C3Sun = Tf;
            Ci_C3Sun = Ci;
            SIF_C3Sun = SIF;
            gs_C3Sun = 1./rs;
        elseif i == 2
            An_C3Sh = An;
            Rn_C3Sh = Rn;
            LE_C3Sh = LE;
            H_C3Sh = H;
            Tf_C3Sh = Tf;
            Ci_C3Sh = Ci;
            SIF_C3Sh = SIF;
            gs_C3Sh = 1./rs;
%         elseif i == 3
%             An_C4Sun = An;
%             Rn_C4Sun = Rn;
%             LE_C4Sun = LE;
%             H_C4Sun = H;
%             Tf_C4Sun = Tf;
%             Ci_C4Sun = Ci;
%             SIF_C4Sun = SIF;
%             gs_C4Sun = 1./rs;
%         else
%             An_C4Sh = An;
%             Rn_C4Sh = Rn;
%             LE_C4Sh = LE;
%             H_C4Sh = H;
%             Tf_C4Sh = Tf;
%             Ci_C4Sh = Ci;
%             SIF_C4Sh = SIF;
%             gs_C4Sh = 1./rs;
        end
        
    end

    % Composite
    An_Sun = An_C3Sun;
    Rn_Sun = Rn_C3Sun;
    LE_Sun = LE_C3Sun;
    H_Sun = H_C3Sun;
    Tf_Sun = Tf_C3Sun;
    SIF_Sun = SIF_C3Sun;
    Ci_Sun = Ci_C3Sun;
    gs_Sun = gs_C3Sun;
%     An_Sun(mskC4) = An_C3Sun(mskC4).*(1-fC4(mskC4)) + An_C4Sun(mskC4).*fC4(mskC4); 
%     Rn_Sun(mskC4) = Rn_C3Sun(mskC4).*(1-fC4(mskC4)) + Rn_C4Sun(mskC4).*fC4(mskC4); 
%     LE_Sun(mskC4) = LE_C3Sun(mskC4).*(1-fC4(mskC4)) + LE_C4Sun(mskC4).*fC4(mskC4); 
%     H_Sun(mskC4) = H_C3Sun(mskC4).*(1-fC4(mskC4)) + H_C4Sun(mskC4).*fC4(mskC4); 
%     Tf_Sun(mskC4) = Tf_C3Sun(mskC4).*(1-fC4(mskC4)) + Tf_C4Sun(mskC4).*fC4(mskC4); 
%     SIF_Sun(mskC4) = SIF_C3Sun(mskC4).*(1-fC4(mskC4)) + SIF_C4Sun(mskC4).*fC4(mskC4); 
%     Ci_Sun(mskC4) = Ci_C3Sun(mskC4).*(1-fC4(mskC4)) + Ci_C4Sun(mskC4).*fC4(mskC4); 
%     gs_Sun(mskC4) = gs_C3Sun(mskC4).*(1-fC4(mskC4)) + gs_C4Sun(mskC4).*fC4(mskC4); 
    An_Sh = An_C3Sh;
    Rn_Sh = Rn_C3Sh;
    LE_Sh = LE_C3Sh;
    H_Sh  = H_C3Sh;
    Tf_Sh = Tf_C3Sh;
    SIF_Sh = SIF_C3Sh;
    Ci_Sh = Ci_C3Sh;
    gs_Sh = gs_C3Sh;
%     An_Sh(mskC4) = An_C3Sh(mskC4).*(1-fC4(mskC4)) + An_C4Sh(mskC4).*fC4(mskC4); 
%     Rn_Sh(mskC4) = Rn_C3Sh(mskC4).*(1-fC4(mskC4)) + Rn_C4Sh(mskC4).*fC4(mskC4); 
%     LE_Sh(mskC4) = LE_C3Sh(mskC4).*(1-fC4(mskC4)) + LE_C4Sh(mskC4).*fC4(mskC4); 
%     H_Sh(mskC4) = H_C3Sh(mskC4).*(1-fC4(mskC4)) + H_C4Sh(mskC4).*fC4(mskC4); 
%     Tf_Sh(mskC4) = Tf_C3Sh(mskC4).*(1-fC4(mskC4)) + Tf_C4Sh(mskC4).*fC4(mskC4); 
%     SIF_Sh(mskC4) = SIF_C3Sh(mskC4).*(1-fC4(mskC4)) + SIF_C4Sh(mskC4).*fC4(mskC4); 
%     Ci_Sh(mskC4) = Ci_C3Sh(mskC4).*(1-fC4(mskC4)) + Ci_C4Sh(mskC4).*fC4(mskC4); 
%     gs_Sh(mskC4) = gs_C3Sh(mskC4).*(1-fC4(mskC4)) + gs_C4Sh(mskC4).*fC4(mskC4); 

%     Tf = Tf_Sun.*fSun + Tf_Sh.*(1-fSun);
    out.TcSu = Tf_Sun;
    out.TcSh = Tf_Sh;
    Tf = Tf_Sun.*fSun + Tf_Sh.*(1-fSun);
    %%
    % m_Soil
    %%
    
%     ra = 0.5 * ra; % BESS original
    ra = 1.5 * ra; % add by me
    
    % Stefan_Boltzmann_constant
    sigma = 5.670373e-8;    % [W m-2 K-4] (Wiki)
    
    % Long wave radiation flux densities from air, soil and leaf
    La = epsa * sigma .* Ta.^4;
    Ls = epss * sigma .* Ts.^4;
    Lf = epsf * sigma .* Tf.^4;

    % Long wave radiation flux densities absorbed by soil
    ALW = (1-exp(-kd.*LAI)).*Lf + exp(-kd.*LAI).*La;    % Eq. (41) in Kowalczyk et al 2006
    % Net radiation on soil surface
    Rn = APAR_Soil/4.5 + ANIR_Soil/4.5 + ALW - Ls - 4*epsa*sigma.*Ta.^3.*(Ts-Ta);

    % Ground heat
    G = 0.35 * Rn;    % Eq. (39) in Ryu et al 2011
    G(G<0) = 0;
    % Latent heat
    LE = desTa./(desTa+gamma) .*(Rn-G) .* RH.^(VPD/1000);    % Eq. (37) in Ryu et al 2011
    LE(LE>Rn) = Rn(LE>Rn);
    LE(Rn<0) = 0;
    LE(LE<0) = 0;
    LE(Ta<273.15) = 0;
    % Sensible heat
    H = Rn - G - LE;

    % Update temperature
    dT = ra./(rhoa.*Cp) .* H;
    dT(dT>10) = 10;
%     dT(dT<0) = 0;
    Ts = Ta + dT;
    
    Rn_Soil = Rn;
    LE_Soil = LE;
    H_Soil = H;

    %%

    % Composite components
    Rn = Rn_Sun + Rn_Sh + Rn_Soil;    % [W m-2] 
    LE = LE_Sun + LE_Sh + LE_Soil;    % [W m-2] 
    H = H_Sun + H_Sh + H_Soil;    % [W m-2]  rhoa
    GPP = An_Sun + An_Sh;    % [umol m-2 s-1]
    ET = LE ./ lambda;    % [g m-2 s-1] = [um s-1]
    
    %% 
    % m_Aerodynamics
    %%

    % von Karman's constant 
    k = 0.4;    % [-]
    % gravitational acceleration
    g = 9.8;    % [m s-2]
    % Monin–Obukhov length
    L = -rhoa.*Cp.*ustar.^3.*Ta./(k.*g.*H);    % [m]
    
    % Businger-Dyer function (Dyer, 1974)
    z_L = 10 ./ L;
    z_L(z_L>0.25)=0.25;
    z_L(z_L<-3)=-3;
    phim1 = nan(size(z_L),'single');
    phim1(z_L>=0) = 1 + 5*z_L(z_L>=0);
    phim1(z_L<0) = (1-16*z_L(z_L<0)).^(-0.5);
    z_L = z_L.*(z0./10);
    phim2 = nan(size(z_L),'single');
    phim2(z_L>=0) = 1 + 5*z_L(z_L>=0);
    phim2(z_L<0) = (1-16*z_L(z_L<0)).^(-0.5);

    % u(z) ~ u*
    ustar = WS.*k./(log(10./z0)-phim1+phim2);    % [m s-1] (Louis, 1979)
    ra = WS./ustar.^2+2./(k.*ustar);
    ra(ra>100) = 100;        

end  

%%

%%
% m_TemporalUpscaling
% %%
% 
% %SZA(SZA>90) = nan;
% 
% % Solar constant
% Ssc = 1362;    % [W/m-2] (Wiki)
% td = DOY;
% 
% % Time zone
% UTC = LON / 15;
% % Decimal time with half an hour interval
% t = (0+0.25):0.5:(24-0.25);
% 
% % Actual top of atmosphere solar radiation
% RgTOA_Snapshot = Ssc * (1+0.033*cos(2*pi*td/365)) .* cosdSZA;    % Eq. (1) in Ryu et al 2008
% 
% % Potential top of atmosphere solar radiation for each half an hour
% RgPOT_t = nan(size(SZA,1),size(SZA,2),length(t));
% for i = 1:length(t)
%     % Solar zenith angle
%     beta = sun_angle(YEAR, DOY, t(i), LAT, LON, UTC);
%     % Calculate for current t
%     RgPOT_t(:,:,i) = Ssc * (1+0.033*cos(2*pi*td/365)) .* cosd(beta);    % Eq. (1) in Ryu et al 2008
% end
% % Constraints
% RgPOT_t(RgPOT_t<0) = 0;
% 
% % Ratio of half-hourly sum to daily sum
% SFd = 1800*RgTOA_Snapshot ./ squeeze(nanmean(RgPOT_t,3)*60*60*24);    % Eq. (2) in Ryu et al 2008
% SFd(SFd<0) = nan;
% 
% % Upscale from snapshot to daily
% GPP_Daily = 1800 * GPP ./ SFd * 1e-6*12;    % Eq. (3) in Ryu et al 2008
% GPP_Daily(SFd<0.01) = 0;
% GPP_Daily(SZA>=90) = 0;
% LE_Daily = 1800 * LE ./ SFd * 1e-6;    % Eq. (3) in Ryu et al 2008
% LE_Daily(SFd<0.01) = 0;
% LE_Daily(SZA>=90) = 0;
% ET_Daily = 1800 * ET ./ SFd * 1;    % Eq. (3) in Ryu et al 2008
% ET_Daily(SFd<0.01) = 0;
% ET_Daily(SZA>=90) = 0;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     
% %%
% % Make outputs
%     
% GPP_8day = reshape(nanmean(reshape([reshape(GPP_Daily,365,length(years)*length(name));nan(3,length(years)*length(name))],8,46*length(years)*length(name)),1),46*length(years),length(name));    % [gC m-2 d-1]
% LE_8day = reshape(nanmean(reshape([reshape(LE_Daily,365,length(years)*length(name));nan(3,length(years)*length(name))],8,46*length(years)*length(name)),1),46*length(years),length(name));    % [MJ m-2 d-1]
% ET_8day = reshape(nanmean(reshape([reshape(ET_Daily,365,length(years)*length(name));nan(3,length(years)*length(name))],8,46*length(years)*length(name)),1),46*length(years),length(name));    % [mm d-1]
% 
% data = GPP_8day; save('Output/BESSGPP.FLUXNET2015.mat','data');
% data = LE_8day; save('Output/BESSLE.FLUXNET2015.mat','data');
% data = ET_8day; save('Output/BESSET.FLUXNET2015.mat','data');
out.Tc = Tf;
out.GPP = GPP;
out.VPD = VPD;
out.ET = ET;
out.Ts = Ts; 
