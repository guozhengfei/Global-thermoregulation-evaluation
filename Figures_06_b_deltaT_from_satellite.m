igbp = read(Tiff('D:\Data\Global Thermoregulation\igbpLandCover.tif','r'));
load('D:\Data\Global Thermoregulation\satellite data_v4\background.mat');
dT = read(Tiff('D:\Data\Global Thermoregulation\dT_8d_v2.tif','r'));
LAI_gsmean = read(Tiff('D:\Data\Global Thermoregulation\LAI_gs_mean.tif','r'));
plantfraction = read(Tiff('D:\Data\global plant fraction\global_plant_fraction.tif','r'));
waterfraction = read(Tiff('D:\Data\Global Thermoregulation\water_fraction.tif','r'));
%% background display as grey
igbp0 = im2single(igbp);
igbp0(igbp0==0 | igbp0>0.0666)= nan;
fun = @(block_struct) mean(block_struct.data(:),'omitnan');
igbp1 = blockproc(igbp0, [30 30], fun);
background = igbp1;
background(~isnan(background))=1;

sizes_data = size(LAI_gsmean);
LAI_gsmean(round(sizes_data(1)*0.4):round(sizes_data(1)*0.6),:)=nan;
plantfraction(plantfraction<0.95)=nan;
plantfraction(~isnan(plantfraction))=1;
LAI_gsmean(LAI_gsmean<22)=nan;
LAI_gsmean(~isnan(LAI_gsmean))=1;
waterfraction(waterfraction>5)=nan;
waterfraction(~isnan(waterfraction))=1;
remain_area = plantfraction.*LAI_gsmean.*single(waterfraction);
%% dT
dT = dT.*remain_area;
figure;histogram(dT(~isnan(dT)),40, 'Normalization','PDF','FaceColor', [0.09 0.79 0.78])

Y = blockproc(dT, [30 30], fun);
figure;histogram(Y(~isnan(Y)), 20, 'Normalization','PDF','FaceColor', [0.09 0.79 0.78])
ylim([0 0.35]);
set(gcf,'position',[500,500,65*2,65*2])
print(gcf, 'D:\Data\Global Thermoregulation\For New Phytologist\Figures\figure_05_dT_satellite.jpg', '-djpeg', '-r600');

background(end-round(size(Y,1)/7):end,:)=[];
Y(end-round(size(Y,1)/7):end,:)=[];

figure;ax1 = axes;
A = imagesc(background,[0,1.1]);
set(A,'AlphaData',~isnan(background))
ax2 = axes;
B = imagesc(Y,[0 6]); 
set(B,'AlphaData',~isnan(Y))
% print(gcf, 'D:\Data\Global Thermoregulation\For RSE\Figures\figure_05_dT_satellite_spatial_colorbar.jpg', '-djpeg', '-r600');
linkaxes([ax1,ax2]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'gray');
colormap(ax2,'parula');
set(ax2,'color','none','visible','off');
set(gcf,'position',[500,500,650*1.2,300*1.2])
print(gcf, 'D:\Data\Global Thermoregulation\For New Phytologist\Figures\figure_05_dT_satellite_spatial.jpg', '-djpeg', '-r600');

%% igbp
igbp0=single(igbp);
dT2 = dT;
igbp2 = igbp0+dT2-dT2;
dT3 = dT2+igbp0-igbp0;% same nan value
igbp_vec = igbp2(~isnan(igbp2));
dT_vec = dT3(~isnan(dT3));
unique(round(igbp_vec))
% Needle forest:1,3; Broad forest:2,4; 5:Mixed; Shrub:6,7; Sav: 8,9; Gras:10; CRO:12
dT_igbp = cat(2,dT_vec,igbp_vec);
% mean
NF_mean = mean(dT_vec(igbp_vec==1 | igbp_vec==3));
NF_std = std(dT_vec(igbp_vec==1 | igbp_vec==3));

BF_mean = mean(dT_vec(igbp_vec==2 | igbp_vec==4));
BF_std = std(dT_vec(igbp_vec==2 | igbp_vec==4));

MF_mean = mean(dT_vec(igbp_vec==5));
MF_std = std(dT_vec(igbp_vec==5));

SAV_mean = mean(dT_vec(igbp_vec==8 | igbp_vec==9));
SAV_std = std(dT_vec(igbp_vec==8 | igbp_vec==9));

SHR_mean = mean(dT_vec(igbp_vec==6 | igbp_vec==7));
SHR_std = std(dT_vec(igbp_vec==6 | igbp_vec==7));

GRA_mean = mean(dT_vec(igbp_vec==10));
GRA_std = std(dT_vec(igbp_vec==10));

CRO_mean = mean(dT_vec(igbp_vec==12));
CRO_std = std(dT_vec(igbp_vec==12));

dT_satellite_PFT = [BF_mean BF_std; NF_mean NF_std; MF_mean MF_std; SAV_mean SAV_std; SHR_mean SHR_std; GRA_mean GRA_std; CRO_mean CRO_std]; 
