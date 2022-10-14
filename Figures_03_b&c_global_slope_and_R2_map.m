scale_dn = read(Tiff('D:\Data\Global Thermoregulation\scale_8d_v2.tif','r'));
igbp = read(Tiff('D:\Data\Global Thermoregulation\igbpLandCover.tif','r'));
corr_dn = read(Tiff('D:\Data\Global Thermoregulation\corr_8d_v2.tif','r'));
load('D:\Data\Global Thermoregulation\satellite data_v4\background.mat');
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

% remove the tropical and sparse vegetation region
sizes_data = size(LAI_gsmean);
LAI_gsmean(round(sizes_data(1)*0.4):round(sizes_data(1)*0.6),:)=nan;
plantfraction(plantfraction<0.95)=nan;
plantfraction(~isnan(plantfraction))=1;
LAI_gsmean(LAI_gsmean<22)=nan;
LAI_gsmean(~isnan(LAI_gsmean))=1;
waterfraction(waterfraction>5)=nan;
waterfraction(~isnan(waterfraction))=1;
remain_area = plantfraction.*LAI_gsmean.*single(waterfraction);

scale_dn(scale_dn<0.2)=nan;
scale_dn(scale_dn>2)=nan;
scale_dn = scale_dn.*remain_area;

figure;histogram(scale_dn(~isnan(scale_dn)), 50, 'Normalization','PDF')
Y0 = scale_dn(1:30:sizes_data(1),1:30:sizes_data(2));
Y = blockproc(scale_dn, [30 30], fun);

figure;histogram(Y(~isnan(Y)), 40, 'Normalization','PDF','FaceColor',[0.09 0.79 0.78])
xlim([0.5 1.5]);
set(gcf,'position',[500,500,65*2,65*2])
print(gcf, 'D:\Data\Global Thermoregulation\For New Phytologist\Figures\figure_02_slope_hist.jpg', '-djpeg', '-r600');

background(end-round(size(Y0,1)/7):end,:)=[];
Y(end-round(size(Y0,1)/7):end,:)=[];

figure;ax1 = axes;
A = imagesc(background,[0,1.1]);
set(A,'AlphaData',~isnan(background))
ax2 = axes;
B = imagesc(Y,[0.5 1.5]); 
set(B,'AlphaData',~isnan(Y))
linkaxes([ax1,ax2]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'gray');
colormap(ax2,'parula');
set(ax2,'color','none','visible','off');
set(gcf,'position',[500,500,650*1.2,300*1.2])
print(gcf, 'D:\Data\Global Thermoregulation\For New Phytologist\Figures\figure_02_slope_colorbar.jpg', '-djpeg', '-r600');

%% corr
figure;
histogram(corr_dn(~isnan(corr_dn)), 50, 'Normalization','PDF')

r2 = corr_dn.^2;
r2 = r2.* remain_area;
r2_fill = blockproc(r2, [30 30], fun);

figure;histogram(r2_fill(~isnan(r2_fill)), 40, 'Normalization','PDF','FaceColor',[0.09 0.79 0.78])
xlim([0.6 1.05]);
set(gcf,'position',[500,500,65*2,65*2])
print(gcf, 'D:\Data\Global Thermoregulation\For New Phytologist\Figures\figure_02_R2_hist.jpg', '-djpeg', '-r600');

figure; r2_plot = imagesc(r2_fill,[0 1]);
set(r2_plot,'AlphaData',~isnan(r2_fill))
colorbar

r2_fill(end-round(size(r2_fill,1)/7):end,:)=[];
figure;ax1 = axes;
A = imagesc(background,[0,1.1]);
set(A,'AlphaData',~isnan(background))
ax2 = axes;
B = imagesc(r2_fill,[0.6 1.0]); 
set(B,'AlphaData',~isnan(r2_fill))
linkaxes([ax1,ax2]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1,'gray');
colormap(ax2,'parula');
set(ax2,'color','none','visible','off');
set(gcf,'position',[500,500,650*1.2,300*1.2])
print(gcf, 'D:\Data\Global Thermoregulation\For New Phytologist\Figures\figure_02_R2_colorbar.jpg', '-djpeg', '-r600');

%% igbp
igbp0=single(igbp);

%% Needle forest:1,3; Broad forest:2,4; 5:Mixed; Shrub:6,7; Sav: 8,9; Gras:10; CRO:12
slope = scale_dn;
figure;histogram(slope(~isnan(slope)),40, 'Normalization','PDF','FaceColor',[0.09 0.79 0.78])

igbp2 = igbp0+slope+slope;
slope2 = slope+igbp0+igbp0;% same nan value
igbp_vec = igbp2(~isnan(igbp2));
slope_vec = slope2(~isnan(slope2));
unique(round(igbp_vec))

% mean
NF_mean = mean(slope_vec(igbp_vec==1 | igbp_vec==3));
NF_std = std(slope_vec(igbp_vec==1 | igbp_vec==3));

BF_mean = mean(slope_vec(igbp_vec==2 | igbp_vec==4));
BF_std = std(slope_vec(igbp_vec==2 | igbp_vec==4));

MF_mean = mean(slope_vec(igbp_vec==5));
MF_std = std(slope_vec(igbp_vec==5));

SAV_mean = mean(slope_vec(igbp_vec==8 | igbp_vec==9));
SAV_std = std(slope_vec(igbp_vec==8 | igbp_vec==9));

SHR_mean = mean(slope_vec(igbp_vec==6 | igbp_vec==7));
SHR_std = std(slope_vec(igbp_vec==6 | igbp_vec==7));

GRA_mean = mean(slope_vec(igbp_vec==10));
GRA_std = std(slope_vec(igbp_vec==10));

CRO_mean = mean(slope_vec(igbp_vec==12));
CRO_std = std(slope_vec(igbp_vec==12));

slope_satellite_PFT = [BF_mean BF_std; NF_mean NF_std; MF_mean MF_std; SAV_mean SAV_std; SHR_mean SHR_std; GRA_mean GRA_std; CRO_mean CRO_std];
