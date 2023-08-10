Ta_10_01 = read(Tiff('Y:\Download\Ta_10-0000000000-0000000000.tif','r'));
Ta_10_02 = read(Tiff('Y:\Download\Ta_10-0000000000-0000001280.tif','r'));

Ta_11_01 = read(Tiff('Y:\Download\Ta_11-0000000000-0000000000.tif','r'));
Ta_11_02 = read(Tiff('Y:\Download\Ta_11-0000000000-0000001280.tif','r'));

Ta_12_01 = read(Tiff('Y:\Download\Ta_12-0000000000-0000000000.tif','r'));
Ta_12_02 = read(Tiff('Y:\Download\Ta_12-0000000000-0000001280.tif','r'));

Ta_13_01 = read(Tiff('Y:\Download\Ta_13-0000000000-0000000000.tif','r'));
Ta_13_02 = read(Tiff('Y:\Download\Ta_13-0000000000-0000001280.tif','r'));

Ta_14_01 = read(Tiff('Y:\Download\Ta_14-0000000000-0000000000.tif','r'));
Ta_14_02 = read(Tiff('Y:\Download\Ta_14-0000000000-0000001280.tif','r'));

% TcNoon_01_info = geotiffinfo('Y:\Download\Ta_10-0000000000-0000000000.tif');
% TcNoon_02_info = geotiffinfo('Y:\Download\Ta_10-0000000000-0000001280.tif');

Ta_10 = single(cat(2,Ta_10_01,Ta_10_02));
Ta_11 = single(cat(2,Ta_11_01,Ta_11_02));
Ta_12 = single(cat(2,Ta_12_01,Ta_12_02));
Ta_13 = single(cat(2,Ta_13_01,Ta_13_02));
% Ta_13 = (Ta_12+Ta_14)/2;
Ta_14 = single(cat(2,Ta_14_01,Ta_14_02));

latlim = TcNoon_01_info.SpatialRef.LatitudeLimits;
lonlim = [min(TcNoon_01_info.SpatialRef.LongitudeLimits), max(TcNoon_02_info.SpatialRef.LongitudeLimits)];
rasterSize = [size(Ta_10,1), size(Ta_10,2)];
R = georefcells(latlim,lonlim,rasterSize,'ColumnsStartFrom','north');

% geotiffwrite('Ta_10.tif',Ta_10,R)

clear Ta_10_01 Ta_10_02 Ta_11_01 Ta_11_02 Ta_12_01 Ta_12_02 Ta_13_01 Ta_13_02 Ta_14_01 Ta_14_02

Ta_10_half = (Ta_10+Ta_11)*0.5;
Ta_11_half = (Ta_11+Ta_12)*0.5;
Ta_12_half = (Ta_12+Ta_13)*0.5;
Ta_13_half = (Ta_13+Ta_14)*0.5;

%read MODIS Pass time
MOD_time = read(Tiff('Y:\Download\MOD_passTime.tif','r'));
MYD_time = read(Tiff('Y:\Download\MYD_passTime.tif','r'));
MOD_info = geotiffinfo('Y:\Download\MOD_passTime.tif');

bottom = (min(TcNoon_01_info.SpatialRef.LatitudeLimits)-min(MOD_info.SpatialRef.LatitudeLimits) ...
    )/TcNoon_01_info.SpatialRef.CellExtentInLatitude;

top = abs(max(TcNoon_01_info.SpatialRef.LatitudeLimits)-max(MOD_info.SpatialRef.LatitudeLimits) ...
    )/TcNoon_01_info.SpatialRef.CellExtentInLatitude;
MOD_time(1:top,:)=[];
MOD_time(end-bottom+1:end,:)=[];

MYD_time(1:top,:)=[];
MYD_time(end-bottom+1:end,:)=[];

% imagesc(Ta_10(:,:,1));
% figure;plot(reshape(Ta_10(151,1641,:),[460,1]));hold on;
% plot(reshape(Ta_10_half(151,1641,:),[460,1]));
% plot(reshape(Ta_11(151,1641,:),[460,1]));
% plot(reshape(Ta_11_half(151,1641,:),[460,1]));
% plot(reshape(Ta_12(151,1641,:),[460,1]));
% plot(reshape(Ta_12_half(151,1641,:),[460,1]));
% plot(reshape(Ta_13(151,1641,:),[460,1]));
% plot(reshape(Ta_13_half(151,1641,:),[460,1]));
% plot(reshape(Ta_14(151,1641,:),[460,1]));


Ta_10_reshape = reshape(Ta_10,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_10half_reshape = reshape(Ta_10_half,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_11_reshape = reshape(Ta_11,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_11half_reshape = reshape(Ta_11_half,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_12_reshape = reshape(Ta_12,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_12half_reshape = reshape(Ta_12_half,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_13_reshape = reshape(Ta_13,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_13half_reshape = reshape(Ta_13_half,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);
Ta_14_reshape = reshape(Ta_14,[size(Ta_10,1)*size(Ta_10,2),size(Ta_10,3)]);

Ta_all = cat(3,Ta_10_reshape,Ta_10half_reshape,Ta_11_reshape,Ta_11half_reshape,Ta_12_reshape,Ta_12half_reshape, ...
    Ta_13_reshape,Ta_13half_reshape,Ta_14_reshape);
MOD_time_reshape = reshape(MOD_time,[size(MOD_time,1)*size(MOD_time,2),1]);
MYD_time_reshape = reshape(MYD_time,[size(MOD_time,1)*size(MOD_time,2),1]);

clear Ta_10_half Ta_10_reshape Ta_10half_reshape Ta_11 Ta_11_reshape Ta_11_half Ta_11half_reshape
clear Ta_12 Ta_12_half Ta_12_reshape Ta_12half_reshape Ta_13 Ta_13_reshape Ta_13_half Ta_13half_reshape
clear Ta_14 Ta_14_reshape

MOD_time_int = 1+(floor(MOD_time_reshape/10)-10)*2+ round(mod(MOD_time_reshape,10)/5);
MOD_time_int(MOD_time_int<1)=nan;MOD_time_int(MOD_time_int>9)=nan;
MYD_time_int = 1+(floor(MYD_time_reshape/10)-10)*2+ round(mod(MYD_time_reshape,10)/5);
MYD_time_int(MYD_time_int<1)=nan;MYD_time_int(MYD_time_int>9)=nan;

% figure;histogram(MOD_time_int(~isnan(MOD_time_int)))
% imagesc(reshape(MYD_time_int,[724,1994]))

parfor i=1:length(MOD_time_reshape)
    ta_slice = Ta_all(i,:,:);
    index1 = MOD_time_int(i,1);
    index2 = MYD_time_int(i,1);
    if isnan(ta_slice(1,1,1)*index1*index2)
        Ta_final(i,:) = ta_slice(1,:,1); 
    else
        ta2 = (ta_slice(1,:,index1)+ta_slice(1,:,index2))/2;
        Ta_final(i,:) = ta2;
    end
    if mod(i,1994*7)==0
        i/(1994*7)
    end
end

Ta_match = reshape(Ta_final,[724,1994,460]);
figure; imagesc(Ta_match(:,:,23))
Ta_match = int16(Ta_match*100);
geotiffwrite('TaNoon_GEE.tif',Ta_match,R)