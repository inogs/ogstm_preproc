% Script to create netCDF4 profiles of WOA2018 data
% (N1p, N3n, N5s, O2o)

clear
clc
close all

%%%%%% EDIT
strPathWOA='DATA_WOA2018';
OUTPUT_DIR='PROFILES'; 
do_graphics=1; % 0
%%%%%%

strDATAwoa  = [strPathWOA,'/woa18_all_'];
strTempData = [strPathWOA,'/woa18_decav_t00_01.nc'];
strSalData =  [strPathWOA,'/woa18_decav_s00_01.nc'];

if ~exist(OUTPUT_DIR, 'dir')
       mkdir(OUTPUT_DIR)
end

if do_graphics && ~exist('PLOTS', 'dir')
       mkdir('PLOTS')
end

varWoaList={'p','n','i','o'};
varNameList={'N1p', 'N3n', 'N5s', 'O2o'};

%%%%%%%% CYCLE ON ALL VARIABLES

for ivar=1:4;

varWoa=varWoaList{1,ivar}; 
varName=varNameList{1,ivar}; 

%%% read data
WOA_lon = double(ncread([strDATAwoa,varWoa,'00_01.nc'],'lon'));
WOA_lat = double(ncread([strDATAwoa,varWoa,'00_01.nc'],'lat'));
depth = double(ncread([strDATAwoa,varWoa,'00_01.nc'],...
    'depth'));
WOA_ann = double(ncread([strDATAwoa,varWoa,'00_01.nc'],...
   [varWoa,'_mn']));

% selection of the lon-lat boxes for the computation
WOA_lonMin=172;
WOA_lonMax=174;
WOA_latMin=123;
WOA_latMax=130;

%%% computation of profiles in all boxes

WOA_nlonAtl=WOA_lonMax-WOA_lonMin+1;
WOA_nlatAtl=WOA_latMax-WOA_latMin+1;
WOA_ndepth=size(depth,1);
WOA_ann_atl=NaN(WOA_ndepth,WOA_nlonAtl*WOA_nlatAtl);
count=1;
for ilon=1:WOA_nlonAtl
   for ilat=1:WOA_nlatAtl
       WOA_ann_atl(:,count)=...
           squeeze(WOA_ann(ilon+WOA_lonMin-1,ilat+WOA_latMin-1,:));
       count=count+1;
   end
end

clear WOA_ann 

%%% oxygen case: units in micromol/kg has to be converted in mmol/m^3
if strcmp(varWoa,'o')     
    
% 1) read temperature and salinity data
    
T_ann = double(ncread(strTempData,'t_mn'));
S_ann = double(ncread(strSalData,'s_mn'));
  
% 2) computation in all boxes

T_ann_atl=NaN(WOA_ndepth,WOA_nlonAtl*WOA_nlatAtl);
S_ann_atl=T_ann_atl;
count=1;
for ilon=1:WOA_nlonAtl
   for ilat=1:WOA_nlatAtl
       T_ann_atl(:,count)=...
           squeeze(T_ann(ilon+WOA_lonMin-1,ilat+WOA_latMin-1,:));
       S_ann_atl(:,count)=...
           squeeze(S_ann(ilon+WOA_lonMin-1,ilat+WOA_latMin-1,:));
       count=count+1;
   end
end

% 3) rho computation in all boxes

rho_pot_ann=zeros(WOA_ndepth,WOA_nlonAtl*WOA_nlatAtl);
%rho_situ_ann=zeros(WOA_ndepth,WOA_nlonAtl*WOA_nlatAtl);
depth_mat=repmat(depth,1,WOA_nlonAtl*WOA_nlatAtl);
for iz=1:WOA_ndepth
rho_pot_ann(iz,:)=densjmd95_sq(S_ann_atl(iz,:),...
    T_ann_atl(iz,:),zeros(size(T_ann_atl(iz,:))));
%rho_situ_ann(iz,:)=densjmd95_sq(S_ann_atl(iz,:),...
%    T_ann_atl(iz,:),depth_mat(iz,:));
end

clear S_ann* T_ann* S_seas* T_seas*
    
% 4) conversion micromol/kg -> mmol/m^3

WOA_ann_atl=WOA_ann_atl.*rho_pot_ann/1000;     
end
%%% end case of oxygen


%%% mean on all boxes
WOA_atlProf(:,1)=nanmean(WOA_ann_atl,2);

fclose('all');


%%%%%%%%% INTERPOLATION

nav_lev=double(ncread('/share/data/vdibiagio/COPERNICUS/reanalisi/forICandBC/masks/REA24/meshmask.nc',...
    'nav_lev'));

profile=zeros(125,1);
x=WOA_atlProf;
nanx = isnan(x);
t    = 1:numel(x);
if nansum(~nanx)>1 % at least 2 points
x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
end
profile = interp1(depth,x,nav_lev,'linear','extrap');   

% find last not-NaN value
BBB = ~isnan(profile); 
if nansum(profile)>0
 indx = arrayfun(@(x) find(~isnan(profile), 1, 'last'), 1);
 Index=indx;
 Value=profile(indx,1);
end

% Final profiles
profile_Def=profile;
profile_Def(Index:125, 1)=Value;

unitsVar=' [mmol m^-3]';

%%%%%%%%%% GRAPHICS (OPTIONAL)

if do_graphics

hf=figure('Position',[100   100   660*3   520*2.1]);
plot(profile_Def,-1.0*nav_lev,'*b','MarkerSize',9);  
hold on
plot(WOA_atlProf,-1.0*depth,'MarkerEdgeColor','r',...
        'MarkerFaceColor','r','Marker','o');
set(gca,'FontSize',16);  
title([varName,unitsVar], 'FontSize',18);
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
print(gcf,'PLOTS/WOA_profile','-djpeg','-r500');

end

%%%%%%% WRITING OF netCDF4 PROFILES

ncid = netcdf.create([OUTPUT_DIR,'/',varName, '.nc'], 'NETCDF4');
nz = netcdf.defDim(ncid,'nav_lev',125);
varid = netcdf.defVar(ncid,varName,'NC_DOUBLE', nz);
netcdf.putAtt(ncid,varid,'units',unitsVar);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,profile_Def);
netcdf.close(ncid);

end 
%%%%%%%% END CYCLE ON ALL VARIABLES



