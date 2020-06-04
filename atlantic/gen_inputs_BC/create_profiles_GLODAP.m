%%%%%%% Script to create netCDF4 profiles of GLODAP data
% (O3c, O3h)

clc
clear
close all

%%%%%%%% EDIT
GLODAP_DIR_DATA='DATA_GLODAP_GIB'; 
OUTPUT_DIR='PROFILES'; 
do_graphics=1; % 0
%%%%%%%%

if ~exist(OUTPUT_DIR, 'dir')
       mkdir(OUTPUT_DIR)
end

if do_graphics && ~exist('PLOTS', 'dir')
       mkdir('PLOTS')
end

%%%%%%% READ GLODAP DATA

formatSpec = '%f %f';
sizeMat=[2 Inf];

varList={'O3c','O3h'};
unitsList={' [mgC m^-3]',' [mmol m^-3]'};
nvar=size(varList,2);

depth=[];
allVars=[];
for iv=1:nvar
varName=varList{1,iv};
filename=[GLODAP_DIR_DATA,'/',varName,'.txt'];
fileID = fopen(filename,'r');
Mat = fscanf(fileID,formatSpec,sizeMat);
Mat(Mat==-999)=NaN;
Mat=Mat';
depth=[depth Mat(:,1)];
allVars=[allVars Mat(:,2)];
fclose(fileID);
end

depth_gl=depth(:,1);
nLevs=size(depth_gl,1);


%%%%%%% INTERPOLATION ON RAN24 VERTICAL LEVELS

nav_lev=double(ncread('/share/data/vdibiagio/COPERNICUS/reanalisi/forICandBC/masks/REA24/meshmask.nc',...
    'nav_lev'));

profInterp=NaN(125,nvar);
MAT_x=zeros(nLevs,nvar);

for ip=1:nvar
    %%% interpolation of inner NANs
    x= allVars(:,ip);
    nanx = isnan(x);
    t    = 1:numel(x);
    if nansum(~nanx)>1 % at least 2 points
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    end
    MAT_x(:,ip)=x';
    %%% linear interpolation on RAN24 levels 
    if nansum(~nanx)>1 
    pItp = interp1(depth(:,ip),x,nav_lev,'linear','extrap');                             
    profInterp(:,ip)=pItp;
    end
end

%%% Correction on deeper levels: last not-NaNs values applied down to max depth
maxZ=125;
BBB = ~isnan(MAT_x); 
Indexes=NaN(1,nvar);
Values=NaN(1,nvar);
for ip=1:nvar
	aaa=MAT_x(:,ip);
 	if nansum(aaa)>0
 		indx = arrayfun(@(x) find(~isnan(aaa), 1, 'last'), 1);% last not-NaN value 
 		Indexes(1,ip)=indx;
 		Values(1,ip)=aaa(indx,1);
 	end
end

Indices_REA=125*ones(1,nvar);
% Find correspondence between Indexes of MAT_x and indexes of nav_lev of RAN 
for ip=1:nvar
	depth_gl=depth(:,ip);
    if ~isnan(Indexes(1,ip))
	depthLC=depth_gl(Indexes(1,ip));
	[vv, idxREA] = min(abs(nav_lev-depthLC)); 
	Indices_REA(1,ip)=idxREA;
    end
end
profInterp_Def=profInterp;
for icol=1:nvar
    if Indices_REA(1,icol)~=125 
        profInterp_Def(Indices_REA(1,icol):maxZ, icol)=Values(1,icol);   
    end
end

%%% Smoothing procedure 
profInterp_Def_Smooth=zeros(125,nvar);
for ip=1:nvar
profInterp_Def_Smooth(:,ip) = smooth(nav_lev,profInterp_Def(:,ip),0.1,'moving');
end

%%%%%%%%% GRAPHICS (OPTIONAL)

if do_graphics
hf=figure('Position',[100   100   660*3   520*2.1]);
for ip=1:nvar
varName=varList{1,ip};
unitsVar=unitsList{1,ip};
values_orig=allVars(:,ip);
values_interp=profInterp_Def(:,ip);
values_interp_smoothed=profInterp_Def_Smooth(:,ip);
subplot(1,nvar,ip);
plot(values_interp, -1.0*nav_lev,'*b');
hold on
plot(values_interp_smoothed, -1.0*nav_lev, 'og');
hold on
plot(values_orig,-1.0*depth(:,ip),'or','MarkerSize',9,'MarkerFaceColor','r')
title([varName,unitsVar], 'FontSize',18);
set(gca,'FontSize',16)        
end
set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
print(gcf,'PLOTS/GLODAP_profile','-djpeg','-r500');
end

%%%%%%%% WRITING OF netCDF4 PROFILES

for ip=1:nvar
   
    varName=varList{1,ip};
    unitsVar=unitsList{1,ip};
	
    ncid = netcdf.create([OUTPUT_DIR,'/',varName, '.nc'], 'NETCDF4');
    nz = netcdf.defDim(ncid,'nav_lev',125);
    varid = netcdf.defVar(ncid,varName,'NC_DOUBLE', nz);
    netcdf.putAtt(ncid,varid,'units',unitsVar);
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid,profInterp_Def_Smooth(:,ip));
    netcdf.close(ncid);

end


 
