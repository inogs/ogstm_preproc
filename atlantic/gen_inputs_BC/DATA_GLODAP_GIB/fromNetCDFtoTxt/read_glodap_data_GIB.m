% read station extracted from GLODAP v2 dataset
% ref: olsen et al., 2019

clear all
close all
addpath([ '/Users/gcossarini/MEDITERRANEO/dataset_carbsys/CO2sys_v2_1/CO2sys_MAT/']); % use the library of seawater properties
addpath('/Users/gcossarini/matlab/proprieties_of_sea_water/')
fontsize=16;
% read global coast LOW RESOLUTION
s =  gshhs('/Users/gcossarini/MEDITERRANEO/data_coastline/gshhs_1.3/gshhs_l.b',[30 46],[-5 36])
% background color
BKG = [1 1 1];
% coast line color
CST = [.2 .2 .2];

nomi={'lat','lon','dic','alk','pho','sil','no3','temp','sal','dens','time','yyyy_mm_dd','pH','depth'}
NVAR=length(nomi);
com=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','longitude');
stat_lon=-360+com; % now in LON W (stations have values from -10W to -5.5W)

stat_lat=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','latitude');
com=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','date_time');
% date_time:units = "days since 1981-01-01 00:00:00 UTC"
stat_time=datenum(1981,1,1,0,0,0)+com;
[stat_yyyy_mm_dd]=datevec(stat_time);
%dic umol/kg
stat_dic=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var11');
%dic QC
stat_dic_qc=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var11_qc');
%alkalinity umol/kg
stat_alk=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var10');
%alk qc
stat_alk_qc=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var10_qc');
%silicat umol/kg
stat_sil=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var7');
%sil qc
stat_sil_qc=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var7_qc');
%phosph umol/kg
stat_pho=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var6');
%pho qc
stat_pho_qc=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var6_qc');
%nitraumol/kg
stat_no3=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var8');
%no8 qc
stat_no3_qc=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var8_qc');
%pH  % pHT [p,T,S] (pH in total scale at insitu condition)
stat_pH=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var30');
%pH qc
stat_pH_qc=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var30_qc');
%gamma "NEUTRAL DENSITY" KG/M**3" ;
stat_gamma=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var37');
% depth
stat_depth=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var2');
% temperatura degC
stat_temp=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var3');
% salinity
stat_sal=ncread('data_from_GLODAPv2.2019_atlanticGIB.nc','var4');
% compute density
stat_dens=sw_dens(stat_sal,stat_temp,stat_depth);


NS=length(stat_lon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qualezona=1; % ATLANTICO
O3hlim=[2380 2500]; O3clim=[25300 27500]; DIClim=[2060 2220];ALKlim=[2300 2450];
qualezona=2; % GIB
O3hlim=[2420 2680]; O3clim=[25500 29000]; DIClim=[2090 2340];ALKlim=[2360 2600];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtri
% 1111  escludo stazioni con profondita' piu' bassa di 150m (considero costiere)
canc=zeros(NS,1);
for i=1:1:NS
    [imax com] = max(stat_depth(:,i))
    if imax<150
        canc(i)=canc(i)+1;
    end
    % 222  escludo le stazioni vicino a gibilterra:: lon@>-6.5
    if qualezona==1
        ETICHETTA='ATL';
        if stat_lon(i)>-6.5
            canc(i)=canc(i)+1;
        end
    else
        ETICHETTA='GIB';
        if stat_lon(i)<-6.5
            canc(i)=canc(i)+1;
        end
    end
    % 3333  escludo le stazioni campionate prima del 1990
    if stat_yyyy_mm_dd(i,1)<1990
        canc(i)=canc(i)+1;
    end
end
canc(canc>=1)=2;
%
colori=jet(NS);
figure('position',[100 100 1000 400],'paperpositionmode','auto','renderer','zbuffer');
subplot(1,4,1);
for i=1:1:NS
    plot(stat_lon(i),stat_lat(i),'*','color',colori(i,:));hold on;
end
%=============== Add coastline ============================
for i_pol=1:size(s)
    if(s(i_pol).Area > 5000.)
        %       disp(i_pol)
        %        patchm(s(i_pol).Lat,s(i_pol).Lon,BKG);
        %        patch(s(i_pol).Lon,s(i_pol).Lat,BKG);
        %        hold on
        %         [h]=plotm(s(i_pol).Lat,s(i_pol).Lon);
        [h]=plot(s(i_pol).Lon,s(i_pol).Lat);
        set(h,'linewidth',1.2,'color',CST);
    end
end
set(gca,'xlim',[-10 -4],'ylim',[34 37.5]);
T1=title('lon,lat stations');set(T1,'fontsize',16);
set(gca,'fontsize',fontsize);

subplot(1,4,2);
for i=1:1:NS
    plot(stat_alk(:,i),stat_depth(:,i),'-','color',colori(i,:));hold on;
end
T1=title('ALK, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');

subplot(1,4,3);
for i=1:1:NS
    plot(stat_dic(:,i),stat_depth(:,i),'-','color',colori(i,:));hold on;
end
T1=title('DIC, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');

subplot(1,4,4);
for i=1:1:NS
    plot(stat_pH(:,i),stat_depth(:,i),'-','color',colori(i,:));hold on;
end
T1=title('pH, -');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');
print tuttestaz_staz_profiliDIC_ALK_pH.png -dpng -r300



%%
figure('position',[100 100 1000 400],'paperpositionmode','auto','renderer','zbuffer');
subplot(3,2,1);
for i=1:1:NS
    plot(stat_yyyy_mm_dd(i,1),stat_yyyy_mm_dd(i,2),'*','color',colori(i,:));hold on;
end
hold on;  T1=title('year vs month');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);

subplot(3,2,2);
for i=1:1:NS
    plot(i,1,'*','color',colori(i,:));hold on;
    plot(i,canc(i),'+k');
end
T1=title(' stations on map and time [CANCELLED DATA=2]');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);

%plot(-360+stat_lon,'*k');hold on;
%plot(stat_lat-30,'xr');hold on;
%title('lonW(BLACK) ,lat (+30N RED) of stations');
set(gca,'xlim',[.5 NS+.5]);

subplot(3,2,3);pcolor(stat_depth);colorbar;hold on;caxis([0 600]);T1=title('depth'); set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,4);pcolor(stat_temp);colorbar;hold on;T1=title('temp');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,5);pcolor(stat_dens);colorbar;hold on;T1=title('dens');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,6);pcolor(stat_sal);colorbar;hold on;T1=title('sal');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
print tuttestaz_var_TSdens.png -dpng -r300

figure('position',[100 100 1000 400],'paperpositionmode','auto','renderer','zbuffer');
subplot(3,2,1);pcolor(stat_no3);colorbar;hold on;T1=title('nitr umol/kg');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,2);pcolor(stat_pho);colorbar;hold on;T1=title('pho');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,3);pcolor(stat_sil);colorbar;hold on;T1=title('silic umol/kg');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,4);pcolor(stat_alk);colorbar;hold on;T1=title('alk');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,5);pcolor(stat_dic);colorbar;hold on;T1=title('dic');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,6);pcolor(stat_pH);colorbar;hold on;T1=title('pH');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
print tuttestaz_var_carb.png -dpng -r300

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%applico il filtro di esclusione dei dati
for iv=1:1:NVAR
    eval(['COM=stat_' nomi{iv} ';']);
    [a b]=size(COM);
    if a==NS & b==1
        COM(canc==2)=[];
    elseif a==NS & b>1
            COM(canc==2,:)=[];
else
        COM(:,canc==2)=[]; 
    end
    
    eval(['SI' nomi{iv} '=COM;']);
end
SINS=length(SIlat);
% rifaccio i profili
SIcolori=jet(SINS);
figure('position',[100 100 1000 400],'paperpositionmode','auto','renderer','zbuffer');
subplot(1,4,1);
for i=1:1:SINS
    plot(SIlon(i),SIlat(i),'*','color',SIcolori(i,:));hold on;
end
%=============== Add coastline ============================
for i_pol=1:size(s)
    if(s(i_pol).Area > 5000.)
        [h]=plot(s(i_pol).Lon,s(i_pol).Lat);
        set(h,'linewidth',1.2,'color',CST);
    end
end
set(gca,'xlim',[-10 -4],'ylim',[34 37.5]);
T1=title('lon,lat stations');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
subplot(1,4,2);
for i=1:1:SINS
    plot(SIalk(:,i),SIdepth(:,i),'-','color',SIcolori(i,:));hold on;
end
T1=title('ALK, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');

subplot(1,4,3);
for i=1:1:SINS
    plot(SIdic(:,i),SIdepth(:,i),'-','color',SIcolori(i,:));hold on;
end
T1=title('DIC, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');

subplot(1,4,4);
for i=1:1:SINS
    plot(SIpH(:,i),SIdepth(:,i),'-','color',SIcolori(i,:));hold on;
end
T1=title('pH, -');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');
eval(['print ' ETICHETTA '_staz_profiliDIC_ALK_pH.png -dpng -r300 ']);

%%
figure('position',[100 100 1000 400],'paperpositionmode','auto','renderer','zbuffer');
subplot(3,2,1);
for i=1:1:SINS
    plot(SIyyyy_mm_dd(i,1),SIyyyy_mm_dd(i,2),'*','color',SIcolori(i,:));hold on;
end
hold on;  T1=title('year vs month');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);

subplot(3,2,2);
for i=1:1:SINS
    plot(i,1,'*','color',SIcolori(i,:));hold on;
end
T1=title(' stations on map and time');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);

%plot(-360+stat_lon,'*k');hold on;
%plot(stat_lat-30,'xr');hold on;
%title('lonW(BLACK) ,lat (+30N RED) of stations');
set(gca,'xlim',[.5 SINS+.5]);

subplot(3,2,3);pcolor(SIdepth);colorbar;hold on;T1=title('depth'); caxis([0 600]);set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,4);pcolor(SItemp);colorbar;hold on;T1=title('temp');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,5);pcolor(SIdens);colorbar;hold on;T1=title('dens');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,6);pcolor(SIsal);colorbar;hold on;T1=title('sal');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
eval(['print ' ETICHETTA '_var_TSdens.png -dpng -r300 ']);




figure('position',[100 100 1000 400],'paperpositionmode','auto','renderer','zbuffer');
subplot(3,2,1);pcolor(SIno3);colorbar;hold on;T1=title('nitr umol/kg');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,2);pcolor(SIpho);colorbar;hold on;T1=title('pho');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,3);pcolor(SIsil);colorbar;hold on;T1=title('silic umol/kg');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,4);pcolor(SIalk);colorbar;hold on;T1=title('alk');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,5);pcolor(SIdic);colorbar;hold on;T1=title('dic');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,6);pcolor(SIpH);colorbar;hold on;T1=title('pH');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
eval(['print ' ETICHETTA '_var_carb.png -dpng -r300 ']);


%%
% ricostruico i profili di DIC partendo da quelli ALK e pH (Total scale and p t e s condition)
% uso un valore di
% calcolo un valore medio di pho e sil per lo strato 0-200 e per 200-fondo
a=SIpho(SIdepth<200); pho200=nanmean(a); a=SIpho(SIdepth>200); phobot=nanmean(a);
a=SIsil(SIdepth<200); sil200=nanmean(a); a=SIsil(SIdepth>200); silbot=nanmean(a);

SIdicric=zeros(size(SIdic))*NaN;
SIdicnew=SIdicric;
for i=1:1:SINS
    for id=1:1:24
        if ~isnan(SIdepth(id,i))
            
            k1k2c=10; % Lueker et al, 2000
            kso4c=3; % Dickson & TB of Lee 2010
            par1=SIalk(id,i);  par1type=1;
            par2=SIpH(id,i); par2type=3; pHscale=1; % total scale
            tempin=SItemp(id,i); % temp in situ
            tempout=tempin;
            sal=SIsal(id,i);
            presin=SIdepth(id,i); % uso pressione in situ
            presout=presin; % uso pressione in situ
            if isnan(SIsil(id,i))
                if SIdepth(id,i)<200, sil=sil200;else sil=silbot;end;
            else
                sil=SIsil(id,i);
            end
            if isnan(SIpho(id,i))
                if SIdepth(id,i)<200, po4=pho200;else po4=phobot;end;
            else
                po4=SIpho(id,i);
            end
            
            [Result1]=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
            SIdicric(id,i)=Result1(2);
            
            if ~isnan(SIdic(id,i)), SIdicnew(id,i)=SIdic(id,i); else, SIdicnew(id,i)=SIdicric(id,i);end;
            
        end
    end
end
% figura con i valor ricostruisti di DIc e i valori sovrapposti
figure;
subplot(3,2,1)
plot(SIdic,SIdicric,'*k'); grid on; set(gca,'xlim',[2080 2220],'ylim',[2080 2220]);
dicricerr=nanstd(SIdic(:)-SIdicric(:)); hold on; T1=title(['dic vs dicric, err=' num2str(dicricerr) 'umol/kg']); set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
subplot(3,2,2);pcolor(SIalk);colorbar;hold on;T1=title('alk');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,3);pcolor(SIpH);colorbar;hold on;T1=title('pH');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,4);pcolor(SIdic);colorbar;hold on;T1=title('dic');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,5);pcolor(SIdicric);colorbar;hold on;T1=title('dic reconstructed');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
subplot(3,2,6);pcolor(SIdicnew);colorbar;hold on;T1=title('dic ORI+RIC');set(T1,'fontsize',fontsize); set(gca,'fontsize',fontsize);
eval(['print ' ETICHETTA '_ricostruzione_DIC.png -dpng -r300']);


%%
%subplot(3,1,2)
% costruisco la colomma media delle variabili volumetriche e
SIO3c=SIdic.*SIdens/1000*12;    % umol/kg * kd/m3  /1000 * 12g
SIO3cnew=SIdicnew.*SIdens/1000*12;
SIO3h=SIalk.*SIdens/1000;


layerUP=  [00 050 100 150 200 250 300 400 0500 1000 ];
layerDOWN=[50 100 150 200 250 300 400 500 1000 2000 ];
layerCEN=(layerUP+layerDOWN)/2;
NL=length(layerUP);
for il=1:1:NL
    ind=find(SIdepth>=layerUP(il) & SIdepth<layerDOWN(il) );
    DICnew(il)=nanmean(SIdicnew(ind));
    DIC(il)=nanmean(SIdic(ind));
    ALK(il)=nanmean(SIalk(ind));
    
    O3cnew(il)=nanmean(SIO3cnew(ind));
    O3c(il)=nanmean(SIO3c(ind));
    O3h(il)=nanmean(SIO3h(ind));
    
    
end
%Figura del profilo medio di DIC e ALK e di O3c e O3h
figure('position',[100 100 1000 700],'paperpositionmode','auto','renderer','zbuffer');
subplot(1,3,1);
for i=1:1:SINS
    plot(SIlon(i),SIlat(i),'*','color',SIcolori(i,:));hold on;
end
%=============== Add coastline ============================
for i_pol=1:size(s)
    if(s(i_pol).Area > 5000.)
        [h]=plot(s(i_pol).Lon,s(i_pol).Lat);
        set(h,'linewidth',1.2,'color',CST);
    end
end
set(gca,'xlim',[-10 -4],'ylim',[34 37.5]);
T1=title('lon,lat stations'); set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'fontsize',fontsize);

subplot(2,3,2);
for i=1:1:SINS
    plot(SIalk(:,i),SIdepth(:,i),'.-','color',SIcolori(i,:));hold on;
end
plot(ALK,layerCEN,'o-k','linewidth',2);
set(gca,'ylim',[0 1600],'xlim',ALKlim); grid on;
T1=title('ALK, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gca,'ydir','reverse');
set(gca,'fontsize',fontsize);

subplot(2,3,3);
for i=1:1:SINS
    plot(SIdicnew(:,i),SIdepth(:,i),'.-','color',SIcolori(i,:));hold on;
end
T1=title('DIC, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'fontsize',fontsize);
plot(DIC,layerCEN,'x-b','linewidth',2);
plot(DICnew,layerCEN,'o-k','linewidth',2);
set(gca,'ydir','reverse');
set(gca,'ylim',[0 1600],'xlim',DIClim); grid on;
T1=title('DIC, umol/kg');set(T1,'fontsize',fontsize);
set(gca,'ydir','reverse');
set(gca,'fontsize',fontsize);

subplot(2,3,5);
for i=1:1:SINS
    plot(SIO3h(:,i),SIdepth(:,i),'.-','color',SIcolori(i,:));hold on;
end
plot(O3h,layerCEN,'o-k','linewidth',2);
T1=title('O3h, mmol/m3');set(T1,'fontsize',fontsize);
set(gca,'ylim',[0 1600],'xlim',O3hlim); grid on;
set(gca,'ydir','reverse');
set(gca,'fontsize',fontsize);

subplot(2,3,6);
for i=1:1:SINS
    plot(SIO3cnew(:,i),SIdepth(:,i),'.-','color',SIcolori(i,:));hold on;
end
title('O3c, umol/kg');
p(1)=plot(O3c,layerCEN,'x-b','linewidth',2); leg{1}='ori';
p(2)=plot(O3cnew,layerCEN,'o-k','linewidth',2);   leg{2}='ori+rec';
set(gca,'ydir','reverse');
set(gca,'ylim',[0 1600],'xlim',O3clim); grid on;
T1=title('O3c, mg/m3'); set(T1,'fontsize',fontsize);
set(gca,'ydir','reverse');set(gca,'fontsize',fontsize);
L1=legend(p,leg); set(L1,'fontsize',16);
eval(['print ' ETICHETTA '_O3c_O3h.png -dpng -r300']);

%%
% salvo i valori di O3c e O3h
COM=[layerCEN; O3h; O3cnew]';
save([ETICHETTA '_depth_O3h_O3c.txt'],'-ascii','COM');





