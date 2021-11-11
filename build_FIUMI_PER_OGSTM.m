

% read excel file of runoff (produced by INGV)


% contain the mean runoff monthly
RUNOFF_INGV=xlsread('runoff_eas2_v12.xlsx','Foglio1');

NOMI_INGV={'Nile','Nile','Ebro','Po Maistra','Po Tramontana','Po Dritta','Po Scirocco + Po Bonifazi','Po Bastimento','Po Bocca Tolle','Po Gnocca','Po Goro','Po Levante','Po Volano', ...
           'Grand Rhone','Petit Rhone','Vjiose','Seman','Buna/Bojana','Buna/Bojana','Dardanelles','Dardanelles','Dardanelles','PIAVE','TAGLIAMENTO','SOCA/ISONZO','LIVENZA','BRENTA-BACCHIGLIONE', ...
'ADIGE','LIKA','RENO','KRKA','ARNO','NERETVA','AUDE','TREBISJNICA','TEVERE/TIBER','TEVERE/TIBER','MATI','VOLTURNO','SHKUMBINI','STRUMA/STRYMONAS','MERIC/EVROS/MARITSA','AXIOS/VARDAR','ARACHTOS','PINIOS',...
'ACHELOOS','GEDIZ','BUYUK MENDERES','KOPRU','MANAVGAT','SEYHAN','CEYHAN','GOKSU','MEDJERDA','ASI/ORONTES'};
BAS_INGV={'lev3','lev3','nwm','adr1','adr1','adr1','adr1','adr1','adr1','adr1','adr1','adr1','adr1','nwm','nwm','adr2','adr2','adr2','adr2','dard','dard','dard','adr1','adr1','adr1','adr1','adr1','adr1','adr1', ...
    'adr1','adr1','tyr1','adr1','nwm','adr2','tyr1','tyr1','adr2','tyr2','adr2','adr2','aeg','aeg','aeg','aeg','aeg','aeg','aeg','lev2','lev2','lev4','lev4','lev4','tyr2','lev4'};
isal=1; ilon=2; ilat=3; im_ingv=4:4+11;


BASIN={'adr1','adr2','aeg','lev2','lev3','lev4','tyr1','tyr2','nwm','dard'};

%% row with nan are equally subdiveded along the river
for im=1:12
RUNOFF_INGV(1:2,3+im)=RUNOFF_INGV(1,3+im)/2; % Nile (2 points)
RUNOFF_INGV(18:19,3+im)=RUNOFF_INGV(18,3+im)/2; % Buna        2 points
RUNOFF_INGV(20:22,3+im)=RUNOFF_INGV(20,3+im)/3; % Dardanelles 3 points
RUNOFF_INGV(36:37,3+im)=RUNOFF_INGV(36,3+im)/2; % Tevere      2 points  
end
RUNOFF_INGV_TOT=mean(RUNOFF_INGV(:,4:15),2);
NF_INGV=55;
INGV2OGS=[40,40,24,8,8,8,8,8,8,8,8,8,8,11,11,26,25,21,21,41,41,41,1,2,3,4,5,6,7,9,10,12,13,14,15,16,16,17,18,19,20,22,23,27,29,30,31,32,33,34,35,36,37,38,39];

% read excel file extract from BAU 4.6

OGSmonth=xlsread('input_obc_ALKcorretta.xlsx','monthly');

% name of hte river in the BAU 4.6 file 
NOMI_OGS={'PIAVE','TAGLIAMENTO','SOCA/ISONZO','LIVENZA','BRENTA-BACCHIGLIONE','ADIGE','LIKA','PO','RENO','KRKA','RHONE','ARNO','NERETVA','AUDE','TREBISJNICA','TEVERE/TIBER','MATI','VOLTURNO','SHKUMBINI','STRUMA/STRYMONAS',...
    'BUNA-DRINI','MERIC/EVROS/MARITSA','AXIOS/VARDAR','EBRO','SEMAN','VJOSE/AOOS','ARACHTOS','SUSURLUK/SIMAV/KOCASU','PINIOS','ACHELOOS','GEDIZ','BUYUK MENDERES','KOPRU','MANAVGAT','SEYHAN','CEYHAN','GOKSU','MEDJERDA','ASI/ORONTES','NILE','DARDANELLES'};

NF_OGS=41;
% one OGS river ( SUSURLUK/SIMAV/KOCASU n=28 ) is in the Marmara sea, therefore it is not considered
% therefore, there are 39 rivers + dardanelles
% build index from ingv to ogs
for i=1:1:NF_OGS
ii=find(INGV2OGS==i);
OGS2INGV{i}=ii;
end
%verify indexes are correct :: name of rivers should be aligned
for i=1:1:NF_OGS
    disp(' ')
  disp(NOMI_OGS{i})
  ii=OGS2INGV{i};
  for j=1:1:length(ii)
      disp(NOMI_INGV{ii(j)});
  end
end

% read the excel file containing the concentration values of O3h, O3c and
% O2o to be used to calculated the loads from rivers.
% this file come after the doc file: source_of_info_for_boundary_condition.doc
% with a few modification after the eas2_v12_8 run

O3hO3cO2o=xlsread('Concentrazione_O3h_O3c_O2o_fiumi.xlsx','Sheet1');
% 55 lines: rivers already in the INGV order
% correct values of O3h are for column 5   mmol/m3
% correct values of O3c are for column 8   mgC/m3
% correct values of Oo2 are for column 10  mmol/m3


% START PROCEDURE FOR BUILDING THE EXCEL RIVER FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILD FOGLIO "monthly"
monthly(1:NF_INGV+1,1:24)=-999;
monthly(1:NF_INGV,1)=[1:1:NF_INGV];                  % ID
monthly(1:NF_INGV,4)=RUNOFF_INGV(1:NF_INGV,ilon)-222; % i
monthly(1:NF_INGV,5)=RUNOFF_INGV(1:NF_INGV,ilat);     % j
for i=1:1:NF_INGV
monthly(i,9)= sum(INGV2OGS==INGV2OGS(i));   % OBS
for im=1:1:12
    monthly(i,9+im)=RUNOFF_INGV(i,3+im)/RUNOFF_INGV_TOT(i)/12; 
end
monthly(i,22)= 100;   % 100
monthly(i,23)= RUNOFF_INGV_TOT(i);
end
[status,message] = xlswrite('monthly.xls',monthly,'monthly')


%%% For rivers with more than 1 mouths (Po,
%%% dardanelles...), a factor is computed based on annual total run off for
%%% each of the INGV rivers.

for i=1:1:NF_INGV
ii=find(INGV2OGS==INGV2OGS(i));
FRAZ(i)=RUNOFF_INGV_TOT(i)/sum(RUNOFF_INGV_TOT(ii));
end

asse_anni=[1980:2020];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list of excel sheet to be created
list={'DIP_KTperYR_NOBLS','DIN_KTperYR_NOBLS','DIS_KTperYR_NOBLS', ... 
    'KM3perYR_NOBLS', 'm3_s_NOBLS', ...
    'ALK_GmolperYR_NOBLS','DIC_KTperYR_NOBLS','O2O_GmolperYR_NOBLS'};

%%% BUILD SHEETs for the list of variables
for ilist=1:1:length(list)
if ilist<=5 % for DIP,DIN,DIS,KM3perY and m3_s    FROM BAUMAN data
disp([ list{ilist} ' --> from input_obc_ALKcorretta.xlsx'])
ORI=xlsread('input_obc_ALKcorretta.xlsx',list{ilist});
ORI(1,:)=[]; % first row is deleted to be consistent with the index ingv 
NEW=[];
NEW(1:NF_INGV+1,1:51)=-999;
NEW(1:NF_INGV,1)=[1:1:NF_INGV];  
NEW(1:NF_INGV,4)=RUNOFF_INGV(1:NF_INGV,ilon)-222; % i
NEW(1:NF_INGV,5)=RUNOFF_INGV(1:NF_INGV,ilat);     % j
for i=1:1:NF_INGV
for ia=1:41
    NEW(i,9+ia)=ORI(INGV2OGS(i),9+ia)*FRAZ(i);
end
end

elseif ilist==6  % ALK_GmolperYR_NOBLS
    disp([list{ilist}  ' --> conc from new excel table'])
NEW=[]; NEW(1:NF_INGV+1,1:51)=-999;
NEW(1:NF_INGV,1)=[1:1:NF_INGV];  
NEW(1:NF_INGV,4)=RUNOFF_INGV(1:NF_INGV,ilon)-222; % i
NEW(1:NF_INGV,5)=RUNOFF_INGV(1:NF_INGV,ilat);     % j
for i=1:1:NF_INGV
for ia=1:41
    NEW(i,9+ia)=RUNOFF_INGV_TOT(i)*O3hO3cO2o(i,5)* 86400*365/1000/1000/1000/1000 ; %[m3/s] [mmol O3h / m3]  * 86400 *365 / 1000 /1000 / 1000 /1000 --> Gmol/Y
end
end
elseif ilist==7  % DIC_KTperYR_NOBLS 
    disp([list{ilist}  ' --> conc from new excel table'])
NEW=[]; NEW(1:NF_INGV+1,1:51)=-999;
NEW(1:NF_INGV,1)=[1:1:NF_INGV];  
NEW(1:NF_INGV,4)=RUNOFF_INGV(1:NF_INGV,ilon)-222; % i
NEW(1:NF_INGV,5)=RUNOFF_INGV(1:NF_INGV,ilat);     % j
for i=1:1:NF_INGV
for ia=1:41
    NEW(i,9+ia)=RUNOFF_INGV_TOT(i)*O3hO3cO2o(i,8)* 86400*365/1000/1000/1000/1000 ; %[m3/s] [mg O3C / m3]  * 86400 *365 / 1000 /1000 / 1000 /1000 --> KT O3C/Y
end
end
elseif ilist==8  % O2O_GmolperYR_NOBLS
    disp([list{ilist}  ' --> conc from new excel table'])
NEW=[]; NEW(1:NF_INGV+1,1:51)=-999;
NEW(1:NF_INGV,1)=[1:1:NF_INGV];  
NEW(1:NF_INGV,4)=RUNOFF_INGV(1:NF_INGV,ilon)-222; % i
NEW(1:NF_INGV,5)=RUNOFF_INGV(1:NF_INGV,ilat);     % j
for i=1:1:NF_INGV
for ia=1:41
    NEW(i,9+ia)=RUNOFF_INGV_TOT(i)*O3hO3cO2o(i,10)* 86400*365/1000/1000/1000/1000 ; %[m3/s] [mmol Oo2 / m3]  * 86400 *365 / 1000 /1000 / 1000 /1000 --> Gmol O2/m3
end
end

end

% check that 
ia=5;
  for iogs=1:1:NF_OGS
%      disp([NOMI_OGS{iogs} ' ' num2str(ORI(iogs,9+ia)) ' --> ' num2str(sum(NEW(OGS2INGV{iogs},9+ia)))]);
  end

  %
% comute mean 2000-2010 discharges and mean concentratino and runoff
% considering INGV river

  
for ib=1:1:length(BASIN)
ind=strmatch(BASIN{ib},BAS_INGV);

BAS_runoff(ib,1)=sum(RUNOFF_INGV_TOT(ind));
% average discharges for the period 2000-2010 --> index of NEW matrix 9+21:9+30 
DISCH=[];
for i=1:1:length(ind)
    DISCH(i)=mean(NEW(ind(i),9+21:9+30));
end
BAS_discharge(ib,ilist)=sum(DISCH)';
end


% calcolo la concentrazione media
if strmatch('DIP_KTperYR_NOBLS',list{ilist})
FACTCONV= 1000*1000*1000*1000/31 / (365*86400)   % KT DIP /Y / (m3/s)  * (1000 *1000 * 1000*1000 /31) / ( 365 *86400) --> mmol/m3
elseif strmatch('DIN_KTperYR_NOBLS',list{ilist})
FACTCONV= 1000*1000*1000*1000/16 / (365*86400)   % KT DIN /Y / (m3/s)  * (1000 *1000 * 1000*1000 /16) / ( 365 *86400) --> mmol/m3
elseif strmatch('DIS_KTperYR_NOBLS',list{ilist})
FACTCONV= 1000*1000*1000*1000/14 / (365*86400)   % KT DIS /Y / (m3/s)  * (1000 *1000 * 1000 *1000/14) / ( 365 *86400) --> mmol/m3
elseif strmatch('DIC_KTperYR_NOBLS',list{ilist})
FACTCONV= 1000*1000*1000*1000 / (365*86400)   % KT DIC /Y / (m3/s)  * (1000 *1000 * 1000 *1000) / ( 365 *86400) --> mgC/m3
elseif strmatch('ALK_GmolperYR_NOBLS',list{ilist})
FACTCONV= 1000*1000*1000*1000 / (365*86400)   % GIGA mol alk /Y / (m3/s)  * (1000 *1000 * 1000 *1000) / ( 365 *86400) --> mmol/m3
elseif strmatch('O2o_GmolperYR_NOBLS',list{ilist})
FACTCONV= 1000*1000*1000*1000 / (365*86400)   % GIFG mol O2 /Y / (m3/s)  * (1000 *1000 * 1000 *1000) / ( 365 *86400) --> mmol/m3
else
    FACTVONV=1
end
BAS_conc(:,ilist)=BAS_discharge(:,ilist)./BAS_runoff(:)*FACTCONV


%disp(BASIN')
%disp(num2str(BAS_runoff'))
%disp(num2str(BAS_discharge'))
  
  
  
  

% PRINT EXCEL PROCEDURE ::: add 1 row of 0 that is the one with the header that was prior cancelled 
NEW1(2:NF_INGV+2,1:51)=NEW(1:NF_INGV+1,1:51);
[status,message] = xlswrite(list{ilist},NEW1,list{ilist});
end


% save file of conc discharge and runoff
[status,message] = xlswrite('BAS_runoff.xls', BAS_runoff);
[status,message] = xlswrite('BAS_conc.xls', BAS_conc);
[status,message] = xlswrite('BAS_discharge.xls', BAS_discharge);

