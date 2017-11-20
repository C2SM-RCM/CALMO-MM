function [odata odata_s sound_exist]=read_calmo_obs(vars,date_lim,avg_type,size_vars_sound)

% NAME
%   read_calmo_obs
% PURPOSE
%   Read observations data from Switzerland and north Italy, as well as soundings data, for specified period
% INPUTS
%   vars - calibrated fields groups. See namelist.m
%   date_lim - Structure which includes the dates range of simulations to be read
%   avg_type - Interpolation method of observations data to model grid. Can be: 'near_neighb','simple_mean','weight_mean','clever_mean'
%   size_vars_sound - length(vars_sound) - number of soundings fields
% OUTPUTS
%   odata - Data matrix with dimensions [Field,Day,1,Lon,Lat] (field can be Tmax,Tmin,Pr)
%   odatas - Data matrix for profiles with dimensions [Field,Day,1,Hour,Sounding location]
%   sound_exist - binary matrix [Day,Hour,Sounding location] with ones where the sounding data exist
% AUTHOR
%   Pavel Khain (pavelkh_il@yahoo.com)

global obsdir extdir;    % the pathes used (chosen in ReadData_and_MetaModel.m)
odata=[];    % total regular obs data
odata_s=[]; % soundings data

%--------------------------------------------------------------------
% DEFINITIONS of variable metadata
%--------------------------------------------------------------------

if strcmp(avg_type,'clever_mean')
    filename_Tmax='obs_Tmax_clever'; filename_Tmin='obs_Tmin_clever'; 
elseif strcmp(avg_type,'near_neighb')
	filename_Tmax='obs_Tmax_neighb'; filename_Tmin='obs_Tmin_neighb'; 
elseif strcmp(avg_type,'simple_mean')
	filename_Tmax='obs_Tmax_simple'; filename_Tmin='obs_Tmin_simple'; 
elseif strcmp(avg_type,'weight_mean')
	filename_Tmax='obs_Tmax_weight'; filename_Tmin='obs_Tmin_weight'; 
else
    display(['Wrong avg_type, exitting.']);
    return
end
display('Checking Swiss obs')   
if ~exist([obsdir '/obs_switzerland_1km/temperature/' filename_Tmax '.mat'],'file')
    display('Reading original Swiss Tmax data'); 
    build_temp_obs(avg_type,'Tmax'); 
end
if ~exist([obsdir '/obs_switzerland_1km/temperature/' filename_Tmin '.mat'],'file')
    display('Reading original Swiss Tmin data');
    build_temp_obs(avg_type,'Tmin'); 
end
if ~exist([obsdir '/obs_switzerland_1km/rain/obs_rain_neighb.mat'],'file')
    display('Reading original Swiss rain data');
    build_rain_obs(); 
end
display('Swiss obs files exist. Go on.')    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% You can look at the corrected observations temperature fields by:
% figure(1); filename_Tmax='obs_Tmax_clever'; load([obsdir '/obs_switzerland_1km/temperature/' filename_Tmax '.mat']); surf(obs_T_fixed(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); caxis([-20 10]); colorbar;
% figure(2); filename_Tmax='obs_Tmax_neighb'; load([obsdir '/obs_switzerland_1km/temperature/' filename_Tmax '.mat']); surf(obs_T_fixed(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); caxis([-20 10]); colorbar;
% figure(3); filename_Tmax='obs_Tmax_simple'; load([obsdir '/obs_switzerland_1km/temperature/' filename_Tmax '.mat']); surf(obs_T_fixed(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); caxis([-20 10]); colorbar;
% figure(4); filename_Tmax='obs_Tmax_weight'; load([obsdir '/obs_switzerland_1km/temperature/' filename_Tmax '.mat']); surf(obs_T_fixed(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); caxis([-20 10]); colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% CORRECT: %%%%%%%%%%%%%%%%%%%%%%%%%%%
[varname,ofact1,ofact2,mfact1,mfact2,unit]=var_meta_calmo(vars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------
% DEFINITIONS of date constants
%--------------------------------------------------------------------
dinm=[31 28 31 30 31 30 31 31 30 31 30 31];
yearlimc=date_lim.dmax(8:11); % the year in char variable
yearlim=str2num(yearlimc); yearesidual=yearlim-floor(yearlim/4)*4; if (yearesidual==0) dinm(2)=29; end 
date_max_num=datenum(date_lim.dmax, 'dd-mmm-yyyy'); date_max_vec=datevec(date_max_num);
date_min_num=datenum(date_lim.dmin, 'dd-mmm-yyyy'); date_min_vec=datevec(date_min_num);
day_max=date_max_vec(3); day_min=date_min_vec(3);
month_max=date_max_vec(2); month_min=date_min_vec(2);
totaldays=date_max_num-date_min_num+1;
%%%%%%%%%%% CORRECT: %%%%%%%%%%%%%%%%%%%%%%%%%%%
days_from_1jan=date_min_num-datenum('01-Jan-2013', 'dd-mmm-yyyy');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------
% CONTROL PARAMETERS :
%--------------------------------------------------------------------
Pmax=200; % Limit for precipitaton intensity values for Italian stations: the max rain intensity in mm per day
Tmaxabs=50; % Limit for temperature values for Italian stations: the absloute maximum temperature in Celsius
Tminabs=-50; % Limit for temperature values for Italian stations: the absloute mainimum temperaturethe in Celsius
soundingcord=load([extdir '/radiosondes_metadata.txt']);    % 'soundingcord' = stations number, lat and lon
odata_s=zeros(size_vars_sound,totaldays,1,8,size(soundingcord,1))-999.9;
sound_exist=zeros(totaldays,8,size(soundingcord,1));
height_step=10; % The interpolation height step in meters. Be careful if decided to enlarge this value
v1u=500; % the upper pressure level for wshear1
v1d=700; % the bottom pressure level for wshear 1
v2u=700; % as mentioned above but for wshear 2
v2d=850; % as mentioned abobe but for wshear 2
v3u=850; % as mentioned above but for wshear 3
v3d=1100; % as mentioned above but for wshear 3, whereas 1100 is the surface level or the lowest level ("below surface")
windshear=[v1u,v1d,v2u,v2d,v3u,v3d]; % wind shear between upper pressure level up ("u") to bottom pressure level ("d") for each of the WSHEAR(1,2,3), ex: if v1u=500 and v1d=700 than the WSHEAR1 is between 500mb to 700mb'
%Choose:
%windshearopt='scalar'; % windshear calculation
windshearopt='vector';    % windshear calculation

%--------------------------------------------------------------------
% READ observations data
%--------------------------------------------------------------------

%--------------------
display('PART 1 - SWISS DATA');
%--------------------
for i=1:length(vars)
    clear data;
	if strcmp(vars{i},'t2m_avg')
        display([' Reading Swiss observations of Tavg']);
        load([obsdir '/obs_switzerland_1km/temperature/' filename_Tavg '.mat']);
        data=obs_T_fixed;
    elseif strcmp(vars{i},'t2m_max')
        display([' Reading Swiss observations of Tmax']);
        load([obsdir '/obs_switzerland_1km/temperature/' filename_Tmax '.mat']);
        data=obs_T_fixed;
    elseif strcmp(vars{i},'t2m_min')
        display([' Reading Swiss observations of Tmin']);
        load([obsdir '/obs_switzerland_1km/temperature/' filename_Tmin '.mat']);
        data=obs_T_fixed;
    elseif strcmp(vars{i},'pr')
        display([' Reading Swiss precipitation']);
        load([obsdir '/obs_switzerland_1km/rain/obs_rain_neighb.mat']);
        data=obs_rain_fixed;
    elseif strcmp(vars{i},'sound') 
        display([' Reading soundings']);
    end
    
    indd=0;
    for m=month_min:month_max % loop over months
        if m<10, ms=['0' num2str(m)]; else ms=num2str(m); end
        if m<month_max, day_max_tmp=dinm(m); else day_max_tmp=day_max; end
        if m>month_min, day_min_tmp=1; else day_min_tmp=day_min; end
        
        for d=day_min_tmp:day_max_tmp % loop over days
            if d<10, ds=['0' num2str(d)]; else ds=num2str(d); end
            indd=indd+1; %index of the day
    
            if strcmp(vars{i},'pr')||strcmp(vars{i},'t2m_max')||strcmp(vars{i},'t2m_min')||strcmp(vars{i},'t2m_avg')
                %%%%%%%%%%% CORRECT: %%%%%%%%%%%%%%%%%%%%%%%%%%%
                odata(i,indd,1,:,:)=squeeze(data(:,:,days_from_1jan+indd));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif(strcmp(vars{i},'sound'))
                display([yearlimc ms ds]);
                clear tmp1;
                tmp1=read_sounding_obs(yearlimc,m,d,height_step,windshear,windshearopt);
                odata_s(1:size_vars_sound,indd,1,:,:)=tmp1;
                tmp2=squeeze(tmp1(3,:,:)); tmp2(tmp2~=-999.9)=1; tmp2(tmp2==-999.9)=0;
                sound_exist(indd,:,:)=tmp2;  
            else
                display(['Wrong field name' vars{i}]);
            end % variables
        end % if month
    end % if var
end % if day

%----------------------
display('PART 2 - ITALIAN DATA');
%----------------------
for i=1:length(vars)
    clear data;
    if ~strcmp(vars{i},'sound')
        if strcmp(vars{i},'t2m_avg'), display([' Reading Italian t2m_avg']);
        elseif strcmp(vars{i},'t2m_max'), display([' Reading Italian t2m_max']);
        elseif strcmp(vars{i},'t2m_min'), display([' Reading Italian t2m_min']);
        elseif strcmp(vars{i},'pr'), display([' Reading Italian precipitation']);
        end

        indd=0;
        for m=month_min:month_max % loop over months
            if m<10, ms=['0' num2str(m)]; else ms=num2str(m); end
            if m<month_max, day_max_tmp=dinm(m); else day_max_tmp=day_max; end
            if m>month_min, day_min_tmp=1; else day_min_tmp=day_min; end

            for d=day_min_tmp:day_max_tmp % loop over days
                if d<10, ds=['0' num2str(d)]; else ds=num2str(d); end
                yy_prev=datestr(datenum([ds '-' ms '-2013'], 'dd-mm-yyyy')+1,'yy');
                ms_prev=datestr(datenum([ds '-' ms '-2013'], 'dd-mm-yyyy')+1,'mm');
                ds_prev=datestr(datenum([ds '-' ms '-2013'], 'dd-mm-yyyy')+1,'dd');
                indd=indd+1; %index of the day

                if strcmp(vars{i},'t2m_avg')               varnameitaly='TAVE_2M';   fname=['tavg_1_20' yy_prev ms_prev ds_prev '_0.grib.nc'];   % The name of the field in the Italian NedCDF Data 
                elseif strcmp(vars{i},'t2m_max')     varnameitaly='TMAX_2M';   fname=['tmax_1_20' yy_prev ms_prev ds_prev '_0.grib.nc'];   % The name of the field in the Italian NedCDF Data
                elseif strcmp(vars{i},'t2m_min')      varnameitaly='TMIN_2M';   fname=['tmin_1_20' yy_prev ms_prev ds_prev '_0.grib.nc'];   % The name of the field in the Italian NedCDF Data 
                elseif strcmp(vars{i},'pr')                   varnameitaly='TOT_PREC';   fname=['tp_1_20' yy_prev ms_prev ds_prev '_0.grib.nc'];   % The name of the field in the Italian NedCDF Data     
                end

                if exist(strcat(obsdir,'/obs_Italy_1km/st_data_nc/',fname),'file')
                    strcat(obsdir,'/obs_Italy_1km/st_data_nc/',fname)
                    ncid = netcdf.open([obsdir '/obs_Italy_1km/st_data_nc/' fname],'NC_NOWRITE');
                    varid = netcdf.inqVarID(ncid,varnameitaly);
                    data = netcdf.getVar(ncid,varid);
                    m_fv=netcdf.getAtt(ncid,varid,'_FillValue');
                    data(data==m_fv)=NaN;
                    %%%data=data*ofact1(i) + ofact2(i);
                    netcdf.close(ncid);
                    M2=data(:); % RAW data from Italy
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    data1=squeeze(odata(i,indd,1,:,:)); % Matrix data from SWISS
                    M1=data1(:); % RAW data from SWISS
                    Mold=M1; % saving backup of SWISS data
                    if strcmp(vars{i},'pr') % checking the identical grid points location in ITALY and SWISS domains.
                        k2=find((M2>=0 & M2<=Pmax) &  ((Mold<0 | Mold>Pmax) | isnan(Mold))); % The avaliable points in italy that indentical to SWISS location grid points and the data there do not exist
                        M1(k2)=M2(k2);  % assign Italy data in SWISS data !
                    else
                        % the temperature data in SWISS NetCDF files are in Celcuis, whereas the temperature data in the Italin NetCDF files are in Kelvin.
                        k2=find((M2>=Tminabs+273.15 & M2<=Tmaxabs+273.15) &  ((Mold<Tminabs | Mold>Tmaxabs) | isnan(Mold))); % The avaliable points in italy that indentical to SWISS location grid points and the datatom there do not exist.                    M1(k2)=M2(k2)-273.15;  % assign Italy data in SWISS data (The italy data are in Kelvin but in SWISS in Celsius
                        M1(k2)=M2(k2)-273.15;  % assign Italy data in SWISS data after transfer from Kelvin to Celcuis degree
                    end
                    data3=reshape(M1,[size(data1,1),size(data1,2)]);
                    odata(i,indd,1,:,:)=data3;      % here odata contains both Swiss and Italian data !
                else
                    display(['missing ' obsdir '/obs_Italy_1km/st_data_nc/' fname]);
                end % if exists
            end % days
        end % months
    end
end % end of variable

%%%%%%%%%%%%%%%%% Check if needed: %%%%%%%%%%%%%%%%%%%%
%size(odata)
%kday_s=size(odata_s,2);
%kday=size(odata,2);
%if(kday-kday_s>0.001)
%    odatatemp=odata;
%    clear mdata
%    AA=size(odatatemp);
%    if(length(AA)==4)
%        odata=odatatemp(:,1:kday_s,:,:);
%    elseif(length(AA)==5)
%        odata=odatatemp(:,1:kday_s,:,:,:);
%    elseif(length(AA)==3)
%        odata=odatatemp(:,1:kday_s,:);
%    end
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







