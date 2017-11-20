function [mdata mdata_s]=read_calmo_sim(vars,sims,date_lim,sound_exist,size_vars_sound)

% NAME
%   read_calmo_sim
% PURPOSE
%   Read simulations data from, for specified period
% INPUTS
%   vars - calibrated fields groups. See namelist.m
%   sims - simulations names to be read
%   date_lim - Structure which includes the dates range of simulations to be read
%   sound_exist - binary matrix [Day,Hour,Sounding location] with ones where the sounding data exist
%   size_vars_sound - length(vars_sound) - number of soundings fields
% OUTPUTS
%   mdata - Data matrix with dimensions [Field,Day,simulation,Lon,Lat] (field can be Tmax,Tmin,Pr)
%   mdatas - Data matrix with dimensions [Field,Day,simulation,Hour,Sounding location]
% AUTHOR
%   Pavel Khain (pavelkh_il@yahoo.com)


global maindir curdir simuldir extdir;    % the pathes used (chosen in ReadData_and_MetaModel.m)

mdata=[];       % regular mod data
mdata_s=[];   % soundings mod data

%--------------------------------------------------------------------
% DEFINITIONS of variable metadata
%--------------------------------------------------------------------

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

%--------------------------------------------------------------------
% CONTROL PARAMETERS :
%--------------------------------------------------------------------
soundingcord=load([extdir '/radiosondes_metadata.txt']);    % 'soundingcord' = stations number, lat and lon
mdata_s=zeros(size_vars_sound,totaldays,length(sims),8,size(soundingcord,1))-999.9;
height_step=10; % The interpolation height step in meters. Be careful of decided to we enlarge this value
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------
% READ mod data
%--------------------------------------------------------------------
for n=1:length(sims)
    display(['Reading data of simulation ' sims{n}])
    for i=1:length(vars)
        simdir=[simuldir '/LF_RING_' char(sims{n})]; % Path to the file

        indd=0;
        for m=month_min:month_max % loop over months
            if m<10, ms=['0' num2str(m)]; else ms=num2str(m); end
            if m<month_max, day_max_tmp=dinm(m); else day_max_tmp=day_max; end
            if m>month_min, day_min_tmp=1; else day_min_tmp=day_min; end

            for d=day_min_tmp:day_max_tmp
                if d<10, ds=['0' num2str(d)]; else ds=num2str(d); end
                
                indd=indd+1; %index of the day
                ndate=['13' ms ds '00'];         
                if strcmp(vars{i},'pr')||strcmp(vars{i},'t2m_max')||strcmp(vars{i},'t2m_min')||strcmp(vars{i},'t2m_avg')
                    fname=['aggregated_' char(sims{n}) '_20' ndate '.nc'];
                    fullfname=[simdir '/' fname];
                    MyFolderInfo = dir(fullfname);
                    if exist(fullfname,'file') && MyFolderInfo.bytes > 20000000
                        display(['reading the forecasting file ' fullfname])
                        ncid = netcdf.open(fullfname,'NC_NOWRITE');
                        varid = netcdf.inqVarID(ncid,varname{i});
                        data  = netcdf.getVar(ncid,varid);
                        m_fv=netcdf.getAtt(ncid,varid,'_FillValue');
                        data(data==m_fv)=NaN;
                        data=data*mfact1(i) + mfact2(i);
                        
                        if strcmp(vars{i},'pr')
                            data(data<0 | data>300)=NaN;
                        elseif strcmp(vars{i},'t2m_max')||strcmp(vars{i},'t2m_min')||strcmp(vars{i},'t2m_avg')
                            data(data<-30 | data>50)=NaN;
                        end
                        mdata(i,indd,n,:,:)=data;
                        netcdf.close(ncid);
                    else
                        display(['missing ' fullfname]);
                        %%%%%%%%%%% Define parameters for the numbers: %%%%%%%%%%%%
                        mdata(i,indd,n,:,:)=NaN(1158,774);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    
                elseif strcmp(vars{i},'sound')
                    display([yearlimc ms ds]);            
                    size(mdata_s)
                    display(['the ',num2str(n),' and d ',num2str(d),' and yearlim ',num2str(yearlim),' and m ',num2str(m)])
                    mdata_s(1:size_vars_sound,indd,n,:,:)=read_profiles_mod(simdir,yearlim,m,d,sims{n},height_step,squeeze(sound_exist(indd,:,:)),windshear,windshearopt);
                else
                    display(['Wrong field name' vars{i}]);
                end
            end % days
        end     % months
        %%%%%%%%%%%%%%%%% Check if needed: %%%%%%%%%%%%%%%%%%%%
        if isempty(mdata_s)
            check_break=~isnan(squeeze(mdata(i,:,n,:,:)));
            if sum(check_break(:))==0
                display(['WARNING: MISSING DATA FOR simulation=' num2str(i) ' var=' num2str(n)]);
                return
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % variables
end % sims




