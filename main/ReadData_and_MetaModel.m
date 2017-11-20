function ReadData_and_MetaModel(date_min,date_max)

% NAME 
%   ReadData_and_MetaModel
% PURPOSE 
%   Read observations and simulations data, then fit the Meta-Models
% INPUTS 
%   time period: from date_min 'dd-mmm-yyyy' to date_max 'dd-mmm-yyyy'
% OUTPUTS 
%   saved (in .mat format) observations and simulations fields as well as
%   the Meta-Models coefficients
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)


date_lim=struct('dmax',date_max,'dmin',date_min);
tmp_str=datestr(date_min);
mkdir(tmp_str);

%------------------------------------------------------------
%-1- namelist:
%------------------------------------------------------------
%clear global comp maindir curdir simuldir;
global comp maindir curdir simuldir obsdir extdir;

[maindir simuldir obsdir extdir vars vars_2d avg_T vars_sound sims_opt ml score w_user lhacc iterations_num best_percent date_min_tot date_max_tot]=namelist();
%------------------------------------------------------------
%-4- Define "parameters" structure
%------------------------------------------------------------

[paramn,paramnt,range,default,expval,valval,simval,sims_reg,sims_inter,sims_con,valcon,param_log,date_min_check,date_max_check]=sims_def(sims_opt);
if datenum(date_min)<datenum(date_min_check) || datenum(date_max)>datenum(date_max_check)
    stop
end
parameters=struct('name',paramn,'range',range','default',default,'experiments', ...
    expval,'constrain',valcon,'validation',valval,'name_tex',paramnt);
    
    
    %------------------------------------------------------------
    %-5- Define "datamatrix" structure
    %------------------------------------------------------------
    
    % (1) Read observations
    display('Read obs')
    [datamatrix.obsdata datamatrix_s.obsdata datamatrix_s.sound_exist]=read_calmo_obs(vars,date_lim,'clever_mean',length(vars_sound));
    
    % (2) Read reference simulation
    sims={'DEF'};
    [datamatrix.refdata datamatrix_s.refdata]=read_calmo_sim(vars,sims,date_lim,datamatrix_s.sound_exist,length(vars_sound));
    
    % (3) Read parameter experiments
    [datamatrix.moddata datamatrix_s.moddata]=read_calmo_sim(vars,sims_reg,date_lim,datamatrix_s.sound_exist,length(vars_sound));
    %size: 3-CAPE,SIN,TCWC    29-number of days    st 8- 3h steps in a day  59-number of r.s. station

    % (4) Read interaction parameter experiments
    par_num=length(paramn);
    if ~isempty(sims_inter) % IN CASE THERE ARE INTERACTION TERMS:
        [datamatrix.moddata(:,:,2*par_num+1:2*par_num+length(sims_inter),:,:) datamatrix_s.moddata(:,:,2*par_num+1:2*par_num+length(sims_inter),:,:)]=read_calmo_sim(vars,sims_inter,date_lim,datamatrix_s.sound_exist,length(vars_sound));
    end
    
    % (5) Read validation experiments
    %if ~strcmp(simval,'MISSING')
    %    datamatrix.valdata=read_calmo_sim(vars,simval,date_lim);
    %end
    
    % (6) Read constrain experiments
    if ~strcmp(sims_con,'MISSING') % IN CASE THERE ARE CONSTRAIN (INTERMEDIATE) SIMULATIONS:
        [datamatrix.constrain datamatrix_s.constrain]=read_calmo_sim(vars,sims_con,date_lim,datamatrix_s.sound_exist,length(vars_sound));
    end
    
    if ~isempty(strfind([vars{:}],'pr')) || ~isempty(strfind([vars{:}],'t2m_avg')) || ~isempty(strfind([vars{:}],'t2m_max')) || ~isempty(strfind([vars{:}],'t2m_min'))
        
        % if the surf. obs. number for given field,i,j (usually per month) is <3, cancel the entire month
        for f=1:size(datamatrix.obsdata,1)
            for i=1:size(datamatrix.obsdata,4)
                for j=1:size(datamatrix.obsdata,5)
                    if sum(~isnan(squeeze(datamatrix.obsdata(f,:,1,i,j))))<3
                        datamatrix.obsdata(f,1:end,1,i,j)=NaN;
                    end
                end
            end
        end
       
        % cancel all the missing Italian grid-points:
        [a]=isnan(datamatrix.obsdata);
        datamatrix.refdata(a==1)=NaN;
        if isfield(datamatrix, 'valdata')
            datamatrix.valdata(a==1)=NaN;
        end
        for i=1:size(datamatrix.moddata,3)
            clear tmp;
            tmp=datamatrix.moddata(:,:,i,:,:);
            tmp(a==1)=NaN;
            datamatrix.moddata(:,:,i,:,:)=tmp;
        end
        if isfield(datamatrix, 'constrain')
            for i=1:size(datamatrix.constrain,3)
                clear tmp;
                tmp=datamatrix.constrain(:,:,i,:,:);
                tmp(a==1)=NaN;
                datamatrix.constrain(:,:,i,:,:)=tmp;
            end
        end
        
        
        % get domain lon/lat:  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        file=[extdir '/aggregated_LTUR_2013011000.nc'];
        ncid = netcdf.open(file,'NC_NOWRITE');
        varid_lon = netcdf.inqVarID(ncid,'lon_1');
        varid_lat = netcdf.inqVarID(ncid,'lat_1');
        lat  = netcdf.getVar(ncid,varid_lat);
        lon  = netcdf.getVar(ncid,varid_lon);
        lon(lon>180)=lon(lon>180)-360;
        netcdf.close(ncid);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %--------------------
        display('Redefine datamatrix according C. Frei regions division of Switzerland');
        %--------------------
        unify_regions=[3 5];    % which regions to unify
        [datamatrix_new] = Frei_regions(datamatrix,lat,lon,vars,avg_T,unify_regions);
        
        
        if avg_T==0
            datamatrix_pr=datamatrix_new;
            clear datamatrix_new
            if ~isempty(strfind([vars{:}],'t2m_max'))
                datamatrix_tmax=gpts_series(datamatrix,vars,'t2m_max');
            end
            if ~isempty(strfind([vars{:}],'t2m_min'))
                datamatrix_tmin=gpts_series(datamatrix,vars,'t2m_min');
            end
        end
    end
    
    
    if ~isempty(strfind([vars{:}],'sound'))
        datamatrix_so=gpts_series(datamatrix_s,vars,'sound');
        
        %%% Assign NaN in Sonde fields that have less than ml observations per month or period. The Assigment is for all days in the specific month and reigons for that specific month and fields. 
        datatemp=del_bad_sounding(datamatrix_so,ml,sims_reg,sims_inter,sims_con,simval,vars_sound);
        clear datamatrix_so
        datamatrix_so=datatemp;
        clear datatemp
	end
   
    
    %%%%%%%%%%%%%%%%%%%%%% SAVE Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%
	display('saving matrices');
    if ~isempty(strfind([vars{:}],'pr')) || ~isempty(strfind([vars{:}],'t2m_avg')) || ~isempty(strfind([vars{:}],'t2m_max')) || ~isempty(strfind([vars{:}],'t2m_min'))
        clear tmptmp;
        tmptmp=datamatrix.moddata;
        save([tmp_str '/moddata.mat'],'tmptmp','-v7.3')
        save([tmp_str '/datamatrix.mat'],'-struct','datamatrix');
        
        if avg_T==0
            save([tmp_str '/datamatrix_pr.mat'],'-struct','datamatrix_pr');
            save([tmp_str '/datamatrix_tmax.mat'],'-struct','datamatrix_tmax');
            save([tmp_str '/datamatrix_tmin.mat'],'-struct','datamatrix_tmin');
        else
            save([tmp_str '/datamatrix_new.mat'],'-struct','datamatrix_new');
        end
    end
        
    if ~isempty(strfind([vars{:}],'sound'))
        clear tmptmp;
        tmptmp=datamatrix_s.moddata;
        save([tmp_str '/moddata_s.mat'],'tmptmp','-v7.3')
        save([tmp_str '/datamatrix_s.mat'],'-struct','datamatrix_s');
        save([tmp_str '/datamatrix_so.mat'],'-struct','datamatrix_so'); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % FIT META-MODELS :
    if avg_T==0
        if ~isempty(strfind([vars{:}],'t2m_max'))
            metamodel_tmax=neelin_e(parameters,datamatrix_tmax,{'t2m_max'}); % structure which contains the fitted metamodel parameters
            save([tmp_str '/metamodel_tmax.mat'],'-struct','metamodel_tmax');
        end
        if ~isempty(strfind([vars{:}],'t2m_min'))
            metamodel_tmin=neelin_e(parameters,datamatrix_tmin,{'t2m_min'}); % structure which contains the fitted metamodel parameters
            save([tmp_str '/metamodel_tmin.mat'],'-struct','metamodel_tmin');
        end
        if ~isempty(strfind([vars{:}],'pr'))
            metamodel_pr=neelin_e(parameters,datamatrix_pr,{'pr'}); % structure which contains the fitted metamodel parameters
            save([tmp_str '/metamodel_pr.mat'],'-struct','metamodel_pr');
        end
    else
        metamodel_new=neelin_e(parameters,datamatrix_new,vars_2d); % structure which contains the fitted metamodel parameters
        save([tmp_str '/metamodel_new.mat'],'-struct','metamodel_new');
    end
    if ~isempty(strfind([vars{:}],'sound'))
        metamodel_so=neelin_e(parameters,datamatrix_so,vars_sound); % structure which contains the fitted metamodel parameters for 18 parmeters including T850,T700mb,T500mb,
        % RH850mb, RH700mb, RH500mb , U850mb, U700mb, U500mb , V850mb, V700mb and V500mb.
        save([tmp_str '/metamodel_so.mat'],'-struct','metamodel_so');
    end
        
