function build_temp_obs(avg_type,maxminavg)

% NAME 
%   build_temp_obs
% PURPOSE 
%   This function interpolates the gridded temperature observations (by C. Frei) to the model grid. 
% METHOD
%   When comparing smoothed topography model-grid 2m-temperature whith the observed 2m temperature, one shoud "correct" the observed
%   2m-temperature to correspond the model grid elevation. The correction may be performed using the neighboring grid points 2m-temperature profile,
%   according the recommendation of C. Frei
%   Steps:
%   1. Read any simulation file to obtain the model grid lat/lon (ex:aggregated_LTUR_2013011000.nc)
%   2. Read any simulation file to obtain the model grid altitude (ex:laf2013111600_filtered.nc)
%   3. Read any observations file to obtain the observations grid lat/lon (ex: TmaxD_ch01r.swisscors_201301010000_201302010000.nc) 
%   4. Read gridded observations altitude ( ex: topo.swiss1_ch01r.swisscors.nc)
%   5. Interpolate the gridded observations to the model grid using one of the possible methods 
% INPUTS 
%   avg_type - one of the interpolation methods: 'near_neighb','simple_mean','weight_mean','clever_mean'
%   maxminavg - Which field to interpolate: 'Tmax','Tmin','Tavg'
% OUTPUTS 
%   saved (in .mat format) interpolated temperature observations
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)

%clear all; close all;
format long;

global obsdir extdir

% 1. Read any simulation file to obtain the model grid lat/lon (ex:aggregated_LTUR_2013011000.nc)
file=[extdir '/aggregated_LTUR_2013011000.nc'];
ncid = netcdf.open(file,'NC_NOWRITE');
varid_lon = netcdf.inqVarID(ncid,'lon_1');
varid_lat = netcdf.inqVarID(ncid,'lat_1');
%varid_T2m_avg = netcdf.inqVarID(ncid,'T2m');
%varid_T2m_max = netcdf.inqVarID(ncid,'T2m_max');
%varid_T2m_min = netcdf.inqVarID(ncid,'T2m_min');
mod_lat  = netcdf.getVar(ncid,varid_lat);
mod_lon  = netcdf.getVar(ncid,varid_lon);
mod_lon(mod_lon>180)=mod_lon(mod_lon>180)-360;
%mod_T_avg = netcdf.getVar(ncid,varid_T2m_avg);
%mod_T_max = netcdf.getVar(ncid,varid_T2m_max);
%mod_T_min = netcdf.getVar(ncid,varid_T2m_min);
%m_fv=netcdf.getAtt(ncid,varid_T2m_max,'_FillValue');
%mod_T_avg(mod_T_avg==m_fv)=NaN;
%mod_T_max(mod_T_max==m_fv)=NaN;
%mod_T_min(mod_T_min==m_fv)=NaN;
netcdf.close(ncid);
clear file ncid varid_lon varid_lat varid_T2m_avg varid_T2m_max varid_T2m_min m_fv
% at this stage we have: mod_lon mod_lat

% 2. Read any simulation file to obtain the model grid altitude (ex:laf2013111600_filtered.nc)
file=[extdir '/laf2013111600_filtered.nc'];
ncid = netcdf.open(file,'NC_NOWRITE');
varid_h = netcdf.inqVarID(ncid,'HSURF');
mod_H= netcdf.getVar(ncid,varid_h);
netcdf.close(ncid);
clear file ncid varid_h
% at this stage we have: mod_lon mod_lat mod_H

% 3. Read any observations file to obtain the observations grid lat/lon (ex: TmaxD_ch01r.swisscors_201301010000_201302010000.nc) 
file_Tmax=[obsdir '/obs_switzeland_1km/temperature/TmaxD_ch01r.swisscors_201301010000_201302010000.nc'];
%file_Tmin=[obsdir '/obs_switzeland_1km/temperature/TminD_ch01r.swisscors_201301010000_201302010000.nc'];
%file_Tavg=[obsdir '/obs_switzeland_1km/temperature/TavgD_ch01r.swisscors_201301010000_201302010000.nc'];
ncid_Tmax = netcdf.open(file_Tmax,'NC_NOWRITE');
%ncid_Tmin = netcdf.open(file_Tmin,'NC_NOWRITE');
%ncid_Tavg = netcdf.open(file_Tavg,'NC_NOWRITE');
varid_Tmax = netcdf.inqVarID(ncid_Tmax,'TmaxD');  %daily max
%varid_Tmin = netcdf.inqVarID(ncid_Tmin,'TminD');  %daily min
%varid_Tavg = netcdf.inqVarID(ncid_Tavg,'TavgD');  %daily avg
varid_lon = netcdf.inqVarID(ncid_Tmax,'lon'); 
varid_lat = netcdf.inqVarID(ncid_Tmax,'lat'); 
%obs_Tmax = netcdf.getVar(ncid_Tmax,varid_Tmax);
%obs_Tmin = netcdf.getVar(ncid_Tmin,varid_Tmin);
%obs_Tavg = netcdf.getVar(ncid_Tavg,varid_Tavg);
m_fv=netcdf.getAtt(ncid_Tmax,varid_Tmax,'_FillValue');
%obs_Tmax(obs_Tmax==m_fv)=NaN;
%obs_Tmin(obs_Tmin==m_fv)=NaN;
%obs_Tavg(obs_Tavg==m_fv)=NaN;
obs_lon = netcdf.getVar(ncid_Tmax,varid_lon);
obs_lat = netcdf.getVar(ncid_Tmax,varid_lat);
netcdf.close(ncid_Tmax);
%netcdf.close(ncid_Tmin);
clear file_Tmax ncid_Tmax varid_Tmax varid_lon varid_lat 
% at this stage we have: mod_lon mod_lat mod_H obs_lon obs_lat m_fv

% 4. Read gridded observations altitude ( ex: topo.swiss1_ch01r.swisscors.nc)
file=[extdir '/topo.swiss1_ch01r.swisscors.nc'];
ncid = netcdf.open(file,'NC_NOWRITE');
varid_H = netcdf.inqVarID(ncid,'height');
obs_H= netcdf.getVar(ncid,varid_H);
netcdf.close(ncid);
clear file ncid varid_H
% at this stage we have: mod_lon mod_lat mod_H obs_lon obs_lat m_fv obs_H

% 5. Read all the temperature observation files to be interpolated:
MyFolderInfo = dir([obsdir '/obs_switzeland_1km/temperature/' maxminavg 'D_ch01r.swisscors_*.nc']);
for f=1:length(MyFolderInfo)
    %display(num2str(f));
    file=[obsdir '/obs_switzeland_1km/temperature/' MyFolderInfo(f).name];
    ncid = netcdf.open(file,'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,[maxminavg 'D']);  %daily rain
    obs_tmp = netcdf.getVar(ncid,varid);
    obs_tmp(obs_tmp==m_fv)=NaN;
    size(obs_tmp)
    obs_T{f}=obs_tmp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% Can look at the readed data by:
% surf(obs_lon,obs_lat,obs_T{1}(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); colorbar;
% surf(mod_lon,mod_lat,obs_T_fixed(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 5. Interpolate the gridded observations to the model grid using one of the possible methods 

% Definitions:
tot_neigh_num=9;	% how many obs grid neighbors around model grid point to find
tot_neigh_plot=9;	  % how many obs grid neighbors around model grid point to use for building the linear profile
i_size=size(mod_lat,1); j_size=size(mod_lat,2);
% Initializations:
%mean_T(1:i_size,1:j_size,1:size(obs_Tmax,3))=NaN;
%meanD_T(1:i_size,1:j_size,1:size(obs_Tmax,3))=NaN;
%gamma(1:i_size,1:j_size,1:size(obs_Tmax,3))=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% Can look at the interpolated example by defining a grid point and a day:
plot_neigh(1:i_size,1:j_size,1:size(obs_T{1},3))=0; 
plot_neigh(491,412,11)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon_reduc=cos(47*pi/180); % distance in longitude direction getting smaller with latitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find corners of obs grid on the mod grid to save computation time:
dist_lower_left=sqrt((mod_lat-obs_lat(1,1)).^2+(mod_lon-obs_lon(1,1)).^2);
dist_lower_right=sqrt((mod_lat-obs_lat(end,1)).^2+(mod_lon-obs_lon(end,1)).^2);
dist_upper_left=sqrt((mod_lat-obs_lat(1,end)).^2+(mod_lon-obs_lon(1,end)).^2);
dist_upper_right=sqrt((mod_lat-obs_lat(end,end)).^2+(mod_lon-obs_lon(end,end)).^2);
[i_lower_left,j_lower_left]=find(dist_lower_left==min(min(dist_lower_left)));
[i_lower_right,j_lower_right]=find(dist_lower_right==min(min(dist_lower_right)));
[i_upper_left,j_upper_left]=find(dist_upper_left==min(min(dist_upper_left)));
[i_upper_right,j_upper_right]=find(dist_upper_right==min(min(dist_upper_right)));
i_min=min([i_lower_left i_lower_right i_upper_left i_upper_right]);
j_min=min([j_lower_left j_lower_right j_upper_left j_upper_right]);
i_max=max([i_lower_left i_lower_right i_upper_left i_upper_right]);
j_max=max([j_lower_left j_lower_right j_upper_left j_upper_right]);

bad_arr.i=[];
bad_arr.j=[];
bad_arr.dd=[];

% Start interpolation:
for f=1:length(MyFolderInfo)
    obs_T_fixed{f}(1:i_size,1:j_size,1:size(obs_T{f},3))=NaN;
end
for i=(i_min-10):(i_max+10)
    display(['i=' num2str(i)]);
    for j=(j_min-10):(j_max+10)
        %display(['i=' num2str(i) 'j=' num2str(j)]);
        ii=i; jj=j;
        clear i_obs_grid j_obs_grid obs_lon_tmp
        obs_lon_tmp=obs_lon;

        k=1;
        while k<=tot_neigh_num      % loop to find "tot_neigh_num" neighbors of a model g.point
            [i_tmp,j_tmp]=find((lon_reduc.*(obs_lon_tmp-mod_lon(ii,jj))).^2+(obs_lat-mod_lat(ii,jj)).^2==min(min((lon_reduc.*(obs_lon_tmp-mod_lon(ii,jj))).^2+(obs_lat-mod_lat(ii,jj)).^2)));

            % insure that the close grid point we found is inside Switzerland:
            break_flag=0;
            if isnan(obs_T{1}(i_tmp(1),j_tmp(1),1)) && k==1
                break_flag=1;
                break; %the model grid point is too close to Swiss border
            end

            if break_flag==0
                i_obs_grid(ii,jj,k)=i_tmp(1);
                j_obs_grid(ii,jj,k)=j_tmp(1);
                k=k+1;
            end
            %change obs_lon_tmp to be large (1000), so that for next k the g.p. we will find will be no more the nearest:
            obs_lon_tmp(i_tmp(1),j_tmp(1))=1000;     
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % At this stage we have the neighbours i_obs_grid(ii,jj,1:tot_neigh_num),j_obs_grid(ii,jj,1:tot_neigh_num)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for f=1:length(MyFolderInfo)
            if MyFolderInfo(f).bytes < 2000000  || break_flag==1 % too small file size, or to close to border of defined area
                obs_T_fixed{f}(ii,jj,1:size(obs_T{f},3))=NaN;
            else
                for dd=1:size(obs_T{f},3)
                    if (f==1 && plot_neigh(ii,jj,dd)==1)     %plot the obs g.p. neighbors of the mod g.p.
                        for k=1:1
                            plot(obs_T{1}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd),obs_H(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k)),'rx','MarkerSize',15); hold on;
                        end
                        for k=2:4
                            plot(obs_T{1}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd),obs_H(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k)),'bx','MarkerSize',15); hold on;
                        end
                        if tot_neigh_plot>4
                            for k=5:6
                                plot(obs_T{1}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd),obs_H(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k)),'bx','MarkerSize',15); hold on;
                            end
                            for k=7:9
                                plot(obs_T{1}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd),obs_H(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k)),'gx','MarkerSize',15); hold on;
                            end
                            if tot_neigh_plot>9
                                for k=10:tot_neigh_plot
                                    plot(obs_T{1}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd),obs_H(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k)),'yx','MarkerSize',15); hold on;
                                end
                            end
                        end
                        %plot(mod_T(ii,jj)-273.15,mod_H(ii,jj),'ko','MarkerSize',15);
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % 4 options to interpolate obs T-grid to model T-grid 
                    % (in addition to the "linear interpolation" which Omar used).
                    % The options are: 
                    % 1. choose nearest neighbor
                    % 2. simple mean of the neighbors
                    % 3. weighted mean of the neighbors
                    % 4. interpolation according the linear T profile of the neighbors
             
                    near_T=obs_T{f}(i_obs_grid(ii,jj,1),j_obs_grid(ii,jj,1),dd);
                    %%%%%%%%%%% 1. NEAREST NEIGHBOR %%%%%%%%
                    if strcmp(avg_type,'near_neighb')
                        obs_T_fixed{f}(ii,jj,dd)=near_T;
                    end
                    %%%%%%%%%%% 2. SIMPLE MEAN %%%%%%%%%%%%%
                    if strcmp(avg_type,'simple_mean') || (f==1 && plot_neigh(ii,jj,dd)==1)
                        mean_T=0;
                        neigh_num=4;
                        neigh_real_num=0;
                        for k=1:neigh_num
                            if ~isnan(obs_T{f}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd))
                                neigh_real_num=neigh_real_num+1;
                                mean_T=mean_T+obs_T{f}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd);
                            end    
                        end
                        if neigh_real_num==0
                            mean_T=NaN;
                        else
                            mean_T=mean_T/neigh_real_num;          %simple mean of obs "neigh_num" 2km-g.p. neighbors of the mod 2km-g.p.
                        end
                        if f==1 && plot_neigh(ii,jj,dd)==1
                            plot(mean_T,mod_H(ii,jj),'co','MarkerSize',15);     
                        end
                        if strcmp(avg_type,'simple_mean')
                            obs_T_fixed{f}(ii,jj,dd)=mean_T;
                        end
                    end
                    %%%%%%%%%%% 3. WEIGHTED MEAN %%%%%%%%%%%%%
                    if strcmp(avg_type,'weight_mean') || (f==1 && plot_neigh(ii,jj,dd)==1)
                        meanD_T=0;
                        neighD_num=min(9,tot_neigh_num);
                        weight_power=2;
                        clear d;
                        for k=1:neighD_num
                            d(k)=sqrt((lon_reduc.*(obs_lon(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k))-mod_lon(ii,jj))).^2+(obs_lat(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k))-mod_lat(ii,jj)).^2);
                        end
                        if d(1)<0.0000001
                            meanD_T=near_T;
                            display(['Perfect match for i=' num2str(ii) ' j=' num2str(jj) ' d=' num2str(d(1))]);
                        else
                            W_tot=0;
                            for k=1:neighD_num
                                if ~isnan(obs_T{f}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd))
                                    W_tot=W_tot+1/d(k)^weight_power;
                                    meanD_T=meanD_T+obs_T{f}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd)./(d(k).^weight_power);
                                end
                            end
                            if W_tot==0
                                meanD_T=NaN;
                            else
                                meanD_T=meanD_T/W_tot;        %weighted mean of "neighD_num" obs 2km-g.p. neighbors of the mod 2km-g.p.
                            end
                        end
                        if (f==1 && plot_neigh(ii,jj,dd)==1)
                            plot(meanD_T,mod_H(ii,jj),'mo','MarkerSize',15);    
                        end
                        if strcmp(avg_type,'weight_mean')
                            obs_T_fixed{f}(ii,jj,dd)=meanD_T;
                        end
                    end
                    %%%%%%%%%%% 4. LINEAR FIT FOR THE PROFILE %%%%%%%%%%%%%
                    if strcmp(avg_type,'clever_mean') || (f==1 && plot_neigh(ii,jj,dd)==1)
                        clear xx yy p;
                        for k=1:tot_neigh_plot
                            xx(k)=obs_H(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k));
                            yy(k)=obs_T{f}(i_obs_grid(ii,jj,k),j_obs_grid(ii,jj,k),dd);
                        end
                        if var(xx)<10 || var(yy)==0 || sum(isnan(yy))>0 % if somehow all the neighbours have the same height (less then 10m difference), use the nearest neighbor.
                            obs_T_fixed{f}(ii,jj,dd)=near_T;
                        else
                            p = polyfit(xx,yy,1);
                            gamma=p(1)*1000;     %the temperature profile
                            if plot_neigh(ii,jj,dd)==1
                                plot(polyval(p,xx),xx);
                            end
                            if strcmp(avg_type,'clever_mean')
                                obs_T_fixed{f}(ii,jj,dd)=polyval(p,mod_H(ii,jj));
                                if abs(obs_T_fixed{f}(ii,jj,dd)-near_T)>7 || abs(gamma)>10 %quality check threshold of 7 degrees or local profile lapse > 10 degrees per km.
                                    obs_T_fixed{f}(ii,jj,dd)=near_T;
                                    bad_arr.i=[bad_arr.i ii];
                                    bad_arr.j=[bad_arr.j jj];
                                    bad_arr.dd=[bad_arr.dd dd];
                                end
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
                end %for f=1:length(MyFolderInfo)
            end
        end
    end
end

clear obs_T_fixed_new
dinm=[31 28 31 30 31 30 31 31 30 31 30 31];
i=0;
for f=1:length(MyFolderInfo)
    f
    if ndims(obs_T_fixed{f})<3
        for d=1:dinm(f)
            i=i+1;
            obs_T_fixed_new(:,:,i)=NaN;
        end
    else
        for d=1:dinm(f)
            i=i+1;
            obs_T_fixed_new(:,:,i)=obs_T_fixed{f}(:,:,d);
        end
    end
end
clear obs_T_fixed
obs_T_fixed=obs_T_fixed_new;
clear obs_T_fixed_new

if strcmp(avg_type,'near_neighb') save([obsdir '/obs_switzeland_1km/temperature/obs_' maxminavg '_neighb'],'obs_T_fixed','-v7.3'); end;
if strcmp(avg_type,'simple_mean') save([obsdir '/obs_switzeland_1km/temperature/obs_' maxminavg '_simple'],'obs_T_fixed','-v7.3'); end;
if strcmp(avg_type,'weight_mean') save([obsdir '/obs_switzeland_1km/temperature/obs_' maxminavg '_weight'],'obs_T_fixed','-v7.3'); end;   
if strcmp(avg_type,'clever_mean') save([obsdir '/obs_switzeland_1km/temperature/obs_' maxminavg '_clever'],'obs_T_fixed','-v7.3'); end; 
  