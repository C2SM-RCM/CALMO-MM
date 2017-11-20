function build_rain_obs()

% NAME 
%   build_rain_obs
% PURPOSE 
%   This function interpolates the gridded rain observations to the model grid. 
% METHOD
%   1. Read any simulation file to obtain the model grid lat/lon (ex:aggregated_LTUR_2013011000.nc)
%   2. Read gridded rain observations you need to interpolate (ex: CPCH_201301080000_01440_c2.nc) 
%   3. Interpolate the gridded observations to the model grid using nearest neighbor
% INPUTS 
%   -
% OUTPUTS 
%   saved (in .mat format) interpolated precipitation observations
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
mod_lat  = netcdf.getVar(ncid,varid_lat);
mod_lon  = netcdf.getVar(ncid,varid_lon);
mod_lon(mod_lon>180)=mod_lon(mod_lon>180)-360;
netcdf.close(ncid);
clear file ncid varid_lon varid_lat
% at this stage we have: mod_lon mod_lat

% 2. Read gridded rain observations you need to interpolate (ex: CPCH_201301080000_01440_c2.nc) 
file_rain=[obsdir '/obs_switzeland_1km/rain/CPCH_201301080000_01440_c2.nc'];
ncid_rain = netcdf.open(file_rain,'NC_NOWRITE');
varid_rain = netcdf.inqVarID(ncid_rain,'TOT_PREC');  %daily rain
varid_lon = netcdf.inqVarID(ncid_rain,'lon_1'); 
varid_lat = netcdf.inqVarID(ncid_rain,'lat_1'); 
%obs_rain = netcdf.getVar(ncid_rain,varid_rain);
m_fv=netcdf.getAtt(ncid_rain,varid_rain,'_FillValue');
%obs_rain(obs_rain==m_fv)=NaN;
obs_lon = netcdf.getVar(ncid_rain,varid_lon);
obs_lat = netcdf.getVar(ncid_rain,varid_lat);
netcdf.close(ncid_rain);
clear file_rain ncid_rain  varid_lon varid_lat
% at this stage we have:mod_lon mod_lat obs_lon obs_lat obs_rain

MyFolderInfo = dir([obsdir '/obs_switzeland_1km/rain/*.nc']);
for f=1:length(MyFolderInfo)
    %display(num2str(f));
    file_rain=[obsdir '/obs_switzeland_1km/rain/' MyFolderInfo(f).name];
    if strfind(file_rain,'missing')
        obs_rain(:,:,f)=NaN;
    else
        ncid_rain = netcdf.open(file_rain,'NC_NOWRITE');
        varid_rain = netcdf.inqVarID(ncid_rain,'TOT_PREC');  %daily rain
        obs_rain_tmp = netcdf.getVar(ncid_rain,varid_rain);
        obs_rain_tmp(obs_rain_tmp==m_fv)=NaN;
        obs_rain(:,:,f)=obs_rain_tmp;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% Can look at the readed data by:
% close all; surf(obs_lon,obs_lat,obs_rain(:,:,10),'EdgeColor','none','LineStyle','none'); colormap('parula(200)');view(0,90); colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 3. Interpolate the gridded observations to the model grid using one of the possible methods 
% Definitions:
i_size=size(mod_lat,1); j_size=size(mod_lat,2);
lon_reduc=cos(47*pi/180); % distance in longitude direction getting smaller with latitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start interpolation:
obs_rain_fixed(1:i_size,1:j_size,1:length(MyFolderInfo))=NaN;
for i=1:i_size
    i
    for j=1:j_size
        [i_tmp,j_tmp]=find((lon_reduc.*(obs_lon-mod_lon(i,j))).^2+(obs_lat-mod_lat(i,j)).^2==min(min((lon_reduc.*(obs_lon-mod_lon(i,j))).^2+(obs_lat-mod_lat(i,j)).^2)));
        if ~isnan(obs_rain(i_tmp(1),j_tmp(1)))
            obs_rain_fixed(i,j,:)=obs_rain(i_tmp(1),j_tmp(1),:);
        end
    end
end
 
save([obsdir '/obs_switzeland_1km/rain/obs_rain_neighb'],'obs_rain_fixed','-v7.3');
