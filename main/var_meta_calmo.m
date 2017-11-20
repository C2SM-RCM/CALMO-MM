function [varname,ofact1,ofact2,mfact1,mfact2,unit]=var_meta_calmo(v)

% Retrieve meta information of model variables that are read in
% read_calmo_sim, read_calmo_sim and read_calmo_obs
% NAME 
%   var_meta_calmo
% PURPOSE 
%   Allocate meta information for netcdf I/O and 
%
% INPUTS 
%   v: Cell array of variables names that are used
%
% OUTUTS 
%   varname:  Variable name in netcdf files
%   ofact1:   Factor for observation input
%   ofcat2:   Addition term for observation input
%   ofact1:   Factor for model input
%   ofcat2:   Addition term for model input
%   unit  :   Unit of variable
%   clat/lon: Latitude/Longitude name of variables
%
% HISTORY 
% First version: 11.10.2013
%
% AUTHOR  
%   Omar Bellprat (omar.bellprat@gmail.com)


nvar=length(v);
for l=1:nvar
  if strcmp(v(l),'t2m_avg')
    varname{l}='T2m';
    ofact1(l) = 1; %ofact1*obs+ofact2
    ofact2(l) = 0;
    mfact1(l) = 1;
    mfact2(l) = -273.15;
    unit{l} = 'degC';
    clat{l} = 'lat_1';
    clon{l} = 'lon_1';
  elseif strcmp(v(l),'t2m_min')
    varname{l}='T2m_min';
    ofact1(l) = 1; %ofact1*obs+ofact2
    ofact2(l) = 0;
    mfact1(l) = 1;
    mfact2(l) = -273.15;
    unit{l} = 'degC';
    clat{l} = 'lat_1';
    clon{l} = 'lon_1';
  elseif strcmp(v(l),'t2m_max')
    varname{l}='T2m_max';
    ofact1(l) = 1; %ofact1*obs+ofact2
    ofact2(l) = 0;
    mfact1(l) = 1;
    mfact2(l) = -273.15;
    unit{l} = 'degC';
    clat{l} = 'lat_1';
    clon{l} = 'lon_1';   
  elseif strcmp(v(l),'pr')
    varname{l}='TOT_PREC';
    ofact1(l) = 1;
    ofact2(l) = 0;
    mfact1(l) = 1;
    %mfact1(l) = 10;     %???
    %mfact1(l) = 24;     %???
    mfact2(l) = 0;
    unit{l} = 'mm/d';
    clat{l} = 'lat_1';
    clon{l} = 'lon_1';
  elseif strcmp(v(l),'sound')
    varname{l}='sounding';
    ofact1(l) = 1;
    ofact2(l) = 0;
    mfact1(l) = 1;
    mfact2(l) = 0;
    unit{l} = '-';
    clat{l} = '-';
    clon{l} = '-';  
  end;
end;
