function [maindir simuldir obsdir extdir vars vars_2d avg_T vars_sound sims_opt ml score w_user lhacc iterations_num best_percent date_min date_max]=namelist()

% NAME 
%   namelist
% PURPOSE 
%   Namelist of the calibration analysis
% INPUTS 
%   -
% OUTPUTS 
%   maindir - main directory 
%   simuldir - "maindir/simuldir": path to simulations files
%   obsdir - "maindir/obsdir": path to observations files
%   extdir - "maindir/extdir": path to "external data" files
%   vars - calibrated fields groups. Can be any combinations of: {'t2m_max','t2m_min','pr','sound'}
%   vars_2d - calibrated 2D fields. Can be any combinations of: {'t2m_max','t2m_min','pr'}
%   avg_T - region average over Precipitation only (avg_T=0), or over Precipitation, Tmax and Tmin (avg_T=1)
%   vars_sound - calibrated profiles fields: {'CAPE','CIN','TCWC','WSHEAR1','WSHEAR2','WSHEAR3','T850mb','T700mb','T500mb','RH850mb','RH700mb','RH500mb','U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'}
%   sims_opt - Choose parameters to callibrate and the simulations to use. The possible values for sims_opt and their meaning appear in sims_def.m file
%   ml - Minimum number of days (during given period) for valid soundings data. If less - current sounding fields are not analyzed
%   score - 'rmse' or 'cosi' for RMSE-type and COSI-type scores, respectively
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields:
%            tmax tmin pr cape cin ws1 ws2 ws3 T850mb T700mb T500mb RH850mb RH700mb R500mb U850mb U700mb U500mb V850 V700mb V500mb
%   lhacc -  Number of experiments to sample parameter space at each iteration
%   iterations_num - Maximum number of iterations
%   best_percent - "winners" percent of lhacc which is used to define the parameters space for the next iteration
%   date_min - beginning of calibration period
%   date_max - end of calibration period
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)


%%%%%%%%%%%%%%%%%%%%%%%%
% set calibration period
%%%%%%%%%%%%%%%%%%%%%%%%

date_min='01-Jan-2013';
date_max='01-Feb-2013';

%------------------------------------------------------------
%-1- Choose the computer and the path:
%------------------------------------------------------------
clear global comp maindir curdir simuldir;
global comp maindir curdir simuldir obsdir extdir;
%comp=1; % HERMON
%comp=2; % UNIV
%comp=3; % HOME
comp=4; % HPC
curdir='main';

% Path to program files and COSMO simulations
if comp==1
	maindir='/home/itsik/CALMO2';
	simuldir='/Research/models/CALMO2/simulations';
	obsdir='/Research/models/CALMO2/obs';
elseif comp==2
	maindir='D:\Pasha\CALMO1km';
	simuldir='D:\Pasha\CALMO1km\simulations';
	obsdir='D:\Pasha\CALMO1km\observations';
	extdir='D:\Pasha\CALMO1km\external_data';
elseif comp==3
	maindir='E:\CALMO1km';
	simuldir='E:\CALMO1km\simulations';
	obsdir='E:\CALMO1km\observations';
	extdir='E:\CALMO1km\external_data';
elseif comp==4
	maindir='/Research/models/CALMO1km';
	simuldir='/Research/models/CALMO1km/simulations';
	obsdir='/Research/models/CALMO1km/observations';
	extdir='/Research/models/CALMO1km/external_data';
end

%------------------------------------------------------------
%-2- Choose COSMO and OBS fields to be read and analyzed
%------------------------------------------------------------

%vars={'t2m_avg'};
%vars={'t2m_max'};
%vars={'t2m_min'};
%vars={'pr'};
%----------------------------------
%vars={'t2m_avg','t2m_max'};
%vars={'t2m_avg','t2m_min'};
%vars={'t2m_avg',pr'};
%vars={'t2m_max','t2m_min'};
%vars={'t2m_max','pr'};
%vars={'t2m_min','pr'};
%----------------------------------
%vars={'t2m_avg','t2m_max','t2m_min'};
%vars={'t2m_avg','t2m_max','pr'};
%vars={'t2m_avg','t2m_min','pr'};
%vars={'t2m_max','t2m_min','pr'};
%----------------------------------
%vars={'t2m_avg','t2m_max','t2m_min','pr'};
%----------------------------------
%vars={'sound'};
%----------------------------------
%vars={'t2m_avg','t2m_max','t2m_min','pr','sound'};
vars={'t2m_max','t2m_min','pr','sound'};
%vars={'t2m_max','t2m_min','pr'};
vars_2d={'t2m_max','t2m_min','pr'};
%avg_T=1;  %region average over Precipitation, Tmax and Tmin
avg_T=0;   %region average over Precipitation only
vars_sound={'CAPE','CIN','TCWC','WSHEAR1','WSHEAR2','WSHEAR3','T850mb','T700mb','T500mb','RH850mb','RH700mb','RH500mb','U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'}; % 18 fileds for SOUNDING
%------------------------------------------------------------
%-3- Choose parameters to callibrate and the simulations to use
%    according "sims_def.m". Look at this file for exact definitions
%------------------------------------------------------------
%sims_opt=31020;%- paramn={'tkhmin','turl','entrsc'}; sims_inter={'LTKHMLTUR','LTKHMLENTR','LTURLENTR'}; sims_con={Vtkhmin,Venter};
%sims_opt=31021;%- paramn={'tkhmin','log scale turl','log scale entrsc'}; sims_inter={'LTKHMLTUR','LTKHMLENTR','LTURLENTR'}; sims_con={Vtkhmin,Venter};
%sims_opt=32020;%- paramn={'tkhmin','turl','csoil'}; sims_inter={'LTKHMLTUR','LTKHMLCSOI','LCSOILENTR'}; sims_con={Vtkhmin,Vcsoil};
%sims_opt=32021;%- paramn={'tkhmin','log scale turl','csoil'}; sims_inter={'LTKHMLTUR','LTKHMLCSOI','LCSOILENTR'}; sims_con={Vtkhmin,Vcsoil};
%sims_opt=33130;%- paramn={'tkhmin','turl','crsm'}; sims_inter={'LTKHMLTUR','LCRSMHTKH','LCRSMHENTR','LCRSMLTKH'}; sims_con={Vtkhmin,Vturl,Vcrsm};
%sims_opt=33131;%- paramn={'tkhmin','log scale turl','crsm'}; sims_inter={'LTKHMLTUR','LCRSMHTKH','LCRSMHENTR','LCRSMLTKH'}; sims_con={Vtkhmin,Vturl,Vcrsm};
%sims_opt=34030;%- paramn={'tkhmin','entrc','csoi'}; sims_inter={'LTKHMLENTR','LTKHMLCSOI','LCSOILENTR'}; sims_con={Vtkhmin,Ventr,Vcsoil};
%sims_opt=34031;%- paramn={'tkhmin','log scale entrc','csoi'}; sims_inter={'LTKHMLENTR','LTKHMLCSOI','LCSOILENTR'}; sims_con={Vtkhmin,Ventr,Vcsoil};
%sims_opt=35130;%- paramn={'tkhmin','entr','crsm'}; sims_inter={'LTKHMLENTR','LCRSMHTKH','LCRSMHENTR','LCRSMLTKH'}; sims_con={Vtkhmin,Ventr,Vcrsm};
%sims_opt=35131;%- paramn={'tkhmin','log scale entr','crsm'}; sims_inter={'LTKHMLENTR','LCRSMHTKH','LCRSMHENTR','LCRSMLTKH'}; sims_con={Vtkhmin,Ventr,Vcrsm};
%sims_opt=36020;%- paramn={'turl','entr','csoil'}; sims_inter={'LTURLENTR','LTURLCSOI','LCSOILENTR'}; sims_con={Ventr,Vsoi};
%sims_opt=36021;%- paramn={'log scale of turl','log scale of entr','csoil'}; sims_inter={'LTURLENTR','LTURLCSOI','LCSOILENTR'}; sims_con={Ventr,Vsoi};
%sims_opt=37030;%- paramn={'turl','entr','crsmin'}; sims_inter={'LTURLENTR','LCRSMHTUR','LCRSMHENTR'}; sims_con={Vturl,Ventr,Vcrsm};
%sims_opt=37031;%- paramn={'log scale of turl','log scale of entr','crsmin'}; sims_inter={'LTURLENTR','LCRSMHTUR','LCRSMHENTR'}; sims_con={Vturl,Ventr,Vcrsm};
%sims_opt=38030;%- paramn={'turl','csoil','crsmin'}; sims_inter={'LTURLCSOI','LCRSMHTUR','LCRSMLSCOI'}; sims_con={Vturl,Vsoi,Vcrsm};
%sims_opt=38031;%- paramn={'log scale turl','csoil','crsmin'}; sims_inter={'LTURLCSOI','LCRSMHTUR','LCRSMLSCOI'}; sims_con={Vturl,Vsoi,Vcrsm};
%sims_opt=39030;%- paramn={'entr','csoil','crsmin'}; sims_inter={'LCSOILENTR','LCRSMHENTR','LCRSMLSCOI'}; sims_con={Ventr,Vsoi,Vcrsm};
%sims_opt=39031;%- paramn={'log scale entr','csoil','crsmin'}; sims_inter={'LCSOILENTR','LCRSMHENTR','LCRSMLSCOI'}; sims_con={Ventr,Vsoi,Vcrsm};
%sims_opt=310030;%- paramn={'tkhmin','csoil','crsmin'}; sims_inter={'LTKHMLCSOI','LCRSMHTKH','LCRSMLSCOI'}; sims_con={Vtkmin,Vsoi,Vcrsm};
%----------------------------
%sims_opt=41030;%- paramn={'tkhmin','turl','entrsc','csoil'}; sims_inter=['LTKHMLTUR','LTKHMLENTR','LTKHMLCSOI','LTURLENTR','LTURLCSOI','LCSOILENTR'}; sims_con={Vtkhmin,Ventrsc,Vsoil}
%sims_opt=41031;%- paramn={'tkhmin','LOG SCALE turl','LOG SCALE entrsc','csoil'}; sims_inter=['LTKHMLTUR','LTKHMLENTR','LTKHMLCSOI','LTURLENTR','LTURLCSOI','LCSOILENTR'}; sims_con={Vtkhmin,Ventrsc,Vsoil}
%sims_opt=42140;%-paramn={'tkhmin','turl','entrsc','crsmin'}; sims_inter={'LTKHMLTUR','LTKHMLENTR','LCRSMLTKH','LCRSMHENTR',,'LCRSMHTUR','LTURLENTR','LCRSMHTKH'};  sims_con={Vtkhmin,Vturl,Ventr,Vcrsm};
%sims_opt=42141;%-paramn={'tkhmin','LOG SCALE turl','LOG SCALE entrsc','crsmin'}; sims_inter={'LTKHMLTUR','LTKHMLENTR','LCRSMLTKH','LCRSMHENTR',,'LCRSMHTUR','LTURLENTR','LCRSMHTKH'};  sims_con={Vtkhmin,Vturl,Ventr,Vcrsm};
%sims_opt=43140;%- paramn={'tkhmin','turl','csoil,'crsm'}; sims_inter={'LTKHMLTUR','LTKHMLCSOI','LCRSMLTKH','LTURLCSOI','LCRSMHTUR','LCRSMLSCOI','LCRSMHTKH'}; sims_con={Vtkhmin,Vturl,Vsoil,Vcrsm}; 
%sims_opt=43141;%- paramn={'tkhmin','log scale turl','csoil,'crsm'}; sims_inter={'LTKHMLTUR','LTKHMLCSOI','LCRSMLTKH','LTURLCSOI','LCRSMHTUR','LCRSMLSCOI','LCRSMHTKH'}; sims_con={Vtkhmin,Vturl,Vsoil,Vcrsm};
%sims_opt=44140;%-paramn={'tkhmin','entrsc','csoil,'crsm'}`ims_inter={'LTKHMLENTR','LTKHMLCSOI','LCRSMLTKH','LCSOILENTR','LCRSMHENTR','LCRSMLSCOI','LCRSMHTKH'}; sims_con={Vtkhmin,Ventr,Vsoil,Vcrsm};
%sims_opt=44141;%-paramn={'tkhmin','log scale entrsc','csoil,'crsm'}`ims_inter={'LTKHMLENTR','LTKHMLCSOI','LCRSMLTKH','LCSOILENTR','LCRSMHENTR','LCRSMLSCOI','LCRSMHTKH'}; sims_con={Vtkhmin,Ventr,Vsoil,Vcrsm};
%sims_opt=45040;%-paramn={'turl','entrsc','csoil,'crsm'}; sims_inter={'LTURLENTR','LTURLCSOI','LCRSMHTUR','LCSOILENTR','LCRSMHENTR','LCRSMLSCOI'}; sims_con={Vturl,Ventr,Vsoil,Vcrsm}; 
%sims_opt=45041;%-paramn={'LOG SCALE turl','LOG SCALE entrsc','csoil,'crsm'}; sims_inter={'LTURLENTR','LTURLCSOI','LCRSMHTUR','LCSOILENTR','LCRSMHENTR','LCRSMLSCOI'}; sims_con={Vturl,Ventr,Vsoil,Vcrsm}; 
%----------------------------
%sims_opt=51150;%- paramn={'tkhmin','turl','entrsc','csoil,'crsm'}; sims_inter={'LTKHMLTUR','LTKHMLENTR','LTKHMLCSOI','LCRSMLTKH','LTURLENTR','LTURLCSOI','LCRSMHTUR','LCSOILENTR','LCRSMHENTR','LCRSMLSCOI','LCRSMHTKH'}; sims_con={Vtkhmin,Vturl,Ventrsc,Vcrsm,Vsoil}
sims_opt=51151;%- paramn={'tkhmin','LOG SCALE turl','LOG SCALE entrsc','csoil,'crsm'}; sims_inter={'LTKHMLTUR','LTKHMLENTR','LTKHMLCSOI','LCRSMLTKH','LTURLENTR','LTURLCSOI','LCRSMHTUR','LCSOILENTR','LCRSMHENTR','LCRSMLSCOI','LCRSMHTKH'}; sims_con={Vtkhmin,Vturl,Ventrsc,Vcrsm,Vsoil} 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% soundings quality control:
ml=6;   %if for given region and sounding field the number of valid days<ml, all the sounding fields are cancelled for that region   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chose RMSE-type and COSI-type score:
% score='rmse';
score='cosi';

% array of user defined weights (for simlicity - from 0 to 1) for calibrated fields
w_user=[1 1 1 0 0 1 0.33 0.33 0.33 0.33 0.33 0.33 0.33 0.33 0.33 0.2 0.2 0.2 0.2 0.2 0.2]; % order: tmax tmin pr cape cin ws1 ws2 ws3 T850mb T700mb T500mb RH850mb RH700mb R500mb U850mb U700mb U500mb V850 V700mb V500mb

%iterative convergence to optimal parameters combination:
lhacc=1000; % Number of experiments to sample parameter space at each iteration
iterations_num=40; % Maximum number of iterations
best_percent=0.1; % "winners" percent of lhacc which is used to define the parameters space for the next iteration
