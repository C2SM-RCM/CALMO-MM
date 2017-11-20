function main()

% NAME 
%   main
% PURPOSE 
%   main program of the CALMO parameters tuning method 
% NOTE
%   Set main definitions at namelist.m file !
% RUN
%   from Bash: "matlab -nodesktop -nosplash -r main"
%   from Matlab: F5 inside main.m
% INPUT
%   -
% OUTPUT 
%   calibration results saved in .mat format
% AUTHORS  
%   Pavel Khain (pavelkh_il@yahoo.com)
%   Itsik Carmona (carmonai@ims.gov.il)
%   Originally: Omar Bellprat (omar.bellprat@gmail.com)


%clear all; 
close all; clc; format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set calibration period and devide it to sub-periods to save MATLAB memory:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[maindir simuldir obsdir extdir vars vars_2d avg_T vars_sound sims_opt ml score w_user lhacc iterations_num best_percent date_min date_max]=namelist();

date_diff=datenum(date_max)-datenum(date_min);
period_len=10;
if date_diff<10
    display('Unreasonably small period');
    stop;
else
    period_last=floor(date_diff/period_len);
    for p=1:period_last-1
        date_min_arr{p}=datestr(datenum(date_min)+(p-1)*period_len,'dd-mmm-yyyy');
        date_max_arr{p}=datestr(datenum(date_min)+p*period_len-1,'dd-mmm-yyyy');
    end
    date_min_arr{period_last}=datestr(datenum(date_min)+period_last*period_len,'dd-mmm-yyyy');
    date_max_arr{period_last}=date_max;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 1: ReadData and fit the MetaModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for mmm=1:length(date_min_arr)
    mmm
    date_min=date_min_arr{mmm};
    date_max=date_max_arr{mmm};
    ReadData_and_MetaModel(date_min,date_max)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 2: Post-Processing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PostProc(date_min_arr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Calibration performed successfully');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%