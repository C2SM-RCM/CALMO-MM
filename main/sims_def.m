function [paramn,paramnt,range,default,expval,valval,simval,sims_reg,sims_inter,sims_con,valcon,param_log,date_min,date_max]=sims_def(sims_opt)

% NAME 
%   sims_def
% PURPOSE 
%   Choose parameters to tune and the simulations to use
% INPUTS 
%   sims_opt - defined by 5-digits number: sims_opt=ABCDE, where:
%              A - number of parameters to calibrate (1,2,3,4,5,6,...)
%              B - serial number of combination for given A
%              C - number of ADDITIONAL (to the minimum required) interaction parameter simulations (interaction terms)
%              D - number of "constrain" 1D simulations (additional simulations where only one parameter is changed from default)
%              E - number of parameters (among A) which are transformed to LOG space
% OUTPUTS 
%   paramn - Parameter names
%   paramnt - Parameter names (for TEX interpreter)
%   range - Parameters ranges (min and max)
%   default - Parameters defaults
%   sims_reg - Names of max-min simulations (where only 1 parameter is shifted to its max/min value)
%   sims_inter - Names of interaction simulations (where 2 parameters are shifted to their max/min values) 
%   expval - Parameters values for max-min and interaction simulations
%   simval - Name of "val" simulation (where all the parameters are shifted from their default values, in order to validate the Meta-Models)
%   valval - Parameters values for "val" simulation
%   sims_con - Name of "constrain" simulations (where each time one parameter is shifted from its default value, but not to its max/min values)
%   valcon - Parameters values for "constrain" simulations
%   param_log - Array of 0/1 numbers (having length of paramn), where ones stand for parameters which are transformed to log space 
%   date_min - The earliest allowed start date ('dd-mmm-yyyy') for chosen sims_opt 
%   date_max - The latest allowed end date ('dd-mmm-yyyy') for chosen sims_opt
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)


%--------------------------------------------------------------------------
sims_opt_original=sims_opt; % The original sims_opt
sims_opt=floor(sims_opt_original/10)*10; % Sims_opt but with zero at the end because I did not wanted to change all if commands with sims_opt , all of them have zero at end digits (the right digit)
sims_log=sims_opt_original-sims_opt; % if sims_log is equal to "0" don't calculted with log scale but if it is equal to "2" than calculated with log scal 'turl' and 'entrsc'

if sims_opt==31020
    paramn={'tkhmin','turl','entrsc'}; % No TEX interpreter
    paramnt={'tkhmin','turl','entr\_sc'}; % Parameter names
    range={[0.1 1];[100 1000];[0.05e-3 2e-3]};    % Parameter ranges (min,default,max)
    default={[0.4 150 0.3e-3]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LENTR','HENTR'};
    sims_inter={'LTKHMLTUR','LTKHMLENTR','LTURLENTR'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VENTR'};
    
    date_min='01-Jan-2013'; date_max='12-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3);
        default{1}(1) default{1}(2) 0.795e-3];
end

%--------------------------------------------------------------------------

if sims_opt==32020
    paramn={'tkhmin','turl','csoil'}; % No TEX interpreter
    paramnt={'tkhmin','turl','c\_soil'}; % Parameter names
    range={[0.1 1];[100 1000];[0 2]};    % Parameter ranges (min,default,max)
    default={[0.4 150 1]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LCSOI','HCSOI'};
    sims_inter={'LTKHMLTUR','LTKHMLCSOI','LTURLCSOI'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VCSOI'};
    
    date_min='01-Jan-2013'; date_max='6-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3);
        default{1}(1) default{1}(2) 0.5];
end

%--------------------------------------------------------------------------

if sims_opt==33130
    paramn={'tkhmin','turl','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','turl','crsmin'}; % Parameter names
    range={[0.1 1];[100 1000];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 150 150]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LCRSM','HCRSM'};
    sims_inter={'LTKHMLTUR','LCRSMHTKH','LCRSMHTUR','LCRSMLTKH'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VTUR','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='1-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3);
        default{1}(1) 316 default{1}(3);
        default{1}(1) default{1}(2) 100];
end

%--------------------------------------------------------------------------

if sims_opt==34030
    paramn={'tkhmin','entrsc','csoil'}; % No TEX interpreter
    paramnt={'tkhmin','entr\_sc','c\_soil'}; % Parameter names
    range={[0.1 1];[0.05e-3 2e-3];[0 2]};    % Parameter ranges (min,default,max)
    default={[0.4 0.3e-3 1]};
    sims_reg={'LTKHM','HTKHM','LENTR','HENTR','LCSOI','HCSOI'};
    sims_inter={'LTKHMLENTR','LTKHMLCSOI','LCSOILENTR'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VENTR','VCSOI'};
    
    date_min='01-Jan-2013'; date_max='12-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3);
        default{1}(1) 0.795e-3 default{1}(3);
        default{1}(1) default{1}(2) 0.5];
end


%--------------------------------------------------------------------------

if sims_opt==35130
    paramn={'tkhmin','entrsc','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','entr\_sc','crsmin'}; % Parameter names
    range={[0.1 1];[0.05e-3 2e-3];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 0.3e-3 150]};
    sims_reg={'LTKHM','HTKHM','LENTR','HENTR','LCRSM','HCRSM'};
    sims_inter={'LTKHMLENTR','LCRSMHTKH','LCRSMHENTR','LCRSMLTKH'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VENTR','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='3-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3);
        default{1}(1) 0.795e-3 default{1}(3);
        default{1}(1) default{1}(2) 100];
end

%--------------------------------------------------------------------------

if sims_opt==36020
    paramn={'turl','entrsc','csoil'}; % No TEX interpreter
    paramnt={'turl','entr\_sc','c\_soil'}; % Parameter names
    range={[100 1000];[0.05e-03 2e-3];[0 2]};    % Parameter ranges (min,default,max)
    default={[150 0.3e-3 1]};
    sims_reg={'LTUR','HTUR','LENTR','HENTR','LCSOI','HCSOI'};
    sims_inter={'LTURLENTR','LTURLCSOI','LCSOILENTR'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VENTR','VCSOI'};
    
    date_min='01-Jan-2013'; date_max='6-Feb-2013'; % choose your dates
    valcon=[default{1}(1) 0.795e-03 default{1}(3);
        default{1}(1) default{1}(2) 0.5];
end

%--------------------------------------------------------------------------

if sims_opt==37030
    paramn={'turl','entrsc','crsmin'}; % No TEX interpreter
    paramnt={'turl','entr\_sc','crsmin'}; % Parameter names
    range={[100 1000];[0.05e-03 2e-3];[50 200]};    % Parameter ranges (min,default,max)
    default={[150 0.3e-3 150]};
    sims_reg={'LTUR','HTUR','LENTR','HENTR','LCRSM','HCRSM'};
    sims_inter={'LTURLENTR','LCRSMHTUR','LCRSMHENTR'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTUR','VENTR','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='1-Feb-2013'; % choose your dates
    valcon=[316 default{1}(2) default{1}(3);
        default{1}(1) 0.795e-3 default{1}(3);
        default{1}(1) default{1}(2) 100];
end


%--------------------------------------------------------------------------

if sims_opt==38030
    paramn={'turl','csoil','crsmin'}; % No TEX interpreter
    paramnt={'turl','c\_soil','crsmin'}; % Parameter names
    range={[100 1000];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[150 1 150]};
    sims_reg={'LTUR','HTUR','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LTURLCSOI','LCRSMHTUR','LCRSMLCSOI'} ;
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTUR','VCSOI','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='1-Feb-2013'; % choose your dates
    valcon=[316 default{1}(2) default{1}(3);
        default{1}(1) 0.5 default{1}(3);
        default{1}(1) default{1}(2) 100];
end

%--------------------------------------------------------------------------

if sims_opt==39030
    paramn={'entrsc','csoil','crsmin'}; % No TEX interpreter
    paramnt={'entr\_sc','c\_soil','crsmin'}; % Parameter names
    range={[0.05e-03  2e-03];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.3e-03 1 150]};
    sims_reg={'LENTR','HENTR','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LCSOILENTR','LCRSMHENTR','LCRSMLCSOI'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VENTR','VCSOI','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='3-Feb-2013'; % choose your dates
    valcon=[0.795e-03 default{1}(2) default{1}(3);
        default{1}(1) 0.5 default{1}(3);
        default{1}(1) default{1}(2) 100];
end

%--------------------------------------------------------------------------

if sims_opt==310030
    paramn={'tkhmin','csoil','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','c\_soil','crsmin'}; % Parameter names
    range={[0 1];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 1 150]};
    sims_reg={'LTKHM','HTKHM','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LTKHMLCSOI','LCRSMHTKH','LCRSMLCSOI'};
    expval=[range{1}(1) default{1}(2) default{1}(3);
        range{1}(2) default{1}(2) default{1}(3);
        default{1}(1) range{2}(1) default{1}(3);
        default{1}(1) range{2}(2) default{1}(3);
        default{1}(1) default{1}(2) range{3}(1);
        default{1}(1) default{1}(2) range{3}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VCSOI','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='19-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3);
        default{1}(1) 0.5 default{1}(3);
        default{1}(1) default{1}(2) 100];
end
%--------------------------------------------------------------------------

if sims_opt==41030
    paramn={'tkhmin','turl','entrsc','csoil'}; % No TEX interpreter
    paramnt={'tkhmin','turl','entr\_sc','c\_soil'}; % Parameter names
    range={[0.1 1];[100 1000];[0.05e-3 2e-3];[0 2]};    % Parameter ranges (min,default,max)
    default={[0.4 150 0.3e-3 1]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LENTR','HENTR','LCSOI','HCSOI'};
    sims_inter={'LTKHMLTUR','LTKHMLENTR','LTKHMLCSOI','LTURLENTR','LTURLCSOI'...
        'LCSOILENTR'};        % List of 8 interactions from 1/1/2013-25/1/2013 for COSMO 1km
    expval=[range{1}(1) default{1}(2) default{1}(3) default{1}(4);
        range{1}(2) default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) range{2}(1) default{1}(3) default{1}(4);
        default{1}(1) range{2}(2) default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) range{3}(1) default{1}(4);
        default{1}(1) default{1}(2) range{3}(2) default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) range{4}(1);
        default{1}(1) default{1}(2) default{1}(3) range{4}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VENTR','VCSOI'};
    
    date_min='01-Jan-2013'; date_max='06-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) 0.795e-3 default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) 0.5];
end

%--------------------------------------------------------------------------

if sims_opt==42140
    paramn={'tkhmin','turl','entrsc','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','turl','entr\_sc','crsmin'}; % Parameter names
    range={[0.1 1];[100 1000];[0.05e-3 2e-3];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 150 0.3e-3 150]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LENTR','HENTR','LCRSM','HCRSM'};
    sims_inter={'LTKHMLTUR','LTKHMLENTR','LCRSMLTKH',...
        'LCRSMHENTR','LCRSMHTUR','LTURLENTR','LCRSMHTKH'};
    expval=[range{1}(1) default{1}(2) default{1}(3) default{1}(4);
        range{1}(2) default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) range{2}(1) default{1}(3) default{1}(4);
        default{1}(1) range{2}(2) default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) range{3}(1) default{1}(4);
        default{1}(1) default{1}(2) range{3}(2) default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) range{4}(1);
        default{1}(1) default{1}(2) default{1}(3) range{4}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VTUR','VENTR','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='01-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) 316 default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) 0.795e-3 default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) 100];
end

%--------------------------------------------------------------------------

if sims_opt==43140
    paramn={'tkhmin','turl','csoil','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','turl','c\_soil','crsmin'}; % Parameter names
    range={[0.1 1];[100 1000];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 150 1 150]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LTKHMLTUR','LTKHMLCSOI','LCRSMLTKH',...
        'LTURLCSOI','LCRSMHTUR','LCRSMLCSOI','LCRSMHTKH'};
    expval=[range{1}(1) default{1}(2) default{1}(3) default{1}(4);
        range{1}(2) default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) range{2}(1) default{1}(3) default{1}(4);
        default{1}(1) range{2}(2) default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) range{3}(1) default{1}(4);
        default{1}(1) default{1}(2) range{3}(2) default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) range{4}(1);
        default{1}(1) default{1}(2) default{1}(3) range{4}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VTUR','VCSOI','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='01-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) 316 default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) 0.5 default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) 100];
end

%--------------------------------------------------------------------------

if sims_opt==44140
    paramn={'tkhmin','entrsc','csoil','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','entr\_sc','c\_soil','crsmin'}; % Parameter names
    range={[0.1 1];[0.05e-3 2e-3];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 0.3e-3 1 150]};
    sims_reg={'LTKHM','HTKHM','LENTR','HENTR','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LTKHMLENTR','LTKHMLCSOI','LCRSMLTKH',...
        'LCSOILENTR','LCRSMHENTR','LCRSMLCSOI','LCRSMHTKH'};
    expval=[range{1}(1) default{1}(2) default{1}(3) default{1}(4);
        range{1}(2) default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) range{2}(1) default{1}(3) default{1}(4);
        default{1}(1) range{2}(2) default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) range{3}(1) default{1}(4);
        default{1}(1) default{1}(2) range{3}(2) default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) range{4}(1);
        default{1}(1) default{1}(2) default{1}(3) range{4}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VENTR','VCSOI','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='03-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) 0.795e-3 default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) 0.5 default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) 100];
end

%--------------------------------------------------------------------------

if sims_opt==45040
    paramn={'turl','entrsc','csoil','crsmin'}; % No TEX interpreter
    paramnt={'turl','entr\_sc','c\_soil','crsmin'}; % Parameter names
    range={[100 1000];[0.05e-3 2e-3];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[150 0.3e-3 1 150]};
    sims_reg={'LTUR','HTUR','LENTR','HENTR','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LTURLENTR','LTURLCSOI','LCRSMHTUR',...
        'LCSOILENTR','LCRSMHENTR','LCRSMLCSOI'};
    expval=[range{1}(1) default{1}(2) default{1}(3) default{1}(4);
        range{1}(2) default{1}(2) default{1}(3) default{1}(4);
        default{1}(1) range{2}(1) default{1}(3) default{1}(4);
        default{1}(1) range{2}(2) default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) range{3}(1) default{1}(4);
        default{1}(1) default{1}(2) range{3}(2) default{1}(4);
        default{1}(1) default{1}(2) default{1}(3) range{4}(1);
        default{1}(1) default{1}(2) default{1}(3) range{4}(2)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTUR','VENTR','VCSOI','VCRSM'};
    
    date_min='01-Jan-2013'; date_max='01-Feb-2013'; % choose your dates
    valcon=[316 default{1}(2) default{1}(3) default{1}(4);
        default{1}(1)  0.795e-3 default{1}(3) default{1}(4);
        default{1}(1) default{1}(2) 0.5 default{1}(4)
        default{1}(1) default{1}(2) default{1}(3) 100];
end

%--------------------------------------------------------------------------

if sims_opt==51150
    paramn={'tkhmin','turl','entrsc','csoil','crsmin'}; % No TEX interpreter
    paramnt={'tkhmin','turl','entr\_sc','c\_soil','crsmin'}; % Parameter names
    range={[0.1 1];[100 1000];[0.05e-3 2e-3];[0 2];[50 200]};    % Parameter ranges (min,default,max)
    default={[0.4 150 0.3e-3 1 150]};
    sims_reg={'LTKHM','HTKHM','LTUR','HTUR','LENTR','HENTR','LCSOI','HCSOI','LCRSM','HCRSM'};
    sims_inter={'LTKHMLTUR','LTKHMLENTR','LTKHMLCSOI','LCRSMLTKH',...
        'LTURLENTR','LTURLCSOI','LCRSMHTUR',...
        'LCSOILENTR','LCRSMHENTR','LCRSMLCSOI','LCRSMHTKH'};
    % List of 11 interactions from 1/1/2013-25/1/2013 for COSMO 1km (The 'sims_inter' includes 1 additional interaction)
    expval=[range{1}(1) default{1}(2) default{1}(3) default{1}(4) default{1}(5);
        range{1}(2) default{1}(2) default{1}(3) default{1}(4) default{1}(5);
        default{1}(1) range{2}(1) default{1}(3) default{1}(4) default{1}(5);
        default{1}(1) range{2}(2) default{1}(3) default{1}(4) default{1}(5);
        default{1}(1) default{1}(2) range{3}(1) default{1}(4) default{1}(5);
        default{1}(1) default{1}(2) range{3}(2) default{1}(4) default{1}(5);
        default{1}(1) default{1}(2) default{1}(3) range{4}(1) default{1}(5);
        default{1}(1) default{1}(2) default{1}(3) range{4}(2) default{1}(5);
        default{1}(1) default{1}(2) default{1}(3) default{1}(4) range{5}(1);
        default{1}(1) default{1}(2) default{1}(3) default{1}(4) range{5}(1)];
    
    for runexpval=1:length(sims_inter) % build "expval" matrix
        expvaltemp=expval_inter_combs(sims_inter{runexpval},range,default,paramn); % building the matrix "expval"
        if (length(expvaltemp)==length(paramn)) %check that returned "expvaltemp" is not empty
            expval=[expval;expvaltemp];
        else
            display('wrong sims_inter list');
            stop    % stop execution
        end
    end
    
    valval=[]; % Validation parameter sets
    simval={'MISSING'};
    sims_con={'VTKHM','VTUR','VENTR','VCSOI','VCRSM'};
    date_min='01-Jan-2013'; date_max='01-Feb-2013'; % choose your dates
    valcon=[0.7 default{1}(2) default{1}(3) default{1}(4) default{1}(5);
        default{1}(1) 316 default{1}(3) default{1}(4) default{1}(5);
        default{1}(1) default{1}(2) 0.795e-3 default{1}(4) default{1}(5);
        default{1}(1) default{1}(2) default{1}(3) 0.5 default{1}(5);
        default{1}(1) default{1}(2) default{1}(3) default{1}(4) 100];
end

%%%%%%%%%%%%%%%%%%%%%%%%% LOG SCALE CONVERTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iend=length(paramn); % THE NUMBER OF PARAMETERS
param_log=zeros(1,iend); %A retrun vector that write where are the places of  log scale and were are the places on non log scale
if(sims_log==1) % if we want at least 1 parameters with log scale than sims_log must be equal 2
     nn=0; % the number of parameters that we will find with log scale
    for i=1:iend
        if(strcmp(paramn{i},'turl'))
            nn=nn+1;
            param_log(i)=1; % assiment of log scale in param_log
            jturl=i; % the place in the paramn of the turl variable
        elseif(strcmp(paramn{i},'entrsc'))
            nn=nn+1;
            jentrsc=i; % the place in the paramn of the entrsc variable
            param_log(i)=1; % assiment of log scale in param_log
        end
    end
    if(nn==1 || nn==2) % if we have TURL and/or ENTRSC parameters than nn must be equal two
        
      nexpval=size(expval,1); % THE SIZE OF MATRIX "EXPVAL" 
      Nfirst=iend*2+1; % the first line in "EXPVAL" with the interaction simulation
       for i2=1:iend
           if(strcmp(paramn{i2},'turl'))
               % TURL LOG SCALE XMIN, XMAX and Xdefault
               xminturl=-1.386294361119891; xmaxturl=4.280132326992542; xdefturl=1.446918982936325;
              default{1}(jturl)=xdefturl;
              range{jturl}(1)=xminturl;
              range{jturl}(2)=xmaxturl;
              expval(:,jturl)=xdefturl;
              expval(jturl*2-1,jturl)=xminturl;
              expval(jturl*2,jturl)=xmaxturl;
           elseif(strcmp(paramn{i2},'entrsc'))
              % ENTRSC LOG SCALE XMIN, XMAX and Xdefault
              xminentrsc=5.347107530717469; xmaxentrsc=9.180911561285370;xdefentrsc=7.263994230454068;
              default{1}(jentrsc)=xdefentrsc;
              range{jentrsc}(1)=xminentrsc;
              range{jentrsc}(2)=xmaxentrsc;
              expval(:,jentrsc)=xdefentrsc;
              expval(jentrsc*2-1,jentrsc)=xminentrsc;
              expval(jentrsc*2,jentrsc)=xmaxentrsc;
           end
       end
       for i=Nfirst:nexpval
            if(strcmp(sims_inter{i+1-Nfirst},'LTKHMLTUR'))
                expval(i,jturl)=xminturl;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LTKHMLENTR'))
                expval(i,jentrsc)=xminentrsc;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LTURLENTR'))
                expval(i,jturl)=xminturl;
                expval(i,jentrsc)=xminentrsc;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LCSOILENTR'))
                expval(i,jentrsc)=xminentrsc;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LCRSMHENTR'))
                expval(i,jentrsc)=xmaxentrsc;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LTURLCSOI'))
                expval(i,jturl)=xminturl;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LCRSMHTUR'))
                expval(i,jturl)=xmaxturl;
            elseif(strcmp(sims_inter{i+1-Nfirst},'HTURLCSOI'))
                expval(i,jturl)=xmaxturl;
            elseif(strcmp(sims_inter{i+1-Nfirst},'LTKHMHTUR'))
                expval(i,jturl)=xmaxturl;
            end
        end
        iend2=length(sims_con); % THE NUMBER OF CONSTRAINS PARAMETERS (WE ASSUME THAT EACH PARAMETER HAS ONLY MAXIMUM ONLY ONE CONSTRAIN VALUE)
        nv=0; % the number of contrain paramters that we will find with log scale
        for i=1:iend2
            if(strcmp(sims_con{i},'VTUR'))
                nv=nv+1;
                valcon(:,jturl)=xdefturl; % change the Valcon matrix in the RAW of TURL to log default value
            elseif(strcmp(sims_con{i},'VENTR'))
                nv=nv+1;
                valcon(:,jentrsc)=xdefentrsc; % change the Valcon matrix in the RAW of ENTRSC to log default value
            end
        end
        if(nv==1 || nv==2) % if we have TURL and/or ENTRSC parameters than nn must be equal one or two
            for i=1:iend2
                if(strcmp(sims_con{i},'VTUR'))
                    xc=2.863913698933143;
                    valcon(i,jturl)=xc;
                elseif(strcmp(sims_con{i},'VENTR'))
                    xc=8.253094089655029;
                    valcon(i,jentrsc)=xc;
                end
            end
        end
    end
end
