function ets = ETS (Obs,Model,thresrain)

% NAME 
%   ETS
% PURPOSE 
%   rain forecast. Answers the question: How well did the forecast "yes" events correspond
%   to the observed "yes" events (accounting for hits that would be expected by chance)
%   Range: -1/3 to 1; 0 indicates no skill. Perfect score: 1.
% NOTE
%   be carefull with 0/0 when there is no rain both in nature and model (ex: "summer in Israel") !!!
% INPUTS 
%   parameters - structure parameters (see definitions in ReadData_and_MetaModel.m)
%   datamatrix_tmp - structure datamatrix (see definitions in ReadData_and_MetaModel.m)
%   metamodel_tmp - structure metamodel (see definitions in neelin_e.m)
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields:
%   score - 'rmse' or 'cosi' for RMSE-type and COSI-type scores, respectively
%   fields - field name (can be 't2m_max','t2m_min','pr',vars_sound)
% OUTPUT 
%   ets - rain forecast score
% AUTHOR  
%   Itsik Carmona (carmonai@ims.gov.il)


ndays=length(Model);
Total=0;
Hits=0;
Miss=0;
Falsealr=0;
Correctneg=0;
for indays=1:ndays
    if(~isnan(Model(indays)) && ~isnan(Obs(indays)))
        Total=Total+1;
        if(Model(indays)>=thresrain && Obs(indays)>=thresrain);
            Hits=Hits+1;
        elseif(Model(indays)<thresrain && Obs(indays)>=thresrain);
            Miss=Miss+1;
        elseif(Model(indays)>=thresrain && Obs(indays)<thresrain);
            Falsealr=Falsealr+1;
        else
            Correctneg=Correctneg+1;
        end
    end
end
Hitsrandom=((Hits+Falsealr)*(Hits+Miss))/Total;
ets=(Hits-Hitsrandom)/(Hits+Miss+Falsealr-Hitsrandom);
