function datatemp=del_bad_sounding(datamatrix_so,ml,sims_reg,sims_inter,sims_con,simval,vars_sound)
 
% NAME 
%   del_bad_sounding
% PURPOSE 
%   delete unrealistic sounding/profiles data
% INPUTS 
%   datamatrix_so - Data matrix for profiles. Dimensions: [Field,Day,region,simulation]
%   ml - Minimum number of days (during given period) for valid soundings data. If less - current sounding fields are not analyzed.
%   sims_reg - Names of max-min simulations (where only 1 parameter is shifted to its max/min value)
%   sims_inter - Names of interaction simulations (where 2 parameters are shifted to their max/min values) 
%   sims_con - Name of "constrain" simulations (where each time one parameter is shifted from its default value, but not to its max/min values)
%   simval - Name of "val" simulation (where all the parameters are shifted from their default values, in order to validate the Meta-Models)
%   vars_sound - calibrated profiles fields: {'CAPE','CIN','TCWC','WSHEAR1','WSHEAR2','WSHEAR3','T850mb','T700mb','T500mb','RH850mb','RH700mb','RH500mb','U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'}
% OUTPUTS 
%   datatemp - Data matrix for profiles. Dimensions: [Field,Day,region,simulation]
% AUTHOR  
%   Itsik Carmona (carmonai@ims.gov.il)


%%%%%%%%%%%%% problematic sounding field values in some days, so we turn those values to NaN  
bad_so_dates={'09-Jan-2013','01-Feb-2013','02-Feb-2013','03-Feb-2013','11-Nov-2013'};
bad_so_fields{1}={'T850mb','T700mb','T500mb'};
bad_so_fields{2}={'U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'};
bad_so_fields{3}={'U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'};
bad_so_fields{4}={'U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'};
bad_so_fields{5}={'U850mb','U700mb','U500mb','V850mb','V700mb','V500mb'};

for bd=1:length(bad_so_dates)
    bad_list.day=0;
    if datenum(bad_so_dates{bd})>=datenum(date_min) && datenum(bad_so_dates{bd})<=datenum(date_max)
        bad_list.day=datenum(bad_so_dates{bd})-datenum(date_min)+1;
        bad_list.vars=bad_so_fields{bd};
        if (bad_list.day~=0)
            clear datamatrix_so_new
            datamatrix_so_new=corr_so(datamatrix_so,bad_list,vars_sound,sims_reg,sims_inter,sims_con,simval);
            datamatrix_so=datamatrix_so_new;
            clear datamatrix_so_new
        end
    end
end

%%% Assign NaN in Sonde fields that have less than ml observations per month or period. The Assigment is for all days in the specific month and reigons for that specific month and fields. 

reg=0;
for i=1:size(datamatrix_so.obsdata,3)   % loop over "regions"
    good_fields=0;
    for jj=1:size(datamatrix_so.obsdata,1)  % loop over sounding fields
         if(length(find(~isnan(datamatrix_so.obsdata(jj,:,i))))>=ml)
                  good_fields=good_fields+1;    %count number of "good" fields
         end
    end
    
    good_fields
    size(datamatrix_so.obsdata,1)
    if (good_fields==size(datamatrix_so.obsdata,1))
        reg=reg+1;
        datatemp.obsdata(:,:,reg)=datamatrix_so.obsdata(:,:,i);
        datatemp.refdata(:,:,reg)=datamatrix_so.refdata(:,:,i);
        if ~strcmp(sims_inter,'MISSING') | ~strcmp(sims_reg,'MISSING')
            datatemp.moddata(:,:,reg,:)=datamatrix_so.moddata(:,:,i,:);
        end
        if ~strcmp(sims_con,'MISSING')
            datatemp.constrain(:,:,reg,:)=datamatrix_so.constrain(:,:,i,:);
        end
        if ~strcmp(simval,'MISSING')
            datatemp.val(:,:,reg)=datamatrix_so.val(:,:,i);
        end
    end
    
end
        
        
