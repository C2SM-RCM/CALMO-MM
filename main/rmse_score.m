function score=rmse_score(qfit,obsdata,W,w_user,new_calc)

% NAME 
%   rmse_score
% PURPOSE 
%   calculate RMSE-type score for Meta-Model predictions (regressions estimations)
% INPUTS 
%   qfit - metamodel predictions for given parameter combination
%   obsdata - observations data
%   W - weights for different fields, to equalize their contributions to the final score
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields 
%   new_calc - 0 or 1: 0 by default, when main_data is devided into cells over periods. 1 - otherwise
% OUTPUT
%   score - RMSE-type score
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)

global load_period
Stotal=0; % accumlated sum of S
Ntotal=0; % accumlated sum of number of existed S multiple by w_user
if new_calc==1
    i=0;
    for type=1:length(obsdata)
        std_days{type}=nanstd(obsdata{type},1,2);
        std_days{type}=repmat(std_days{type},1,size(obsdata{type},2));
        for j=1:size(obsdata{type},1)
            i=i+1;
            S(i)=0;
            if w_user(i)~=0
                tmp_mod=squeeze(qfit{type}(j,:,:));
                tmp_obs=squeeze(obsdata{type}(j,:,:));
                tmp_std=squeeze(std_days{type}(j,:,:));
                S(i)=(w_user(i)./W{type}(j)).*nanmean(((tmp_mod(:)-tmp_obs(:))./tmp_std(:)).^2);
                if(~isnan(S(i)))
                    Stotal=Stotal+S(i);
                    Ntotal=Ntotal+w_user(i);
                end
            end
        end
    end
    score=sqrt(Stotal/Ntotal);
else
    for mm=load_period
        i=0;
        for type=1:length(obsdata{mm})
            std_days{mm}{type}=nanstd(obsdata{mm}{type},1,2);
            std_days{mm}{type}=repmat(std_days{mm}{type},1,size(obsdata{mm}{type},2));
            for j=1:size(obsdata{mm}{type},1)
                i=i+1;
                S(mm,i)=0;
                if w_user(i)~=0
                    tmp_mod=squeeze(qfit{mm}{type}(j,:,:));
                    tmp_obs=squeeze(obsdata{mm}{type}(j,:,:));
                    
                    %check that std is not zero!
                    tmp_std=squeeze(std_days{mm}{type}(j,:,:));
                    S(mm,i)=(w_user(i)./W{mm}{type}(j)).*nanmean(((tmp_mod(:)-tmp_obs(:))./tmp_std(:)).^2);
                    if(~isnan(S(mm,i)))
                        Ntotal=Ntotal+w_user(i);
                        Stotal=Stotal+S(mm,i);
                    end
                end
            end
        end
    end
    score=sqrt(Stotal/Ntotal);
end

        

        