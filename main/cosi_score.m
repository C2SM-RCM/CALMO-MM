function score=cosi_score(qfit,obsdata,W,w_user,new_calc)

% NAME 
%   cosi_score
% PURPOSE 
%   calculate COSI score for Meta-Model predictions (regressions
%   estimations). Defined on the basis of the COSI score by Ulrich Damrath (DWD)
% INPUTS 
%   qfit - metamodel predictions for given parameter combination
%   obsdata - observations data
%   W - weights for different fields, to equalize their contributions to the final score
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields 
%   new_calc - 0 or 1: 0 by default, when main_data is devided into cells over periods. 1 - otherwise
% OUTPUT
%   score - COSI score
% AUTHOR  
%   Itsik Carmona (carmonai@ims.gov.il)


global load_period threshold

for mm=load_period % I caclulate it for "load_period" scalar
    nii=0; % number of skindex that is not infinte or NaN
    i=0;
    w_user_sum=0;
    for type=1:length(obsdata{mm})
        
        for j=1:size(obsdata{mm}{type},1)
            i=i+1;
            skindex(i)=0; % every change of field there is a new skindex which is zero at the beginning before the RMSE or COSI calculations
            if w_user(i)~=0
                tmp_mod=squeeze(qfit{mm}{type}(j,:,:));
                tmp_obs=squeeze(obsdata{mm}{type}(j,:,:));
                
                nr=size(tmp_obs,2); % the number of regions
            
                if(type==1 && j==3) % The rain filed is when type==1 & j==3
                    skprec=zeros(1,length(threshold));
                    valid_threshs=0;
                    for i_thresh=1:length(threshold)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        valid_regs=0;
                        for nreg=1:nr                               % loop over regions
                            ets=ETS(tmp_mod(:,nreg)',tmp_obs(:,nreg)',threshold(i_thresh));  % (-1/3)<ETS<1 % the area dimension  is the left and the number of days is in the right dimension .
                            if isfinite(ets)
                                valid_regs=valid_regs+1;
                                skprec(i_thresh)=skprec(i_thresh)+ets;      % after this loop, skprec is the sum(ETS) over regions (representig all the time period), for each threshold separately
                            end
                        end
                        if(valid_regs>0)                                    % valid_regs=0 if ETS=infinite for all the regions for threshold i_thresh
                            skprec(i_thresh)=skprec(i_thresh)/valid_regs;   % region-averaged ETS for given threshold i_thresh
                            valid_threshs=valid_threshs+1;                              % the number of good thresholds
                        end
                    end
                    if(valid_threshs>0)                                       % valid_threshs=0 if ETS=infinite for all the regions for all thresholds
                        %skindex(i)=(w_user(i)./W{mm}{type}(j)).*sum(skprec)/valid_threshs;                 % skindex(ii) is the score (averaged over thresholds and regions) for field ii
                        skindex(i)=(w_user(i)./1).*sum(skprec)/valid_threshs;                 % skindex(ii) is the score (averaged over thresholds and regions) for field ii
                        nii=nii+1;                                  % number of valid sk(field) scores
                        w_user_sum=w_user_sum+w_user(i);
                    end

                else % all other fields except rain
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    rmse1=nanmean((tmp_mod(:)-tmp_obs(:)).^2);
                    clear tmp_obs_prev
                    for nreg=1:nr
                        tmp_obs_prev(:,nreg)=[NaN,tmp_obs(1:end-1,nreg)'];
                    end
                    rmse2=nanmean((tmp_obs(:)-tmp_obs_prev(:)).^2);    
                    skindex(i)=(w_user(i)/1)*(1-rmse1/rmse2);
                    if isnan(skindex(i))
                        display('nan score !!!')
                    else
                        nii=nii+1; % number of valid sk(field) scores
                        w_user_sum=w_user_sum+w_user(i);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
    if(nii>0)
        cosi_tmp(mm)=nansum(skindex)./(nii.*w_user_sum);
    else
        cosi_tmp(mm)=NaN;
    end
end

cosi=nanmean(cosi_tmp);
score=cosi;

