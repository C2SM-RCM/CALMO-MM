function [W_fin]=weights_calc(parameters,datamatrix_tmp,metamodel_tmp,w_user,score,fields)

% NAME 
%   weights_calc
% PURPOSE 
%   calculate weights for different fields, to equalize their contributions to the final score (assuming user defined weights are uniform).
% INPUTS 
%   parameters - structure parameters (see definitions in ReadData_and_MetaModel.m)
%   datamatrix_tmp - structure datamatrix (see definitions in ReadData_and_MetaModel.m)
%   metamodel_tmp - structure metamodel (see definitions in neelin_e.m)
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields
%   score - 'rmse' or 'cosi' for RMSE-type and COSI-type scores, respectively
%   fields - field name (can be 't2m_max','t2m_min','pr',vars_sound)
% OUTPUTS 
%   W_fin - weights array for different fields
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)


global threshold;

lhacc=1000;
N=length(parameters); % Number of model parameters
range={parameters.range}; % Parameter ranges
obsdata=datamatrix_tmp.obsdata;

std_days = nanstd(obsdata,1,2);                           % 1 - is just some flag, 2 is the dimension for which we calculate std (days)
std_days=repmat(std_days,1,size(obsdata,2));	% check that std is not zero!

for i=1:N
  UB(i)=range{i}(2);
  LB(i)=range{i}(1);
end

lh=lhsdesign(lhacc,N,'criterion','maximin');
xstar=repmat(LB,[lhacc,1])+lh.*repmat((UB-LB),[lhacc,1]);


for p=1:length(xstar)
    qfit=neelin_p(metamodel_tmp,parameters,datamatrix_tmp,xstar(p,:));
    for f=1:size(obsdata,1)
        if w_user(f)==0
            W(f,p)=NaN;
        else
            clear tmp_mod; tmp_mod=squeeze(qfit(f,:,:));
            clear tmp_obs; tmp_obs=squeeze(obsdata(f,:,:));
            nr=size(tmp_obs,2); % the number of regions
            if strcmp(score,'rmse')
                clear tmp_std; tmp_std=squeeze(std_days(f,:,:));
                W(f,p)=nanmean(((tmp_mod(:)-tmp_obs(:))./tmp_std(:)).^2);
            elseif strcmp(score,'cosi')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(fields{f},'pr') % The rain field is when type==1 & j==3
                    skprec=zeros(1,length(threshold));
                    valid_threshs=0;
                    for i_thresh=1:length(threshold)
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
                        W(f,p)=sum(skprec)/valid_threshs;                 % skindex(ii) is the score (averaged over thresholds and regions) for field ii
                    end

                else % all other fields except rain
                    rmse1=nanmean((tmp_mod(:)-tmp_obs(:)).^2);
                    clear tmp_obs_prev
                    for nreg=1:nr
                        tmp_obs_prev(:,nreg)=[NaN,tmp_obs(1:end-1,nreg)'];
                    end

                    rmse2=nanmean((tmp_obs(:)-tmp_obs_prev(:)).^2);    

                    W(f,p)=1-rmse1/rmse2;
                    if isnan(W(f,p))
                        error('nan weight')
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end
W_fin=mean(W,2);
        

        