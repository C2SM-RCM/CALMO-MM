function [Sopt, xstar, xopt, sc_stat,UB_next,LB_next]=lhopt(main_data,parameters,w_user,score,new_calc,lhacc,tmp_str,iteration,best_percent,UB_new,LB_new,param_log)

% NAME
%   lhopt
% PURPOSE
% Optimise model parameters using a latin hypercube sampling. See Bellprat (2012)
% METHOD
%   Create a sample of parameters using a latin hypercube design
%   and predict the model performance of the sample using the metamodel.
% INPUTS
%   main_data - big structure, which includes the sub-structures: 
%   main_data.data - datamatrix structure
%   main_data.metamodel metamodel structure
%   main_data.field - field names
%   main_data.W - weights array for different fields
%   parameters structure:
%   parameters.range:
%            Range of values for each paramter to normalize the scale.
%   parameters.default:
%            Default values of parameters to center the scale
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields:
%            tmax tmin pr cape cin ws1 ws2 ws3 T850mb T700mb T500mb RH850mb RH700mb R500mb U850mb U700mb U500mb V850 V700mb V500mb
%   score - 'rmse' or 'cosi' for RMSE-type and COSI-type scores, respectively
%   new_calc - 0 or 1: 0 by default, when main_data is devided into cells over periods. 1 - otherwise
%   lhacc - Number of experiments to sample parameter space at each iteration
%   tmp_str - path to the calibration results
%   iteration - iteration number (of convergence process to the optimal parameters combination)
%   best_percent - "winners" percent of lhacc which is used to define the parameters space for the next iteration
%   UB_new - upper limit of parameters range at given iteration
%   LB_new - lower limit of parameters range at given iteration
%   param_log - Array of 0/1 numbers (having length of paramn), where ones stand for parameters which are transformed to log space 
% OUTPUTS
%   Sopt - Scores for all experiments at given iteration
%   xstar - Latin hypercube parameter experiments at given iteration
%   xopt - Parameter setting with highest score at given iteration
%   sc_stat - score statistics at given iteration
%   UB_next - upper limit of parameters range at NEXT iteration
%   LB_next - lower limit of parameters range at NEXT iteration
% AUTHORS
%   Originally: Omar Bellprat (omar.bellprat@gmail.com)
%   Currently:  Pavel Khain (pavelkh_il@yahoo.com)




global load_period

%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------
N=length(parameters); % Number of model parameters
refp=parameters(1).default; % Default modelparameters
range={parameters.range}; % Parameter ranges

clear obsdata field W
if new_calc==1
    for type=1:length(main_data)
        obsdata{type}=main_data{type}.data.obsdata;
        refdata{type}=main_data{type}.data.refdata;
        field{type}=main_data{type}.field;
        W{type}=main_data{type}.W;
    end
else
    for mm=load_period
        for type=1:length(main_data{mm})
            obsdata{mm}{type}=main_data{mm}{type}.data.obsdata;
            refdata{mm}{type}=main_data{mm}{type}.data.refdata;
            field{mm}{type}=main_data{mm}{type}.field;
            W{mm}{type}=main_data{mm}{type}.W;
        end    
    end
end
%------------------------------------

if iteration==1
    for i=1:N
        UB(i)=range{i}(2);
        LB(i)=range{i}(1);
    end
else
    UB=UB_new;
    LB=LB_new;
end

lh=lhsdesign(lhacc,N,'criterion','maximin');
xstar=repmat(LB,[lhacc,1])+lh.*repmat((UB-LB),[lhacc,1]);

%--------------------------------------------------------------------
% PREDICT Performance of all experiments
%--------------------------------------------------------------------

%Timing variables
tic
cnt2=0;
cnt=0;
st=1000;

clear qfit_ref
if new_calc==1
    for type=1:length(main_data)
        qfit_ref{type}=neelin_p(main_data{type}.metamodel,parameters,main_data{type}.data,refp);     %paraboloid (MM) prediction for each parameters combination
    end
else
    for mm=load_period
        for type=1:length(main_data{mm})
            qfit_ref{mm}{type}=neelin_p(main_data{mm}{type}.metamodel,parameters,main_data{mm}{type}.data,refp);     %paraboloid (MM) prediction for each parameters combination
        end 
    end
end

if strcmp(score,'rmse')
    norm_sc=rmse_score(qfit_ref,obsdata,W,w_user,new_calc);
elseif strcmp(score,'cosi')
    norm_sc=cosi_score(qfit_ref,obsdata,W,w_user,new_calc);
end        

clear Sopt;
for p=1:length(xstar)
    clear qfit
    if new_calc==1
        display(['Iteration ' num2str(iteration) ' point ' num2str(p) ' of ' num2str(length(xstar))]);
        for type=1:length(main_data)
            qfit{type}=neelin_p(main_data{type}.metamodel,parameters,main_data{type}.data,xstar(p,:));     %paraboloid (MM) prediction for each parameters combination
        end
    else
        for mm=load_period
            display(['Iteration ' num2str(iteration) ' point ' num2str(p) ' of ' num2str(length(xstar))]);
            for type=1:length(main_data{mm})
                qfit{mm}{type}=neelin_p(main_data{mm}{type}.metamodel,parameters,main_data{mm}{type}.data,xstar(p,:));     %paraboloid (MM) prediction for each parameters combination
            end 
        end
    end
    
    if strcmp(score,'rmse')
        Sopt(p)=-rmse_score(qfit,obsdata,W,w_user,new_calc)./norm_sc+1;
    elseif strcmp(score,'cosi')
        Sopt(p)=cosi_score(qfit,obsdata,W,w_user,new_calc)./norm_sc-1; % may be negative!
    end
       
    %Timing code
    cnt=cnt+1;
    if cnt==st
      cnt2=cnt+cnt2;
      cnt3=length(xstar)-cnt2;
      cnt=0;
      t(cnt2/st)=toc/cnt2*st;
      display([ num2str(cnt3) ' parameter sets operations left']);
      display(['Approximately ' num2str(round(cnt3*mean(t)/st)) ' seconds left']);
    end
end


%--------------------------------------------------------------------
% FIND Best parameter set
%--------------------------------------------------------------------

max_sc=max(Sopt);     % max_sc is the highest (worst) RMSE , i.e. worst combination score
min_sc=min(Sopt);
avg_sc=mean(Sopt);
sc_stat.max=max_sc;
sc_stat.min=min_sc;
sc_stat.avg=avg_sc;

xopt=xstar(find(Sopt==max(Sopt)),:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sortedX,sortingIndices] = sort(Sopt,'descend');   % sortedX will be in descending order. Therefore, the first N elements will be the N maximum values.
params_i = sortingIndices(1:(best_percent*lhacc));
good_params=xstar(params_i,:);

for i=1:N
    UB(i)=range{i}(2);
    LB(i)=range{i}(1);
end
close all;
fig=figure;

UB_reg=UB; LB_reg=LB; refp_reg=refp;good_params_reg=good_params;xopt_reg=xopt;
for i=1:N    
    UB_next(i)=max(good_params(:,i));
    LB_next(i)=min(good_params(:,i));
    if ((UB(i)-UB_next(i))/(UB(i)-LB(i)))<0.02
        UB_next(i)=UB(i);
    end
	if ((LB_next(i)-LB(i))/(UB(i)-LB(i)))<0.02
        LB_next(i)=LB(i);
    end
    
    UB_next_reg(i)=UB_next(i); LB_next_reg(i)=LB_next(i);
    if param_log(i)==1
        %then we transfer this parameter back from log to normal
        UB_reg(i)=log_turlen_entrsc(UB(i),parameters(i).name);
        LB_reg(i)=log_turlen_entrsc(LB(i),parameters(i).name);
        refp_reg(i)=log_turlen_entrsc( refp(i),parameters(i).name);
        UB_next_reg(i)=log_turlen_entrsc(UB_next(i),parameters(i).name);
        LB_next_reg(i)=log_turlen_entrsc(LB_next(i),parameters(i).name);
        xopt_reg(i)=log_turlen_entrsc(xopt(i),parameters(i).name);
        for gp=1:size(good_params,1)
            good_params_reg(gp,i)=log_turlen_entrsc(good_params(gp,i),parameters(i).name);
        end    
    end
        
    subplot(N,1,i);
    plot(1:length(params_i),good_params_reg(:,i)); hold on;
    line([1 length(params_i)],[UB_reg(i) UB_reg(i)],'color','r'); hold on;
    line([1 length(params_i)],[LB_reg(i) LB_reg(i)],'color','r'); hold on;
    line([1 length(params_i)],[refp_reg(i) refp_reg(i)],'color','k');
    line([1 length(params_i)],[UB_next_reg(i) UB_next_reg(i)],'color','g'); hold on;
    line([1 length(params_i)],[LB_next_reg(i) LB_next_reg(i)],'color','g'); hold on;
    plot(find(sum(good_params_reg,2)==sum(xopt_reg)),xopt_reg(i),'rX','markersize',16); hold on;
    ylim([LB_reg(i)-(UB_reg(i)-LB_reg(i))/10 UB_reg(i)+(UB_reg(i)-LB_reg(i))/10]);
    ax = gca;

    ylabel(parameters(i).name_tex,'fontsize',15);
    set(gcf,'color','w');
end
filename=[tmp_str '/best_0p00' num2str(best_percent*1000) '_iter_' num2str(iteration) '_' num2str(lhacc)];
print(fig,filename,'-dpng');
saveas(fig,filename);




