function PostProc(period)

% NAME 
%   PostProc
% PURPOSE 
%   plot analysis results and calculate the optimal parameters combination
% INPUTS 
%   time periods array (more precisely - initial dates of the periods):
%   {'dd-mmm-yyyy','dd-mmm-yyyy',...}
% OUTPUTS 
%   saved (in .mat format) analysis results
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)

format long; close all;

clear global load_period; global load_period
load_period=1:length(period); %choose which period to run
tmp_str='results';  % path to the calibration results
mkdir(tmp_str);

%------------------------------------------------------------
%-1- Choose the computer and the path:
%------------------------------------------------------------
%clear global comp maindir curdir simuldir;
global comp maindir curdir simuldir obsdir extdir;

[maindir simuldir obsdir extdir vars vars_2d avg_T vars_sound sims_opt ml score w_user lhacc iterations_num best_percent date_min_tot date_max_tot]=namelist();

%------------------------------------------------------------
%-6- Define "parameters" structure
%------------------------------------------------------------

[paramn,paramnt,range,default,expval,valval,simval,sims_reg,sims_inter,sims_con,valcon,param_log,date_min_check,date_max_check]=sims_def(sims_opt);
parameters=struct('name',paramn,'range',range','default',default,'experiments', ...
    expval,'constrain',valcon,'validation',valval,'name_tex',paramnt);


    for mm=load_period
        mm
        datamatrix{mm}=load([period{mm} '/datamatrix.mat']);
        datamatrix_s{mm}=load([period{mm} '/datamatrix_s.mat']);
        if avg_T==0
            datamatrix_tmax{mm}=load([period{mm} '/datamatrix_tmax.mat']);
            datamatrix_tmin{mm}=load([period{mm} '/datamatrix_tmin.mat']);
            datamatrix_pr{mm}=load([period{mm} '/datamatrix_pr.mat']);
        else
            datamatrix_new{mm}=load([period{mm} '/datamatrix_new.mat']);
        end
        if ~isempty(strfind([vars{:}],'sound'))
            datamatrix_so{mm}=load([period{mm} '/datamatrix_so.mat']);
        end
    end


%  NOTE !!!!! Taverage does not exist in Italy data therefore It was fixed on 14/3/2015 the Tmax and Tmin and Precipitation

    for mm=load_period
        mm
        if avg_T==0
            metamodel_tmax{mm}=load([period{mm} '/metamodel_tmax.mat']);
            metamodel_tmin{mm}=load([period{mm} '/metamodel_tmin.mat']);
            metamodel_pr{mm}=load([period{mm} '/metamodel_pr.mat']);
        else
            metamodel_new{mm}=load([period{mm} '/metamodel_new.mat']);
        end
        if ~isempty(strfind([vars{:}],'sound'))
            metamodel_so{mm}=load([period{mm} '/metamodel_so.mat']);
        end
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear global threshold; global threshold
    threshold=[0.1,1,3,7.5,10]; % for rain in case of cosi score
    
  %  score='rmse';
  %  score='cosi';
    
    if strcmp(score,'rmse')
        calc_weight=1;
    else
        calc_weight=0;
    end
  %   w_user=[1 1 1 0 0 1 0.33 0.33 0.33 0.33 0.33 0.33 0.33 0.33 0.33 0.2 0.2 0.2 0.2 0.2 0.2]; %%% ORIGINAL
    % tmax tmin pr cape cin ws1 ws2 ws3 'T850mb','T700mb','T500mb','RH850mb','RH700mb','R500mb','U850mb','U700mb','U500mb','V850','V700mb','V500mb'

    
        clear main_data;
        for mm=load_period
            mm
            i=0;
            if avg_T==0
                if ~isempty(strfind([vars{:}],'t2m_max'))
                    if calc_weight==1
                        [W_tmax{mm}]=weights_calc(parameters,datamatrix_tmax{mm},metamodel_tmax{mm},w_user(1),score,{'t2m_max'});
                    else
                        W_tmax{mm}=1;
                    end
                    i=i+1;
                    main_data{mm}{i}.data=datamatrix_tmax{mm};
                    main_data{mm}{i}.metamodel=metamodel_tmax{mm};
                    main_data{mm}{i}.field='t2m_max';
                    main_data{mm}{i}.W=W_tmax{mm};
                end
                if ~isempty(strfind([vars{:}],'t2m_min'))
                    if calc_weight==1
                        [W_tmin{mm}]=weights_calc(parameters,datamatrix_tmin{mm},metamodel_tmin{mm},w_user(2),score,{'t2m_min'});
                    else
                        W_tmin{mm}=1;
                    end
                    i=i+1;
                    main_data{mm}{i}.data=datamatrix_tmin{mm};
                    main_data{mm}{i}.metamodel=metamodel_tmin{mm};
                    main_data{mm}{i}.field='t2m_min';
                    main_data{mm}{i}.W=W_tmin{mm};
                end
                if ~isempty(strfind([vars{:}],'pr'))
                    if calc_weight==1
                        [W_pr{mm}]=weights_calc(parameters,datamatrix_pr{mm},metamodel_pr{mm},w_user(3),score,{'pr'});
                    else
                        W_pr{mm}=1;
                    end
                    i=i+1;
                    main_data{mm}{i}.data=datamatrix_pr{mm};
                    main_data{mm}{i}.metamodel=metamodel_pr{mm};
                    main_data{mm}{i}.field='pr';
                    main_data{mm}{i}.W=W_pr{mm};
                end
            else
                if calc_weight==1
                    [W_new{mm}]=weights_calc(parameters,datamatrix_new{mm},metamodel_new{mm},w_user(1:3),score,vars_2d);
                else
                    W_new{mm}=1;
                end
                i=i+1;
                main_data{mm}{i}.data=datamatrix_new{mm};
                main_data{mm}{i}.metamodel=metamodel_new{mm};
                main_data{mm}{i}.field=vars_2d;
                main_data{mm}{i}.W=W_new{mm};
            end
            if ~isempty(strfind([vars{:}],'sound'))
                if calc_weight==1
                    [W_so{mm}]=weights_calc(parameters,datamatrix_so{mm},metamodel_so{mm},w_user(4:21),score,vars_sound); %for 18 SOUNDING FILEDS
                else
                    W_so{mm}=1;
                end
                i=i+1;
                main_data{mm}{i}.data=datamatrix_so{mm};
                main_data{mm}{i}.metamodel=metamodel_so{mm};
                main_data{mm}{i}.field='so';
                main_data{mm}{i}.W=W_so{mm};
            end
        end

    save([tmp_str '/main_data.mat'],'main_data','-v7.3');
    

% (1) Check if fitted data is reproduced with 0 errors
%%%%%%ctrlpred(metamodel,parameters,datamatrix);
% (2) Estimate the error of independent simulations and plot scattorplots
%%%%%%errpred=errmeta(metamodel,parameters,datamatrix);
%%%%%%%errpred=errmeta_without_errpred(metamodel,parameters,datamatrix);
% (3) Visualize mean metamodel parameters for linear, quadratic and interaction terms
%%%%%%%metaparam(metamodel,parameters,datamatrix)
% (4) Visualize performance landscape for each parameter pair between all experiments
%planes(metamodel,parameters,datamatrix,date_lim,vars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%planes(main_data,parameters,w_user,score,new_calc,tmp_str)

    new_calc=0;
    planes(main_data,parameters,w_user,score,new_calc,tmp_str,param_log)
    
    %--------------------
    %-14- Paracleearmeters optimization
    % The analyzed metamodel can now be used to find optimal parameter
    % values and to perform a perfect model experiment.
    %--------------------
    
    % (1) Find optimal model parameters using a latin hypercube optimization
    
    %!!!!!!!!!!!!!!!!!!!!! [lhscore, lhexp, popt, sc_stat]=lhopt2(main_data,parameters,w_user,score,new_calc,lhacc);
    
    %lhacc=1000; % Number of experiments to sample parameter space
    %iterations_num=40;
    %best_percent=0.1;
    last_iteration=iterations_num;
    UB_new=0; LB_new=0;
    best=0;
    clear popt_iter UB LB lhscore_save popt_arr_max  popt_arr_min delta_popt
    for iteration=1:iterations_num
        [lhscore,lhexp,popt,sc_stat,UB_next,LB_next]=lhopt(main_data,parameters,w_user,score,new_calc,lhacc,tmp_str,iteration,best_percent,UB_new,LB_new,param_log);
        UB_new=UB_next;
        LB_new=LB_next;
        histplot(lhscore,score,best,tmp_str,iteration)
        
        popt_iter(iteration,:)=popt;
        UB(iteration,:)=UB_next;
        LB(iteration,:)=LB_next;
        lhscore_save(iteration,:)=lhscore;
        popt_arr_max(iteration)=max(lhscore);
        popt_arr_min(iteration)=min(lhscore);
        delta_popt(iteration)=(popt_arr_max(iteration)-popt_arr_min(iteration))/popt_arr_max(iteration);
        display(['iteration=' num2str(iteration) ' delta_popt=' num2str(delta_popt(iteration))]);
        if popt_arr_max(iteration)>0 && popt_arr_min(iteration)>0 && delta_popt(iteration)<10^-5
            last_iteration=iteration;
            break
        end
    end
    
    % lhscore: Modelscore for all experiments;
    % lhexp: Latin hypercube parameter experiments
    % popt: Parameter setting with highest score
    % sc_stat: sc_stat.max (sc_stat.min,sc_stat.avg) is the highest (worst) RMSE , i.e. worst combination score
        
    % find the best 10% : 
    clear delta_popt_goodenough
    iteration_goodenough=last_iteration;
    for it=1:last_iteration
        delta_popt_goodenough(it)=(popt_arr_max(last_iteration)-popt_arr_min(it))/popt_arr_max(last_iteration);
        if popt_arr_max(it)>0 && popt_arr_min(it)>0 && delta_popt_goodenough(it)<10^-1
            iteration_goodenough=it;
            break
        end
    end    
    
	% (3) Plot optimised parameter distributions
    errm=0.001; %error of metamodel
    optparam(parameters,lhscore,lhexp,popt,errm)
    
    % transfer back from log space:
    UB_reg=UB; LB_reg=LB;popt_reg=popt;
    for i=1:length(parameters)
        if param_log(i)==1
            for iteration=1:last_iteration
                %then we transfer this parameter back from log to normal
                UB_reg(iteration,i)=log_turlen_entrsc(UB(iteration,i),parameters(i).name);
                LB_reg(iteration,i)=log_turlen_entrsc(LB(iteration,i),parameters(i).name);
               % popt_arr_max_reg(iteration)=log_turlen_entrsc(popt_arr_max(i),parameters(i).name);
               % popt_arr_min_reg(iteration)=log_turlen_entrsc(popt_arr_min(i),parameters(i).name);
            end
            popt_reg(i)=log_turlen_entrsc(popt(i),parameters(i).name);
        end
    end
    % now we can plot errorbar from LB(iteration_goodenough,:) to UB(iteration_goodenough,:) arround popt
    close all;
    fig=figure();
    errorbar(1:length(popt_reg),popt_reg,popt_reg-LB_reg(iteration_goodenough,:),UB_reg(iteration_goodenough,:)-popt_reg);
    filename=[tmp_str '/best_params'];
    print(fig,filename,'-dpng');
    saveas(fig,filename);

    save([tmp_str '/UB_reg.mat'],'UB_reg','-v7.3');
    save([tmp_str '/LB_reg.mat'],'LB_reg','-v7.3');
    save([tmp_str '/lhscore_save.mat'],'lhscore_save','-v7.3');
    save([tmp_str '/popt_arr_max.mat'],'popt_arr_max','-v7.3');
    save([tmp_str '/popt_arr_min.mat'],'popt_arr_min','-v7.3');
    save([tmp_str '/popt_reg.txt'],'popt_reg','-ascii');

 
