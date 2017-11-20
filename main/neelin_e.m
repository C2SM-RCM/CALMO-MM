function metamodel=neelin_e(parameters, datamatrix_test,vars_2d)

% Quadratic regression metamodel as described in Neelin et al. (2010) PNAS and Bellprat (2012)
% NAME
%   neelin_e
% PURPOSE
%   Fit a mutlivariate quadratic quadratic regressions (metamodels)
% INPUTS
%   From the structure parameters and datamatrix the following fields are processed
%   parameters.experiments:
%            Parameter values for each experiment with the dimension of [N, 2*N+N*(N-1)/2]
%            The structure NEEDS to be as follows. Example for 2 parameters (p1,p2):
%            [p1_l dp2 ] ! Low parameter value for p1 default dp2
%            [p1_h dp2 ] ! High parameter value for p1 default dp2
%            [dp1  p2_l] ! Low parameter value for p2 default dp1
%            [dp1  p2_h] ! Hihg parameter value for p2 default dp1
%            [p1_l p2_h] ! Experiments with interaction (no default)
%                        ! Additional experiments used to constrain interaction terms
%   parameters.range:
%            Range of values for each parameter to normalize the scale
%   parameters.default:
%            Default values of parameters to center the scale
%   datamatrix_test.moddata:
%            Modeldata corresponding to the dimensions of parameter.experiments
%   datamatrix_test.refdata:
%            Modeldata when using default parameter values to center the datamatrix fitted
%   vars_2d - calibrated 2D fields. Can be any combinations of: {'t2m_max','t2m_min','pr'}
% OUTPUT
%   structure metamodel.
%   a: Metamodel parameter for linear terms [N,1]
%   B: Metamodel parameter for quadratic and interaction terms
%      [N,N]. Quadratic terms in the diagonal, interaction terms
%      in the off-diagonal. Matrix symmetric, B(i,j)=B(j,i).
%   c: Metamodel parameter for zero order (constant)
% AUTHORS
%   Originally: Omar Bellprat (omar.bellprat@gmail.com)
%   Currently: Pavel Khain (pavelkh_il@yahoo.com)




N=length(parameters); % Number of model parameters
refp=parameters(1).default; % Default modelparameters
range={parameters.range}; % Parameter ranges
pmatrix=parameters(1).experiments;  % Parameter combinations values 
sd=size(datamatrix_test.moddata);
nd=prod(sd(1:end-1)); % Number of datapoints   ex: 139*152=21128                               
				    

% Liearize all dimensions 
refd=datamatrix_test.refdata(:);
dvector=reshape(datamatrix_test.moddata,[prod(sd(1:end-1)),sd(end)]);    %ex: (21128,5)

%--------------------------------------------------------------------
% CHECK Input consistency
%--------------------------------------------------------------------

dm=2*N; %ex: 4
ds=2*N+N*(N-1)/2; %Number of experiments required to estimate the metamodel
di=N*(N-1)/2; %Number of all possible pairs
dp=size(pmatrix);   %ex: (5,2) - 5 combinations of 2 parameters each.
%rmsest=false; % Least-square estimation of inter-action terms
rmsest=true; % Least-square estimation of inter-action terms
intest=true; % Determination of inter-action terms 
%intest=false;  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if dp(2)~=N
  error('Dimension of pmatrix does not correspond to number of parameters')
elseif dp(1)<dm
  error('Linear set of equations is underdetermined, parameter matrix too short')
elseif dp(1)>dm && dp(1)<ds
  display(['Not enough experiments specified, interaction parameters' ...
	   'not determined'])
  intest=false;
elseif dp(1)>ds
  display(['Linear system of equations overdetermined, addtional experiments'...
           ' used to constrain interaction terms using least-squares'])
 rmsest=true;
end

%in the standard case, all the conditions above are NOT fulfilled

% Check if default value is in the center of the parameter ranges
lwb=false(1,N);upb=false(1,N);

for i=1:N
  if sum(find(pmatrix(:,i)<refp(i)))==0
    lwb(i)=true;
    display(['Default of parameter ' parameters(i).name ...
	     ' is taken at the lower bound'])
  end
  if sum(find(pmatrix(:,i)>refp(i)))==0
    upb(i)=true;
    display(['Default of parameter ' parameters(i).name ...
             ' is taken at the upper bound'])
  end
end
%in the standard case, all the conditions above are NOT fulfilled

% Normalize parameter values by the total range and center them around default - in pmatrix
% and center simulations values around default simulation values - in dvector

for i=1:N
  varp(i)=abs(diff(range{i}));
end

pmatrix=(pmatrix-repmat(refp,[dp(1),1]))./repmat(varp,[dp(1),1]);
dvector=dvector-repmat(refd,[1,sd(end)]);

%--------------------------------------------------------------------
% ALLOCATE Output variables
%--------------------------------------------------------------------

a=zeros(nd,N,1); B=zeros(nd,N,N);

%--------------------------------------------------------------------
% DEFINE Additional needed vectors
%--------------------------------------------------------------------

% Compute index vector for all possible pairs

pqn=allcomb(1:N,1:N);
cnt=1;
if N==1
    max_i=1;
else
    max_i=length(pqn);
end
for i=1:max_i
%for i=1:length(pqn)
  if pqn(i,1)>=pqn(i,2)
   cind(cnt)=i;
   cnt=cnt+1;
  end
end
pqn(cind,:)=[];

%--------------------------------------------------------------------
% DETERMINE PARAMETERS FOR MULTIVARIATE QUADRATIC MODEL
%--------------------------------------------------------------------
c=zeros(nd,1,1);
if isfield(datamatrix_test, 'constrain')
    display(['Linear system of equations overdetermined, addtional experiments'...
           ' used to constrain'])
    pmatrixc=parameters(1).constrain;
    sc=size(datamatrix_test.constrain);
    if size(datamatrix_test.constrain,4)>1
        dvectorc=reshape(datamatrix_test.constrain,[prod(sc(1:end-1)),sc(end)]);
    else
        dvectorc=reshape(datamatrix_test.constrain,[prod(sc(1:end)),1]); %I have added it, check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    dc=size(pmatrixc);
    ptmp=pmatrixc-repmat(refp,[dc(1),1]);
    % Normalize parameter values by the total range and center around default value
    pmatrixc=(pmatrixc-repmat(refp,[dc(1),1]))./repmat(varp,[dc(1),1]);
    if size(datamatrix_test.constrain,4)>1
        dvectorc=dvectorc-repmat(refd,[1,sc(end)]);
    else
        dvectorc=dvectorc-repmat(refd,[1,1]);       %I have added it, check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! not working if there are 2 constrain points!!!!!!!!!!!!!!!!!!!!!
    end
    
    l=length(pmatrix(:,1));
    for simul=1:length(pmatrixc(:,1))
        pmatrix(l+simul,:)=pmatrixc(simul,:);
        dvector(:,l+simul)=dvectorc(:,simul);
    end
end

    
if N==1
    modelterms=[1;2];
end
if N==2
    modelterms=[1 0;2 0;0 1;0 2;1 1];
end
if N==3
    modelterms=[1 0 0;2 0 0;0 1 0;0 2 0;0 0 1;0 0 2;1 1 0;1 0 1;0 1 1];
end
if N==4
    modelterms=[1 0 0 0;2 0 0 0;0 1 0 0;0 2 0 0;0 0 1 0;0 0 2 0;0 0 0 1;0 0 0 2;1 1 0 0;1 0 1 0;1 0 0 1;0 1 1 0;0 1 0 1; 0 0 1 1];
end
if N==5
    modelterms=[1 0 0 0 0 ;2 0 0 0 0;0 1 0 0 0;0 2 0 0 0;0 0 1 0 0;0 0 2 0 0;0 0 0 1 0;0 0 0 2 0;0 0 0 0 1;0 0 0 0 2;...
                1 1 0 0 0;1 0 1 0 0;1 0 0 1 0;1 0 0 0 1;0 1 1 0 0;0 1 0 1 0;0 1 0 0 1;...
                0 0 1 1 0;0 0 1 0 1;...
                0 0 0 1 1];
end
if N==6
    modelterms=[1 0 0 0 0 0;2 0 0 0 0 0;0 1 0 0 0 0;0 2 0 0 0 0;0 0 1 0 0 0;0 0 2 0 0 0;0 0 0 1 0 0;0 0 0 2 0 0;0 0 0 0 1 0;0 0 0 0 2 0;0 0 0 0 0 1;0 0 0 0 0 2;...
                1 1 0 0 0 0;1 0 1 0 0 0;1 0 0 1 0 0;1 0 0 0 1 0;1 0 0 0 0 1;0 1 1 0 0 0;0 1 0 1 0 0;0 1 0 0 1 0;0 1 0 0 0 1;...
                0 0 1 1 0 0;0 0 1 0 1 0;0 0 1 0 0 1;...
                0 0 0 1 1 0;0 0 0 1 0 1;...
                0 0 0 0 1 1];
end


ignore_param=1;    %the option to ignore the noisy parameter for given day/reg
worst_noise=0.7;   %worst allowed R^2 for single dimension

ref_index=length(pmatrix(:,1))+1;
pmatrix_new=pmatrix;
dvector_new=dvector;
pmatrix_new(ref_index,:)=0;
dvector_new(:,ref_index)=0;
for fieldayreg=1:nd %loop over fields/days/regs

    y_tmp=dvector_new(fieldayreg,:)';      %the COSMO forecast for "fieldayreg" for different parameter combinations
    %%%%%%% new: %%%%%%%%%%%%%%
    pmatrix_inloop=pmatrix_new;
    if sum(isnan(y_tmp(ds+1:end,1)))>0
        nan_index=find(isnan(y_tmp));
        nan_index(nan_index<ds+1)=[];
        y_tmp(nan_index)=[];
        pmatrix_inloop(nan_index,:)=[];
    end
    ref_index_inloop=length(y_tmp);
    pmatrix_inloop_nochange=pmatrix_inloop;
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if sum(isnan(y_tmp))>0
        % or if the current field/reg/day is nan
        c(fieldayreg,1,1)=NaN;
        a(fieldayreg,1:N,1)=NaN;
        B(fieldayreg,1:N,1:N)=NaN;      
    elseif y_tmp(:)==0    
    % don't calculate anything for zero precipitation days:
        c(fieldayreg,1,1)=0;
        a(fieldayreg,1:N,1)=0;
        B(fieldayreg,1:N,1:N)=0;     
    elseif length(find(y_tmp(:)~=0))==1    
    % don't calculate anything for zero CAPE days:
        c(fieldayreg,1,1)=0;
        a(fieldayreg,1:N,1)=0;
        B(fieldayreg,1:N,1:N)=0;     
    elseif (var(y_tmp)/mean(abs(y_tmp)))<10^(-6)
    % fit = constant default simulation value when the 1D variations are to small:
        c(fieldayreg,1,1)=0;
        a(fieldayreg,1:N,1)=0;
        B(fieldayreg,1:N,1:N)=0; 
    else                  % regular case:
        %--------------------------------------------------------
        modelterms_new=modelterms;
        y_final=y_tmp;
        
        if ignore_param==1
            noise=zeros(1,N);
            for param=1:N
                %pmatrix_single=pmatrix_new;
                pmatrix_single=pmatrix_inloop;
                dvector_single=y_tmp;
                for i=length(y_tmp):(-1):1
                    if ~(sum(abs(pmatrix_single(i,:)))==abs(pmatrix_single(i,param)))
                        pmatrix_single(i,:)=[]; %select only the one dimenstional parameter "param" changes
                        dvector_single(i)=[];   %the predictions for these parameter values 
                    end
                end
                if length(dvector_single)>3     %if we have constrain runs for single parameter "param"
                    modelterms_single=[0;1;2];
                else
                    modelterms_single=[1;2];
                end
                                
                %fieldayreg
                if sum(abs(dvector_single))~=0
                    err_occured=0;
                    try
                        parab2D = polyfitn(pmatrix_single(:,param),dvector_single,modelterms_single);
                    catch
                        err_occured=1;
                    end
                    if  err_occured==0
                        noise(param)=parab2D.R2;
                    end
                end
            end
        
            if max(noise)<worst_noise
                worst_noise=max(noise);
            end
            
            bad_params=zeros(1,N);
                        
            for param=1:N
                
                if noise(param)<worst_noise       %ignore if the parameter for that field/day/reg is very noisy
                    % --- calculate which field/day/reg is problematic here:
                    if mod(fieldayreg,sd(1)*sd(2))==0
                        reg_disp=(fieldayreg-mod(fieldayreg,sd(1)*sd(2)))/(sd(1)*sd(2));
                    else
                        reg_disp=(fieldayreg-mod(fieldayreg,sd(1)*sd(2)))/(sd(1)*sd(2))+1;
                    end
                    remain=fieldayreg-(reg_disp-1)*sd(1)*sd(2);
                    
                    if mod(remain,sd(1))==0
                        day_disp=(remain-mod(remain,sd(1)))/sd(1);
                    else
                        day_disp=(remain-mod(remain,sd(1)))/sd(1)+1;
                    end
                    remain=remain-(day_disp-1)*sd(1);
                    field_disp=remain;

                    display(['For field/day/reg= ' num2str(fieldayreg) ' i.e. for field= ' vars_2d{field_disp} ' , day= ' num2str(day_disp) ' , reg= ' num2str(reg_disp) ', the R^2 1D correlations are [' num2str(noise) ...
                    ']. So we ignore noisy parameter "' parameters(param).name '"']);
                    % ---
                    bad_params(param)=1;
                    tmp_length1=length(modelterms_new(:,1));
                    for i=tmp_length1:(-1):1     %run over n-polynomial terms (powers)
                        if modelterms_new(i,param)>0
                            modelterms_new(i,:)=[];
                        end
                    end
                    
                    % reduce the parameters number both in pmatrix and dvector:
                    tmp_length2=length(y_final);
                    for i=tmp_length2:(-1):1
                        if (sum(pmatrix_inloop(i,:))==pmatrix_inloop(i,param))&&(pmatrix_inloop(i,param)~=0)||(sum(pmatrix_inloop(i,:))~=pmatrix_inloop(i,param)&&pmatrix_inloop(i,param)~=0)
                            pmatrix_inloop(i,:)=[]; %select only the one dimenstional parameter "param" changes
                            y_final(i)=[];   %the predictions for these parameter values 
                        end
                    end
                end
            end
            %----------------------------------------------
            % here we still have zero column/columns in "modelterms_new". Kill them: 
            modelterms_new=modelterms_new(:,bad_params~=1);
            % here we still have zero column/columns in "pmatrix_final". Kill them: 
            pmatrix_inloop=pmatrix_inloop(:,bad_params~=1);        
            %----------------------------------------------
        end

        modelterms_new=unique(modelterms_new,'rows');
        N_remain=size(modelterms_new,2);
        if length(y_final)==(2*N_remain+N_remain*(N_remain-1)/2+1)         % no additional constrain simulations, then don't look for c in n-polynom
            modelterms_new(find(sum(abs(modelterms_new),2)==0),:)=[];
        end

        parab3D = polyfitn(pmatrix_inloop,y_final,modelterms_new); % Fits a general 2d polynomial regression model in several parameters dimensions

        %------------------------------- ------------------------- 
        % calculate the a,B,c coefficients:
        c(fieldayreg,1,1)=0;
        a(fieldayreg,1:N,1)=0;
        B(fieldayreg,1:N,1:N)=0;  
        sum_pow=sum(parab3D.ModelTerms,2);
        for coef=1:length(parab3D.ModelTerms(:,1)) %rows
            good_param=0;
            for n=1:N %columns
                if bad_params(n)==0
                    good_param=good_param+1;                    
                    
                    if (parab3D.ModelTerms(coef,good_param)==0)&&(sum_pow(coef)==0)
                        c(fieldayreg,1,1)=parab3D.Coefficients(coef);
                    end           
                    if (parab3D.ModelTerms(coef,good_param)==1)&&(sum_pow(coef)==1)
                        a(fieldayreg,n,1)=parab3D.Coefficients(coef);
                    end
                    if (parab3D.ModelTerms(coef,good_param)==2)&&(sum_pow(coef)==2)
                        B(fieldayreg,n,n)=parab3D.Coefficients(coef);
                    end
                    if (parab3D.ModelTerms(coef,good_param)==1)&&(sum_pow(coef)==2)
                        tmp_arr=parab3D.ModelTerms(coef,:);
                        % push zeros back at their places:
                        tmp_arr=zeros(1,N);
                        real_place=0;
                        for place=1:N
                            if bad_params(place)==0
                                real_place=real_place+1;
                                tmp_arr(place)=parab3D.ModelTerms(coef,real_place);                                            
                            end
                        end
                        ind=find(tmp_arr==1);
                        B(fieldayreg,ind(1),ind(2))=0.5*parab3D.Coefficients(coef);
                        B(fieldayreg,ind(2),ind(1))=0.5*parab3D.Coefficients(coef);
                    end
                end
            end
        end
    end
end

% reshape a and B to original data structure
a=reshape(a,[sd(1:end-1),N]);
B=reshape(B,[sd(1:end-1),N,N]);
c=reshape(c,[sd(1:end-1)]);

metamodel.a=a; metamodel.B=B; metamodel.c=c;

