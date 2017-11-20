function planes(main_data,parameters,w_user,score,new_calc,predir,param_log)

% NAME
%   planes
% PURPOSE
%   Plot performance scores for pair-wise parameters cross sections
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
%   datamatrix.reffdata:
%            Modeldata of when using default parameter values to center the datamatrix
%   w_user - array of user defined weights (for simlicity - from 0 to 1) for calibrated fields:
%            tmax tmin pr cape cin ws1 ws2 ws3 T850mb T700mb T500mb RH850mb RH700mb R500mb U850mb U700mb U500mb V850 V700mb V500mb
%   score - 'rmse' or 'cosi' for RMSE-type and COSI-type scores, respectively
%   new_calc - 0 or 1: 0 by default, when main_data is devided into cells over periods. 1 - otherwise
%   predir - path for saving output planes figures
%   param_log - Array of 0/1 numbers (having length of paramn), where ones stand for parameters which are transformed to log space
% OUTPUT
%   saved planes figures
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

%--------------------------------------------------------------------
% DEFINE Additional needed vectors
%--------------------------------------------------------------------
prd=([206 81 77]-50)./255;
pbd=([184 210 237]-100)./255;
% Hotcold colormap
hotcold=[linspace(pbd(1),1,100)' linspace(pbd(2),1,100)' linspace(pbd(3),1,100)';...
         linspace(1,prd(1),100)' linspace(1,prd(2),100)' linspace(1,prd(3),100)'];

%acc=20;  % Number of contourlevel intervals
acc=10;  % Number of contourlevel intervals

% Compute index vector for all possible pairs
pqn=allcomb(1:N,1:N);

cnt=1;
for i=1:length(pqn)
  if pqn(i,1)>=pqn(i,2)
   cind(cnt)=i;
   cnt=cnt+1;
  end
end
pqn(cind,:)=[];

% Division of panels in the plot
di=N*(N-1)/2; %Number of interactions

dv=divisor(di);
pn1=dv(floor(length(dv)/2));
pn2=di/pn1; pm=[pn1,pn2];


clear obsdata field W
if new_calc==1
    for type=1:length(main_data)
        obsdata{type}=main_data{type}.data.obsdata;
        field{type}=main_data{type}.field;
        W{type}=main_data{type}.W;
    end
else
    for mm=load_period
        for type=1:length(main_data{mm})
            obsdata{mm}{type}=main_data{mm}{type}.data.obsdata;
            field{mm}{type}=main_data{mm}{type}.field;
            W{mm}{type}=main_data{mm}{type}.W;
        end    
    end
end

%--------------------------------------------------------------------
% DEFINE Plot characteristics
%--------------------------------------------------------------------
close all;
fig=figure;
set(fig,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

myfontsize=10;

ncomb=size(pqn,1);
for i=1:ncomb
    display(['processing combination ' num2str(i) ' of ' num2str(ncomb)]);
	xgrid=linspace(range{pqn(i,1)}(1),range{pqn(i,1)}(2),acc);
	ygrid=linspace(range{pqn(i,2)}(1),range{pqn(i,2)}(2),acc);
  
	[X Y]=meshgrid(xgrid,ygrid);
	xstarp=repmat(refp,[acc^2,1]);
	xstarp(:,pqn(i,1))=X(:);
	xstarp(:,pqn(i,2))=Y(:);
	for j=1:length(xstarp)

        clear qfit
        if new_calc==1
            display(['point ' num2str(j) ' of ' num2str(length(xstarp))]);
            for type=1:length(main_data)
                qfit{type}=neelin_p(main_data{type}.metamodel,parameters,main_data{type}.data,xstarp(j,:));     %paraboloid (MM) prediction for each parameters combination
            end
        else
            for mm=load_period
                for type=1:length(main_data{mm})
                    qfit{mm}{type}=neelin_p(main_data{mm}{type}.metamodel,parameters,main_data{mm}{type}.data,xstarp(j,:));     %paraboloid (MM) prediction for each parameters combination
                end 
            end
        end

        if strcmp(score,'rmse')
            Stmp{i}(j)=rmse_score(qfit,obsdata,W,w_user,new_calc);
        elseif strcmp(score,'cosi')
            Stmp{i}(j)=cosi_score(qfit,obsdata,W,w_user,new_calc);        %can be negative !
        end
	end
end

for i=1:ncomb
  s=subplot(max(pm),min(pm),i);
  Splanes=reshape(Stmp{i},[acc acc]);
  % For CALMO only deviations shown
  Splanes=Splanes-mean(Splanes(:));
  
	display(['plotting combination ' num2str(i) ' of ' num2str(ncomb)]);
	xgrid=linspace(range{pqn(i,1)}(1),range{pqn(i,1)}(2),acc);
	ygrid=linspace(range{pqn(i,2)}(1),range{pqn(i,2)}(2),acc);
	[X Y]=meshgrid(xgrid,ygrid);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % transfer back from log space:
  nx=size(X,2);
  ny=size(Y,1);
  X_reg=X;
  Y_reg=Y;
  % checking X:
  if param_log(pqn(i,1))==1
      VecX=[];
      for ix=1:nx
          VecX(ix)=log_turlen_entrsc(X(1,ix),parameters(pqn(i,1)).name);
      end
      for iy=1:ny
          X_reg(iy,:)=VecX;
      end
  end
  if param_log(pqn(i,2))==1
      VecY=[];
      for iy=1:ny
          VecY(iy)=log_turlen_entrsc(Y(iy,1),parameters(pqn(i,2)).name);
      end
      for ix=1:nx
          Y_reg(:,ix)=VecY;
      end
  end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [ch ch]=contourf(X_reg,Y_reg,Splanes,acc);
  hold on
  % Plot design points for optimization (optional)
  xl=get(gca,'Xlim');
  yl=get(gca,'ylim');
  xlabel(parameters(pqn(i,1)).name_tex,'Fontsize',myfontsize)
  ylabel(parameters(pqn(i,2)).name_tex,'Fontsize',myfontsize)
  
  cmax=max(abs(min(Stmp{i}-mean(Stmp{i}))),abs(max(Stmp{i}-mean(Stmp{i}))));
  cmin=-cmax;
  caxis([cmin cmax])
  colormap(hotcold)
  
  
  %%%%%%%%%%%%%%%%%%%%%
    sPos = get(s,'position');
    cb=colorbar('location','eastoutside');
    set(s,'position',sPos);
  %%%%%%%%%%%%%%%%%%%%%
  
end
ax=axes('Visible','off');
set(ax,'Position',[0.116 0.09 0.8 0.2],'Fontsize',myfontsize);
set(gca,'fontname','times','fontsize',myfontsize,'fontweight','bold');
set(gcf,'color','w');
print(fig,[predir '/planes'],'-dpng');
saveas(fig,[predir '/planes']);

