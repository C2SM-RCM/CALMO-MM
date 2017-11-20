function histplot(lhscore,score,best,predir,iteration)  

% NAME
%   histplot
% PURPOSE
%   plot histogram of SCORES for meta-models predictions
% INPUTS
%   lhscore - Scores for all experiments at given iteration
%   score - 'rmse' or 'cosi' for RMSE-type and COSI-type scores, respectively
%   best - 0 or 1: o by default, 1 if a special simulation exists and verified
%   predir - path for saving output planes figures
%   iteration - iteration number (of convergence process to the optimal parameters combination)
% OUTPUTS
%   saved scores histogram for specific iteration
% AUTHORS
%   Originally: Omar Bellprat (omar.bellprat@gmail.com)
%   Currently:  Pavel Khain (pavelkh_il@yahoo.com)




%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

if best==1
%    vald=datamatrix.valdata; % Best data (defined as "val")
end

% New colors
pr=([206 81 77])./255; 
pb=([184 210 237])./255;
  

%--------------------------------------------------------------------
% PLOT Metamodel range
%-------------------------------------------------------------------- 

fig=figure;
[hi hx]=hist(lhscore,200);
hi=hi/sum(hi);
fhy=[hi,zeros(1,length(hi))]
fhx=[hx,flipdim(hx,2)]
hhi=fill(fhx,fhy,pb, 'EdgeColor',pb,'Linewidth',2);
hold on
lht=max(hi);
href=plot(ones(1,100)*PSref,linspace(0,lht,100),'Linewidth',2,'color','k')
text(PSref,lht+0.0005,'REF','Rotation',90,'Fontsize',12);


set(gca,'Fontsize',18,'YTick',[],'Layer','top','Box','on','TickDir','in', 'Linewidth',1)
ylabel('Relavtive densitiy','Fontsize',18)
xlabel('Score','Fontsize',18)
title('Objective calibration','Fontsize',18)
ylim([0 lht+.004])
if strcmp(score,'ps')
    xlim([0 1])
elseif strcmp(score,'rmse') || strcmp(score,'cosi')
    % normalization by refd score option:
    xlim([min(hx)-0.05 max(hx)+0.05])
end
hl=legend([href,hhi],'Rerefence','Metamodel Range',2);    
set(hl,'Box','off')
set(gcf,'Paperposition',[1 1 20 5])
set(gcf, 'Renderer', 'painters')

set(gcf,'color','w');

print(fig,[predir '/histr_iter_' num2str(iteration)],'-dpng');
saveas(fig,[predir '/histr_iter_' num2str(iteration)]);

