function xnolog=log_turlen_entrsc(xlog,paramname) 

% NAME 
%   log_turlen_entrsc
% PURPOSE 
%   convert the optimal parameters (tur_len and entr_sc) values from log-space back to the regular space
% INPUTS 
%   xlog - input vector of paramaeters values in log space
%   paramname - parameter names
% OUTPUT
%   xnolog - output vector of paramaeters values (tur_len and entr_sc) in regular space
% AUTHOR  
%   Itsik Carmona (carmonai@ims.gov.il)

if strcmp(paramname,'turl')
    xmaxmin=[100,1000];    
	xdefault=150;
    xlogminmax=[-1.386294361119891,4.280132326992542];
    xlogdef=1.446918982936325;
elseif strcmp(paramname,'entrsc')
    xmaxmin=[0.05e-3,2e-3];
	xdefault=0.3e-3;
    xlogminmax=[5.347107530717469,9.180911561285370];
    xlogdef=7.263994230454068;
end

b=exp(xlogminmax(1));
a=(exp(xlogdef)-b)*(xmaxmin(2)-xmaxmin(1))/(xdefault-xmaxmin(1));
xnolog=((exp(xlog)-b)*(xmaxmin(2)-xmaxmin(1))/a)+xmaxmin(1);
testlog=log(a*(xnolog-xmaxmin(1))/(xmaxmin(2)-xmaxmin(1))+b);   %check 
        