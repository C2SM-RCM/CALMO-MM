function [dmatrix]=neelin_p(metamodel,parameters,datamatrix,pvector)

% NAME 
%   neelin_p
% PURPOSE 
%   Forecast using regression metamodel as described in Neelin et al. (2010) PNAS and Bellprat (2012).
% METHOD
%   Predict data using the metamodel for a parameter matrix
% INPUTS 
%   From the structure metamodel, parameters and datamatrix the following fields are
%   processed (mind the same naming in the input)
%   metamodel.a:
%          Vector of linear terms of the metamodel [...,N,1] additional
%          data dimensions possible (ex:a~[Regions,Variables,Time,N,1])
%   metamodel.B: 
%           Matrix of quadratic and interactions terms [...,N,N] additional
%           data dimensions possible (ex:a~[Regions,Variables,Time,N,N])
%   metamodel.c: Metamodel parameter for zero order (constant)
%   parameters.range:
%            Range of values for each paramter to normalize the scale.
%   parameters.default:
%            Default values of parameters to center the scale
%   datamatrix.reffdata:
%            Modeldata of when using default parameter values to center the datamatrix
%   pvector: Parameter values for one experiment with the
%            dimension of [N,1] N=Number parameters
% OUTPUT 
%   dmatrix: Predicted data for parameter experiment
% AUTHORS
%   Originally: Omar Bellprat (omar.bellprat@gmail.com)
%   Currently: Pavel Khain (pavelkh_il@yahoo.com)


%--------------------------------------------------------------------
% READ Input values from structures
%--------------------------------------------------------------------

a=metamodel.a;
B=metamodel.B;
N=length(parameters); % Number of model parameters
refp=parameters(1).default; % Default modelparameters
range={parameters.range}; % Parameter ranges
refd=datamatrix.refdata; % Reference data


%--------------------------------------------------------------------
% CHECK Input consistency
%--------------------------------------------------------------------
if ~exist('a','var')
  error('Please give metamodel parameter "a" as input')
end

if ~exist('B','var')
  error('Please give metamodel parameter "B" as input')
end

if ~exist('pvector','var')
  error('Please give parameter matrix as input')
end

da=size(a);
nda=length(da); %Number of dimensiones in a
db=size(B);
ndb=length(db); %Number of dimensions in B
dp=size(pvector); 
dd=da(1:end-1); %Dimension of the data

if length(pvector) > N
   error(['Parameter vector wrongly specified, only one set of' ...
	  ' parameters allowed'])
elseif db(end)~=N
  error('Parameter number in metamodel parameters a and B does agree')
elseif length(pvector)<N
  error(['Number of parameters in parameter matrix does not agree' ...
	 ' with a and B'])
end

%Normalize parameter values by the total range

for i=1:N
  varp(i)=abs(diff(range{i}));
end
pvector=(pvector-refp)./varp;


%--------------------------------------------------------------------
% ALLOCATE Output variables
%--------------------------------------------------------------------

dmatrix=NaN(dd);

%--------------------------------------------------------------------
% COMPUTE Data for each parameter experiment
%--------------------------------------------------------------------

x=reshape(pvector,[dd./dd N]);
xa=repmat(x,[dd 1]);

xi=x(:,:,:,1);          %pavel
%%%xi=x;xi(x~=0)=1;
xc=repmat(xi,[dd]);   %pavel
%%%xc=repmat(xi,[dd 1]);
x=reshape(pvector,[dd./dd 1 N]);
xb=repmat(x,[dd N 1]);
x=reshape(pvector,[dd./dd N]);

% Multiple dimension matrix operation allows for
% direct computation of dmatrix with no loop 
% (works only because B is symetric!)


if isfield(metamodel,'c') % Check if constant term present, only if
                         % neelin_c used to additionally constrain a,B
  c=metamodel.c;
  
else
  c=zeros(size(a));
end

%for each day/reg, xa is the vector of normalized parameter shifts of VAL
%(the same for any day/reg)

%for each day/reg, calculate the MM prediction, by summing up all the
%polynomial terms at that day/reg
if ndims(a)==2
  dmatrix=sum(xa.*a,nda)+sum(xa.*sum(xb.*B,ndb)',nda)+xc.*c+refd; 
else
  dmatrix=sum(xa.*a,nda)+sum(xa.*sum(xb.*B,ndb),nda)+xc.*c+refd; 
end
%%%if ndims(a)==2
%%%  dmatrix=sum(xa.*a,nda)+sum(xa.*sum(xb.*B,ndb)',nda)+sum(xc.*c,nda)+refd; 
%%%else
%%%  dmatrix=sum(xa.*a,nda)+sum(xa.*sum(xb.*B,ndb),nda)+sum(xc.*c,nda)+refd; 
%%%end
  
  
if isfield(datamatrix,'limits')
  dd=ndims(datamatrix.refdata);
  if dd>2
    for i=1:dd
      indd{i}=':';
    end
  end

  for i=1:length(datamatrix.variables)-1
    indd{datamatrix.variables{1}}=i;
    % Find values in datamatrix which are small than lower limit
    % and set to lower limit
    tmpdata=dmatrix(indd{:});
    llim=datamatrix.limits{i}(1);
    tmpdata(tmpdata < llim)= llim; 
    dmatrix(indd{:})=tmpdata;

    % Find values in datamatrix which are larger than upper limit
    % and set to upper limit
    tmpdata=dmatrix(indd{:});
    hlim=datamatrix.limits{i}(2);
    tmpdata(tmpdata > hlim)= hlim;
    dmatrix(indd{:})=tmpdata;
  end
end                                            
