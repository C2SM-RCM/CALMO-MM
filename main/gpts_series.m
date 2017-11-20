function datamatrix_t = gpts_series(datamatrix,vars,var)
  
% NAME 
%   gpts_series
% PURPOSE 
%   Redefine data arrays (which were not averaged over regions) to the following structure: (fields,days,regions,simulations)
% INPUTS 
%   datamatrix - Structure which includes observations and simulations data
%                for surface fields ('t2m_max','t2m_min','pr'). Dimensions:
%                [Field,Day,simulation,Lon,Lat]
%   vars - calibrated fields groups. Can be any combinations of: {'t2m_max','t2m_min','pr','sound'}
%   var - specific field group (one of vars)
% OUTPUTS 
%   datamatrix_new - Structure which includes observations and simulations
%                    data for surface fields ('t2m_max','t2m_min','pr'). Dimensions: [Field,Day,region,simulation]
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)    

datamatrix_tmp=datamatrix;

if ~strcmp(var,'sound')
    % only if var~=sound keep only one field
    for i=length(vars):-1:1
        if ~strcmp(vars{i},var)&&~strcmp(vars{i},'sound')
            %delete all the fields which are not "var" from datamatrix:
            datamatrix_tmp.obsdata(i,:,:,:,:)=[];
            datamatrix_tmp.refdata(i,:,:,:,:)=[];
            datamatrix_tmp.moddata(i,:,:,:,:)=[];
            if isfield(datamatrix, 'valdata')
                datamatrix_tmp.valdata(i,:,:,:,:)=[];
            end
            if isfield(datamatrix, 'constrain')
                datamatrix_tmp.constrain(i,:,:,:,:)=[];
            end
        end
    end
else
    datamatrix_tmp.obsdata(datamatrix_tmp.obsdata==-999.9)=NaN;
    datamatrix_tmp.refdata(datamatrix_tmp.refdata==-999.9)=NaN;
    datamatrix_tmp.moddata(datamatrix_tmp.moddata==-999.9)=NaN;
    if isfield(datamatrix, 'valdata')
        datamatrix_tmp.valdata(datamatrix_tmp.valdata==-999.9)=NaN;
    end
    if isfield(datamatrix, 'constrain')
        datamatrix_tmp.constrain(datamatrix_tmp.constrain==-999.9)=NaN;
    end
end


for f=1:size(datamatrix_tmp.moddata,1)
	for d=1:size(datamatrix_tmp.moddata,2)
        clear tmp;
        tmp=squeeze(datamatrix_tmp.obsdata(f,d,1,:,:));
        datamatrix_tmp2.obsdata(f,d,:)=tmp(:);
        clear tmp;
        tmp=squeeze(datamatrix_tmp.refdata(f,d,1,:,:));
        datamatrix_tmp2.refdata(f,d,:)=tmp(:);
        for k=1:size(datamatrix_tmp.moddata,3)
            clear tmp;
            tmp=squeeze(datamatrix_tmp.moddata(f,d,k,:,:));
            datamatrix_tmp2.moddata(f,d,:,k)=tmp(:);
        end
        if isfield(datamatrix, 'valdata')
            clear tmp;
            tmp=squeeze(datamatrix_tmp.valdata(f,d,1,:,:));
            datamatrix_tmp2.valdata(f,d,:)=tmp(:);
        end
        if isfield(datamatrix, 'constrain')
            for k=1:size(datamatrix_tmp.constrain,3)
                clear tmp;
                tmp=squeeze(datamatrix_tmp.constrain(f,d,k,:,:));
                datamatrix_tmp2.constrain(f,d,:,k)=tmp(:);
            end
        end
	end
end

%find the gpts that are always nan:
gpts_num=size(datamatrix.obsdata,4)*size(datamatrix.obsdata,5);
bad_days=[];
for f=1:size(datamatrix_tmp2.moddata,1)
	for d=1:size(datamatrix_tmp2.moddata,2)
        clear arr1; arr1=squeeze(datamatrix_tmp2.obsdata(f,d,:));
        clear arr_nan_tmp; arr_nan_tmp=isnan(arr1);
        if sum(arr_nan_tmp)~=gpts_num
            arr_nan=arr_nan_tmp;
        else
            bad_days=[bad_days d];
        end
        clear arr1; arr1=squeeze(datamatrix_tmp2.refdata(f,d,:));
        clear arr_nan_tmp; arr_nan_tmp=isnan(arr1);
        if sum(arr_nan_tmp)~=gpts_num
            arr_nan(arr_nan_tmp==1)=1;
        else
            bad_days=[bad_days d];
        end
        for m=1:size(datamatrix_tmp2.moddata,4)
            clear arr1; arr1=squeeze(datamatrix_tmp2.moddata(f,d,:,m));
            clear arr_nan_tmp; arr_nan_tmp=isnan(arr1);
            if sum(arr_nan_tmp)~=gpts_num
                arr_nan(arr_nan_tmp==1)=1;
            else
                bad_days=[bad_days d];
            end
        end
	end
end


good_gpts_num=sum(~arr_nan);

for f=1:size(datamatrix_tmp2.moddata,1)
    for d=1:size(datamatrix_tmp2.moddata,2)
        display(['this is d=' num2str(d)])
        if ~isempty(find(bad_days==d, 1))
            datamatrix_t.obsdata(f,d,1:good_gpts_num)=NaN;
            datamatrix_t.refdata(f,d,1:good_gpts_num)=NaN;
            for m=1:size(datamatrix_tmp2.moddata,4)
                datamatrix_t.moddata(f,d,1:good_gpts_num,m)=NaN;
            end
            if isfield(datamatrix, 'valdata')
                datamatrix_t.valdata(f,d,1:good_gpts_num)=NaN;
            end
            if isfield(datamatrix, 'constrain')
                for m=1:size(datamatrix_tmp2.constrain,4)
                    datamatrix_t.constrain(f,d,1:good_gpts_num,m)=NaN;
                end
            end
        else
            clear arr1; arr1=squeeze(datamatrix_tmp2.obsdata(f,d,:));
            arr1=arr1(~arr_nan);
            datamatrix_t.obsdata(f,d,:)=arr1;
            clear arr1; arr1=squeeze(datamatrix_tmp2.refdata(f,d,:));
            arr1=arr1(~arr_nan);
            datamatrix_t.refdata(f,d,:)=arr1;
            for m=1:size(datamatrix_tmp2.moddata,4)
                clear arr1; arr1=squeeze(datamatrix_tmp2.moddata(f,d,:,m));
                arr1=arr1(~arr_nan);
                datamatrix_t.moddata(f,d,:,m)=arr1;
            end
            if isfield(datamatrix, 'valdata')
                clear arr1; arr1=squeeze(datamatrix_tmp2.valdata(f,d,:));
                arr1=arr1(~arr_nan);
                datamatrix_t.valdata(f,d,:)=arr1;
            end
            if isfield(datamatrix, 'constrain')
                for m=1:size(datamatrix_tmp2.constrain,4)
                    clear arr1; arr1=squeeze(datamatrix_tmp2.constrain(f,d,:,m));
                    arr1=arr1(~arr_nan);
                    datamatrix_t.constrain(f,d,:,m)=arr1;
                end
            end
        end
    end
end
      
        