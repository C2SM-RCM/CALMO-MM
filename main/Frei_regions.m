function datamatrix_new = Frei_regions(datamatrix,lat,lon,vars,avg_T,unify_regions)

% NAME 
%   Frei_regions
% PURPOSE 
%   Average part of the surface fields over predefined big regions
% METHOD
%   use image file (regions_italy_swiss_for_matlab.bmp) where each region has its color
% INPUTS 
%   datamatrix - Structure which includes observations and simulations data
%                for surface fields ('t2m_max','t2m_min','pr'). Dimensions: [Field,Day,simulation,Lon,Lat]
%   lat - latitudes of model domain 
%   lon - longitudes of model domain 
%   vars - calibrated fields groups. Can be any combinations of: {'t2m_max','t2m_min','pr','sound'}
%   avg_T - region average over Precipitation only (avg_T=0), or over Precipitation, Tmax and Tmin (avg_T=1)
%   unify_regions - array that defines which regions (out of 1-7) to unify (option to unify several regions into one bigger)
% OUTPUTS 
%   datamatrix_new - Structure which includes observations and simulations
%                    data for surface fields ('t2m_max','t2m_min','pr'). Dimensions: [Field,Day,region,simulation]
% AUTHOR  
%   Pavel Khain (pavelkh_il@yahoo.com)


global extdir
clear newregs_num  newregs_size reg_test

% 1. Reduction of datamatrix_tmp.moddata to inlude only the fields which are averaged over regions

if avg_T==1
    datamatrix_tmp=datamatrix;
    
else
    if isempty(strfind(vars,'pr'))
        return
    else
        index=find(strcmp(vars,'pr'));
        vars={'pr'};
        datamatrix_tmp=datamatrix;
        datamatrix_tmp.obsdata(1,:,:,:,:)= datamatrix.obsdata(index,:,:,:,:);
        datamatrix_tmp.obsdata(2:end,:,:,:,:)=[];
        datamatrix_tmp.refdata(1,:,:,:,:)= datamatrix.refdata(index,:,:,:,:);
        datamatrix_tmp.refdata(2:end,:,:,:,:)=[];
        datamatrix_tmp.moddata(1,:,:,:,:)= datamatrix.moddata(index,:,:,:,:);
        datamatrix_tmp.moddata(2:end,:,:,:,:)=[];
        if isfield(datamatrix, 'valdata')
            datamatrix_tmp.valdata(1,:,:,:,:)= datamatrix.valdata(index,:,:,:,:);
            datamatrix_tmp.valdata(2:end,:,:,:,:)=[];
        end
        if isfield(datamatrix, 'constrain')
            datamatrix_tmp.constrain(1,:,:,:,:)= datamatrix.constrain(index,:,:,:,:);
            datamatrix_tmp.constrain(2:end,:,:,:,:)=[];
        end
    end
end

area_minimum=50;
prec_tresh=400;
newregs_num=7;

reg_test(1:size(datamatrix_tmp.refdata,2),1:size(lat,1),1:size(lat,2))=NaN;
%%%reg_test(1:size(lat,1),1:size(lat,2))=NaN;

img = imread([extdir '/regions_italy_swiss_for_matlab.bmp']);


% 2. For each i,j grid-point define to which region it belongs, and put it into reg_test(day,i,j), in case the field value at that point is not NaN.
size(datamatrix_tmp.refdata)
for day=1:size(datamatrix_tmp.refdata,2)
    day
    for i=1:size(lat,1)
        for j=1:size(lat,2)
            clear check
            check=1;
            for field=1:size(datamatrix_tmp.refdata,1)
                if isnan(datamatrix_tmp.obsdata(field,day,1,i,j)) || isnan(datamatrix_tmp.refdata(field,day,1,i,j))
                    check=NaN;
                end    
                if ~isnan(check)
                    for mod=1:size(datamatrix_tmp.moddata,3)
                        if isnan(datamatrix_tmp.moddata(field,day,mod,i,j))
                            check=NaN;
                        end
                    end
                end
                if ~isnan(check)
                    if strcmp(vars{field},'pr') && ( datamatrix_tmp.obsdata(field,day,1,i,j)>prec_tresh || datamatrix_tmp.obsdata(field,day,1,i,j)<0 ) 
                        check=NaN;
                    end
                end
               
            end
            if isnan(check)
                reg_test(day,i,j)=NaN;
            else
                % reg_test(day,i,j)=define_Frei_regions(lat(i,j),lon(i,j));
                areatmp=regions_bmp(lat(i,j),lon(i,j),img,unify_regions); % temporary area
                if(areatmp<=newregs_num) 
                    reg_test(day,i,j)=areatmp; % only for regions that are smaller or equal 7 otherwise it is undeifend area (NaN)
                end
            end
        end
    end
end 
        
% 3. Define datamatrix_new with the new dimensions (regions instead of grid-points)

num_fields=size(datamatrix_tmp.refdata,1);
datamatrix_new.refdata(1:num_fields,1:size(datamatrix_tmp.refdata,2),1:newregs_num)=0;
datamatrix_new.moddata(1:num_fields,1:size(datamatrix_tmp.refdata,2),1:newregs_num,1:size(datamatrix_tmp.moddata,3))=0;
if isfield(datamatrix_tmp, 'valdata')
    datamatrix_new.valdata(1:num_fields,1:size(datamatrix_tmp.refdata,2),1:newregs_num)=0;
end
if isfield(datamatrix_tmp, 'constrain')
    datamatrix_new.constrain(1:num_fields,1:size(datamatrix_tmp.refdata,2),1:newregs_num,1:size(datamatrix_tmp.constrain,3))=0;
end
datamatrix_new.obsdata(1:num_fields,1:size(datamatrix_tmp.obsdata,2),1:newregs_num)=0;

% 4. For averaging, aggregate the field values for each region

for num_f=1:num_fields
    for day=1:size(datamatrix_tmp.refdata,2)
        display(['averaging field ' num2str(num_f) ' day' num2str(day)]);
        tot_area(1:newregs_num)=0;
        for i=1:size(datamatrix_tmp.refdata,4)
            for j=1:size(datamatrix_tmp.refdata,5)
                if ~isnan(reg_test(day,i,j))
                    cur_reg=reg_test(day,i,j);
                    
                    datamatrix_new.refdata(num_f,day,cur_reg)=datamatrix_new.refdata(num_f,day,cur_reg)+datamatrix_tmp.refdata(num_f,day,1,i,j);
                    if isfield(datamatrix_tmp, 'valdata')
                        datamatrix_new.valdata(num_f,day,cur_reg)=datamatrix_new.valdata(num_f,day,cur_reg)+datamatrix_tmp.valdata(num_f,day,1,i,j);
                    end
                    if isfield(datamatrix_tmp, 'constrain')
                        for mod=1:size(datamatrix_tmp.constrain,3) 
                            datamatrix_new.constrain(num_f,day,cur_reg,mod)=datamatrix_new.constrain(num_f,day,cur_reg,mod)+datamatrix_tmp.constrain(num_f,day,mod,i,j);
                        end
                    end
                    datamatrix_new.obsdata(num_f,day,cur_reg)=datamatrix_new.obsdata(num_f,day,cur_reg)+datamatrix_tmp.obsdata(num_f,day,1,i,j);
                    for mod=1:size(datamatrix_tmp.moddata,3) 
                        datamatrix_new.moddata(num_f,day,cur_reg,mod)=datamatrix_new.moddata(num_f,day,cur_reg,mod)+datamatrix_tmp.moddata(num_f,day,mod,i,j);
                    end
                    tot_area(cur_reg)=tot_area(cur_reg)+1;
                end
            end
        end

        % 5. Complite averaging, deviding by the region area

        for i=1:newregs_num
            if tot_area(i)>=area_minimum
                datamatrix_new.refdata(num_f,day,i)=datamatrix_new.refdata(num_f,day,i)/tot_area(i);
                if isfield(datamatrix_tmp, 'valdata')
                    datamatrix_new.valdata(num_f,day,i)=datamatrix_new.valdata(num_f,day,i)/tot_area(i);
                end
                if isfield(datamatrix_tmp, 'constrain')
                    for mod=1:size(datamatrix_new.constrain,4) 
                        datamatrix_new.constrain(num_f,day,i,mod)=datamatrix_new.constrain(num_f,day,i,mod)/tot_area(i);
                    end
                end
                datamatrix_new.obsdata(num_f,day,i)=datamatrix_new.obsdata(num_f,day,i)/tot_area(i);
                for mod=1:size(datamatrix_new.moddata,4) 
                    datamatrix_new.moddata(num_f,day,i,mod)=datamatrix_new.moddata(num_f,day,i,mod)/tot_area(i);
                end
            else
                datamatrix_new.refdata(num_f,day,i)=NaN;
                if isfield(datamatrix_tmp, 'valdata')
                    datamatrix_new.valdata(num_f,day,i)=NaN;
                end
                if isfield(datamatrix_tmp, 'constrain')
                    for mod=1:size(datamatrix_new.constrain,4) 
                        datamatrix_new.constrain(num_f,day,i,mod)=NaN;
                    end
                end
                datamatrix_new.obsdata(num_f,day,i)=NaN;
                for mod=1:size(datamatrix_new.moddata,4) 
                    datamatrix_new.moddata(num_f,day,i,mod)=NaN;
                end
            end
        end
    end
end
