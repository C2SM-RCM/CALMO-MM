function datatemp=del_bad_sounding_obs(datamatrix_so,ml,sims_reg,sims_inter,sims_con,simval)
 
reg=0;
for i=1:size(datamatrix_so.obsdata,3)   % loop over "regions"
    good_fields=0;
    for jj=1:size(datamatrix_so.obsdata,1)  % loop over sounding fields
         if(length(find(~isnan(datamatrix_so.obsdata(jj,:,i))))>=ml)
                  good_fields=good_fields+1;    %count number of "good" fields
         end
    end
    
    good_fields
    size(datamatrix_so.obsdata,1)
    if (good_fields==size(datamatrix_so.obsdata,1))
        reg=reg+1;
        datatemp.obsdata(:,:,reg)=datamatrix_so.obsdata(:,:,i);
        datatemp.refdata(:,:,reg)=datamatrix_so.refdata(:,:,i);
        if ~strcmp(sims_inter,'MISSING') | ~strcmp(sims_reg,'MISSING')
            datatemp.moddata(:,:,reg,:)=datamatrix_so.moddata(:,:,i,:);
        end
        if ~strcmp(sims_con,'MISSING')
            datatemp.constrain(:,:,reg,:)=datamatrix_so.constrain(:,:,i,:);
        end
        if ~strcmp(simval,'MISSING')
            datatemp.val(:,:,reg)=datamatrix_so.val(:,:,i);
        end
    end
    
end
        
        
