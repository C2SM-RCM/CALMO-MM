function [datamatrix_so]=corr_so(datamatrix_so,bad_list,vars_sound,sims_reg,sims_inter,sims_con,simval)

for iday=1:length(bad_list.day)
    for ivar=1:length(bad_list.vars);
        varindex=find(strcmp(bad_list.vars(ivar),vars_sound));
        datamatrix_so.obsdata(varindex,bad_list.day(iday),:)=NaN;
        datamatrix_so.refdata(varindex,bad_list.day(iday),:)=NaN;
        if ~strcmp(sims_inter,'MISSING') | ~strcmp(sims_reg,'MISSING')
            datamatrix_so.moddata(varindex,bad_list.day(iday),:,:)=NaN;
        end
        if ~strcmp(sims_con,'MISSING')
            datamatrix_so.constrain(varindex,bad_list.day(iday),:,:)=NaN;
        end
        if ~strcmp(simval,'MISSING')
            datamatrix_so.val(varindex,bad_list.day(iday),:,:)=NaN;
        end
    end
end