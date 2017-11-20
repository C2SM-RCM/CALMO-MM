function [VectorValues]=expval_inter_combs(temp_inter,range,default,paramn)

% NAME 
%   expval_inter_combs
% PURPOSE 
%   fill expval matrix for sims_def.m
% INPUTS 
%   temp_inter - one of the interaction simulations
%   range - Parameters ranges (min and max)
%   default - Parameters defaults
%   paramn - Parameter names
% OUTPUTS 
%   expval matrix for sims_def.m
% AUTHOR  
%   Itsik Carmona (carmonai@ims.gov.il)

          VectorValues=[];
          iend=length(range);
          if(strcmp(temp_inter,'LTKHMLTUR')==1)
            paramn1='tkhmin';
            paramn2='turl';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'LTKHMHTUR'))
            paramn1='tkhmin';
            paramn2='turl';
            max1=1;
            max2=2;
          elseif(strcmp(temp_inter,'LTKHMLENTR'))
            paramn1='tkhmin';
            paramn2='entrsc';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'LTKHMLCSOI'))
            paramn1='tkhmin';
            paramn2='csoil';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'LTKHMHCSOI'))
            paramn1='tkhmin';
            paramn2='csoil';
            max1=1;
            max2=2;
          elseif(strcmp(temp_inter,'LCRSMHTKH'))
            paramn1='crsmin';
            paramn2='tkhmin';
            max1=1;
            max2=2;
          elseif(strcmp(temp_inter,'LCRSMLTKH'))
            paramn1='crsmin';
            paramn2='tkhmin';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'LTURLENTR'))
            paramn1='turl';
            paramn2='entrsc';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'LTURLCSOI'))
            paramn1='turl';
            paramn2='csoil';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'HTURLCSOI'))
            paramn1='turl';
            paramn2='csoil';
            max1=2;
            max2=1;
          elseif(strcmp(temp_inter,'LCRSMHTUR'))
            paramn1='crsmin';
            paramn2='turl';
            max1=1;
            max2=2;
          elseif(strcmp(temp_inter,'LCSOILENTR'))
            paramn1='csoil';
            paramn2='entrsc';
            max1=1;
            max2=1;
          elseif(strcmp(temp_inter,'LCRSMHENTR'))
            paramn1='crsmin';       
            paramn2='entrsc';
            max1=1;
            max2=2;
          elseif(strcmp(temp_inter,'LCRSMLCSOI'))
            paramn1='crsmin';
            paramn2='csoil';
            max1=1;
            max2=1;
          end
          for i=1:iend
            if(strcmp(paramn1,paramn{i}))
          VectorValues=[VectorValues,range{i}(max1)];
            elseif(strcmp(paramn2,paramn{i}))
          VectorValues=[VectorValues,range{i}(max2)];
            else
          VectorValues=[VectorValues,default{1}(i)];
            end
          end
