function [area] = regions_bmp(lat,lon,img,unify_regions)

% NAME 
%   regions_bmp
% PURPOSE 
%   Definition of regions over Switzerland and north Italy
% METHOD
%   Analyze image file (regions_italy_swiss_for_matlab.bmp) where each region has its color
% INPUTS 
%   lat - latitude
%   lon - longitude
%   img - image file (regions_italy_swiss_for_matlab.bmp) where each region has its color
%   unify_regions - array that defines which regions (out of 1-7) to unify (option to unify several regions into one bigger)
% OUTPUTS 
%   area - region number to which lat lon belong
% AUTHOR  
%   Itsik Carmona (carmonai@ims.gov.il)


%imagesc(img)
ysize=size(img,1);
xsize=size(img,2);
LAT1=round(1+(48-lat)*(ysize-1)/6);
LON1=round(1+(lon-5)*(xsize-1)/9);
area=8; % (no area at all) 8 it is outof border
if((LAT1>1 && LAT1<size(img,1)) && (LON1>1 && LON1<size(img,2)))
    
    if(img(LAT1,LON1)==113); area=1 ; end % the green area
    if(img(LAT1,LON1)==87); area=2 ; end % the light brown area
    if(img(LAT1,LON1)==251); area=3; end % the yellow area
    if(img(LAT1,LON1)==253); area=4; end % The dark strong pink area
    if(img(LAT1,LON1)==231); area=5; end % The light pink area
    if(img(LAT1,LON1)==232); area=6; end % The blue area
    if(img(LAT1,LON1)==1); area=7; end % The strong brown area
    if(img(LAT1,LON1)==0) % The line area can be in between two areas that smaller than 8. For instance between area 2 to 3. It is not auotmatically area 8 when img=0 [line(black color)].
        if(LAT1+1>1 && LON1+1>1)  % anywave the area is equal 8 in the boundaries of the map.
            totalcolor=[img(LAT1+1,LON1),img(LAT1,LON1-1),img(LAT1-1,LON1),img(LAT1,LON1+1)]; % the order of the 4 grid points is: north, west, south, east.
            kcolor=find(totalcolor>0 & totalcolor<255); % if one grid point or more grid points out of the 4 points that sournding the line grid point (color 0) around the is not in area 8 than find it.
            if(~isempty(kcolor)) % if there is a color who is not black (color=0) or white (color=255) than choose the area by the fi
                color1=kcolor(1);
                if(totalcolor(color1)==113); area=1 ; end % the green area
                if(totalcolor(color1)==87); area=2 ; end % the light brown area
                if(totalcolor(color1)==251); area=3; end % the yellow area
                if(totalcolor(color1)==253); area=4; end % The dark strong pink area
                if(totalcolor(color1)==231); area=5; end % The light pink area
                if(totalcolor(color1)==232); area=6; end % The blue area
                if(totalcolor(color1)==1); area=7; end % The strong brown area
            end
        end  
    end
end

if length(unify_regions)>1
    if (~isempty(find(unify_regions==area, 1)))
        area=unify_regions(1);
    end
end
    