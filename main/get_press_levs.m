function [OUTPUT]=get_press_levs(M,pmin,pmax,dp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation step for the model (int2) is for the
% simulation and int1 is for observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P2=pmax:-dp:pmin';

if(size(M,1)>7.99999)
    height1=interp1(M(:,1),M(:,2),P2);
    tempint1=interp1(M(:,1),M(:,3),P2);
    rhint1=interp1(M(:,1),M(:,4),P2);
    qint1=interp1(M(:,1),M(:,8),P2);
    KK=find((M(:,5)<100 & M(:,5)>=-0.000000001) & (M(:,end)<360.0001 & M(:,end)>-0.001)); % possible wind speed and direction values because there are values with WS=10^7 unreal values
    if(length(KK)>7.99999) % at list 8 levels exist
        v1=-M(KK,5).*cos(M(KK,end)*pi/180); % the North-South wind vector v (if v is positive than the wind ahead to the North, southerlies winds
        u1=-M(KK,5).*sin(M(KK,end)*pi/180); % the westrlies winds is positve  u
        vint1=interp1(M(KK,1),v1,P2);% interpolated v1
        uint1=interp1(M(KK,1),u1,P2); %interpolated u1
    else
        vint1=zeros(length(P2),1)-999.9;
        uint1=zeros(length(P2),1)-999.9;
    end
    k1=find(P2==850);
    k2=find(P2==700);
    k3=find(P2==500);
    OUTPUT=[tempint1(k1),tempint1(k2),tempint1(k3),rhint1(k1),rhint1(k2),rhint1(k3),uint1(k1),uint1(k2),uint1(k3),vint1(k1),vint1(k2),vint1(k3)];
else
    OUTPUT=zeros(1,12)-999.9;
end


                                                    
   
