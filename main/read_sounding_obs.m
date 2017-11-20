function output=read_sounding_obs(year,month,day,height_step,windshear,windshearopt)

% NAME
%   read_sounding_obs
% PURPOSE
%   Read and interpolate soundings data, calculate sounding characteristics  
% INPUTS
%   year
%   month
%   day
%   height_step - The interpolation height step in meters
%   windshear - [v1u,v1d,v2u,v2d,v3u,v3d] - pressure levels for calculating wind shears:
%       v1u - the upper pressure level for wshear1
%       v1d - the bottom pressure level for wshear 1
%       v2u - as mentioned above but for wshear 2
%       v2d - as mentioned abobe but for wshear 2
%       v3u - as mentioned above but for wshear 3
%       v3d - as mentioned above but for wshear 3, whereas 1100 is the surface level or the lowest level ("below surface")
%   windshearopt - windshear calculation method: 'scalar' or 'vector'
% OUTPUTS
%   output - Data matrix with dimensions [Field,Day,1,Hour,Sounding location]
% AUTHOR
%   Itsik Carmona (carmonai@ims.gov.il)



global obsdir extdir

soundingcord=load([extdir '/radiosondes_metadata.txt']);    % 'soundingcord' = stations number, lat and lon
station_list=soundingcord(:,3);
Mread=[]; % The total file matrix with vertical profile data
for st_num=1:length(station_list)
    for h_num=1:8
        Mread{st_num,h_num}=[]; % The total file matrix with vertical profile data
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME OF PREDICTION
time_thresh=1.25/24; % if the time mistake is more than 1h+15min, cancel the sounding
yyyy_mod=num2str(year);
mm_mod=num2str(month); if month<10 mm_mod=['0' mm_mod]; end
dd_mod=num2str(day); if day<10 dd_mod=['0' dd_mod]; end

for h_num=1:8
    hourcalc=(h_num-1)*3;
    hour_2digits=num2str(hourcalc); if hourcalc<10 hour_2digits=['0' num2str(hourcalc)]; end % hour_2digits is the string variable of the hour with two digits
    obsdate_list(h_num)=datenum([yyyy_mod mm_mod dd_mod hour_2digits],'yyyymmddHH');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SONDEDIR=[obsdir '/soundings']; % The directories of the radiosondes
ml=5; % the minimum number of levels in sounding
linestop=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of the output 3D matrix
output=zeros(18,8,length(station_list))-999.9; % definition of OUTPUT 3D MATRIX
datechar2=[yyyy_mod mm_mod dd_mod];
Aopen=[SONDEDIR,'/radiosondes_',datechar2,'.txt'];
fclose('all');
FID=fopen(Aopen,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rd=287.058; % dry air gas constant (Rd) [J/(kg*K)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rv=461.5;   % wate{nstation_place(1),nhour_place(1)}r vapur gas constant (Rv) [J/(kg*K)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsi=Rd/Rv; % The ratio between Rd to Rv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DP=10;      % DP is the layer pressure thickness [mb]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.8076;   % Gravity Accelration [m/s^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cpd=1003.5; % Cpd is the specific heat for dry air [J/(kg*K)]
gammadry=g/Cpd;     % 'Lambada' is the temperature lapse rate [K/m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File example:
%    wmo_id           termin  latitude longitude      elev   stn_name              int_ind  nat_abbr  track_type   prof_type   type      level           742           745           747           748           743
%
%        70 | 20130429000000 |  36.683 |   3.217 |      25 | DAR-EL-BEIDA        |   60390 |         |         8 |        57 |    1 |  1006.00 | 10000000.00 |       13.40 |       11.60 |        8.70 |      310.00

%flag_found=0;   %flag that turns to =1 only when the profile for the input station and date is found in the file
flag_next=0;
flag_skip=0;
if exist(Aopen,'file')  % if sounding file exists
    line1=fgets(FID); clear line1;  % skip 3 header lines
    line2=fgets(FID); clear line2;  % skip 3 header lines
    line3=fgets(FID); clear line3;  % skip 3 header lines
    
    while (linestop==0)                             % go over the file
        if flag_next==0 || flag_skip==1
            line1=fgets(FID);
        end
        if length(line1)>20                         % the line shorter than 20 chars is end of file
            lat1=str2num(line1(30:38));             % read lat
            lon1=str2num(line1(40:48));             % read lon
            
            if ~isempty(str2num(line1(82:90)))
                nstation=str2num(line1(82:90));         % read station number
                flag_skip=0;
            else
                flag_skip=1;                            % when current station lines became undef, flag_skip=1
                nstation=-999.9;
            end
                
            year1=str2num(line1(14:17));            % read year
            month1=str2num(line1(18:19));           % read month
            day1=str2num(line1(20:21));             % read day
            hour1=str2num(line1(22:23));            % read hour
            min1=str2num(line1(24:25));             % read minute
            % reading sounding time:
            yyyy_obs=num2str(year1);
            mm_obs=num2str(month1); if month1<10 mm_obs=['0' mm_obs]; end
            dd_obs=num2str(day1); if day1<10 dd_obs=['0' dd_obs]; end
            hh_obs=num2str(hour1); if hour1<10 hh_obs=['0' hh_obs]; end
            obsdate_tmp=datenum([yyyy_obs mm_obs dd_obs hh_obs],'yyyymmddHH')+min1/1440; %ex: 2013010100 % *also adding the minutes
            time_diff=abs(obsdate_list-obsdate_tmp);
            
            flag_next=0;        % when current station read finished, flag_next=1
            
            if ~isempty(find(station_list==nstation)) && ~isempty(find(time_diff<=time_thresh))
                nstation_place=find(station_list==nstation);
                nhour_place=find(time_diff<=time_thresh);
                
                nstation_tmp=nstation;
                
                while (flag_next==0) && (flag_skip==0)      % go over one specific profile    
                    temp=str2num(line1(158:170));           % read temp.
                    tdew=str2num(line1(172:184));           % read dew temp.
                    pressobs=str2num(line1(133:142));       % read pressure
                    ws=str2num(line1(186:198));            % read wind speed
                    wdir=str2num(line1(200:212));           % read wind dir.
                    height=str2num(line1(144:156));         % read level height
                    L=2500800-2360*temp+1.6*(temp^2)-0.06*(temp^3);     % latent heat L (func. of temp.), assuming -25<T<40
                    es=611.3*exp((L/Rv)*(1/273.15-1/(temp+273.15)));    % water vapor pressure at saturation [Pascal]
                    e=611.3*exp((L/Rv)*(1/273.15-1/(tdew+273.15)));     % water vapor pressure [Pascal]
                    rh=100*(e/es);                                      % rel. humidity
                    r=e*epsi/(pressobs*100-e);                          % mixing ratio = water vapur mass / dry air mass [kg/kg]
                    q=r/(1-r);                                          % specific humidity = water vapur mass / total air mass [kg/kg]
                    rsatu=(es/(pressobs*100))*epsi;                     % Saturated ixing ratio = Saturated water vapur mass / dry air mass [kg/kg]
                    gammaw=g*(1+L*rsatu/(Rd*(temp+273.15)))/(Cpd+(L^2)*rsatu*epsi/(Rd*(temp+273.15)^2));        % The wet lapse rate [K/m]
                    vapdens1=((rh/100)*es)/(Rv*(temp+273.15));                                          % Absolute humidity (kg/m^3), using T,RH
                    vapdens2=((q/(1-q))*pressobs*100/(Rd*(temp+273.15)))*(1+(q/(1-q))*Rv/Rd);           % Absolute humidity (kg/m^3), using T,Tdew
                    ro=(pressobs*100-e+pressobs*100*r)/(Rd*(temp+273.15));               % air density for wet idealized gas: (P-e+Pr)/RT
                    
                    % read line, even if the height is undefined
                    
                    if ((temp<=50 && temp>=-80) && rh<=102 && q<=0.1 && r<=0.1 && (tdew<=temp) && tdew>=-150)
                        Mread{nstation_place(1),nhour_place(1)}=[Mread{nstation_place(1),nhour_place(1)};pressobs,height,temp,rh,ws,vapdens1,vapdens2,q,r,gammaw,ro,lat1,lon1,nstation,tdew,year1,month1,day1,hour1,wdir];
                    end
                    line1=fgets(FID);
                    if length(line1)>20                             % the line shorter than 20 chars is end of file         
                        if ~isempty(str2num(line1(82:90)))
                            nstation_tmp=str2num(line1(82:90));     % read station number
                            if (nstation_tmp~=nstation)
                                flag_next=1;                        % station changed --> exit while over giving station
                            end
                        else
                            flag_skip=1;                            % station number missing --> exit while over giving station
                        end
                    else
                        flag_next=1;                                % reached end of file --> exit while over giving station
                        linestop=1;
                    end
                end
            end
            
            
        else
            linestop=1;     % end of file
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % At this stage matrix Mread contains all the levels data for the input station and time, before missing heights definitions, before interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The next stage contatins 2 loops for the stations and hours in order to calcuate the sounding fileds if they exist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for st_num=1:length(station_list)
        for h_num=1:8
            data_exist=1; % if data is missing for this station then data_exist=0, if not it remains 1.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Next, we replace undef heights by calculated ones according temperature and hydrostatic approx.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear M;
            M=Mread{st_num,h_num};   %  assign the Mread matrix (the total file read Matrix to the loop Matrix named "M"
            if isempty(M)
               output(:,h_num,st_num)=-999.9;
               data_exist=0;
            end
            if(data_exist==1) % because of data_exist==0 than matrix M is empty and the next condition term can be treated well
                if length(find(M(:,2)<50000))<ml
                   output(:,h_num,st_num)=-999.9;
                   data_exist=0;
                end
            end
            if(data_exist==1) % if data_exist=1 (there are more than "ml" levels with data including height data) than it starts to cacluate the sounding meteorlogical fileds such as CAPE and etc'
                for i=1:size(M,1)
                    if M(i,2)>50000                 % if the level height is undef, 2 options:
                        if(i==1)                            % case when the first line is undef
                            iexist=find(M(:,2)<50000);   % find all defined levels above the current one
                            if(length(iexist)>=1)
                                for ii2=iexist(1):-1:2      % iexist(1) is the lowest level where height is defined. Calculating Height climbing down: for i=iexist(1)-1,iexist(1)-2,...
                                    ii1=ii2-1;
                                    M(ii1,2)=M(ii2,2)+(M(ii2,1)-M(ii1,1))*100/(g*(M(ii2,11)+M(ii1,11))/2);
                                end
                            else
                                display('No height information within the sounding data. Setting output=-999.9');
                                output(:,h_num,st_num)=-999.9;
                                return
                            end
                        else                                % case when first line exists
                            M(i,2)=M(i-1,2)-(M(i,1)-M(i-1,1))*100/(g*((M(i,11)+M(i-1,11))/2)); % assign the height in the M matrix = Z1-(P2-P1)/(g*ro)
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % delete all non monotonic heights :
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Mtmp=M; clear M;
                M=Mtmp(1,:);
                i_good=1;
                for i=1:size(Mtmp,1)-1
                    if((Mtmp(i+1,2)-M(i_good,2))>=0.1)  % if the height is increasing at least at 10 cm
                        i_good=i_good+1;
                        M=[M;Mtmp(i+1,:)];
                    end
                end
                % Mp is a matrix with pressure that is monotonic with height and no constant , it is necessary for interpolation needs.
                Mp=[];
                Mp=M(1,:);
                i_good=1;
                for i=1:size(M,1)-1
                    if((M(i+1,1)-Mp(i_good,1))<=-0.1)  % if the pressyre is decreasing at least at -0.1 mb
                        i_good=i_good+1;
                        Mp=[Mp;M(i+1,:)];
                    end
                end
                pmin=50; % min pressure for interpolation
                pmax=1000; %max pressure for interpolation
                dp=50; %pressure intervals for interpolation
                extrafields=get_press_levs(Mp,pmin,pmax,dp);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % At this stage matrix M contains all the levels data for the input station and time, before interpolation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % interpolation step :
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                H11=ceil(M(1,2));
                H22=floor(max(M(:,2)));
                DH=H11:height_step:H22;         % interpolated heights
                DH=DH';
                
                pressint=interp1(M(:,2),M(:,1),DH,'pchip');                     % pressure interpolation
                tempint=interp1(M(:,2),M(:,3),DH,'pchip');                      % temp. interpolation
                tempdint=interp1(M(:,2),M(:,15),DH,'pchip');                    % Tdew interpolation
                Lint=2500800-2360*tempint+1.6*(tempint.^2)-0.06*(tempint.^3);   % latent heat L interpolation
                eint=611.3*exp((Lint./Rv).*(1/273.15-1./(tempdint+273.15)));    % water vapor pressure interpolation
                esint=611.3*exp((Lint./Rv).*(1/273.15-1./(tempint+273.15)));    % saturated water vapor pressure interpolation
                rhint=100*(eint./esint);                                      % rel. humidity interpolation
                %rhint=100*eint./esint;
                rint=eint.*epsi./(pressint*100-eint);                          % mixing ratio interpolation = water vapur mass / dry air mass [kg/kg]
                qint=rint./(1-rint);                                          % specific humidity interpolation = water vapur mass / total air mass [kg/kg]
                vapdens1int=((rhint./100).*esint)./(Rv.*(tempint+273.15));                                                % Absolute humidity interpolation (kg/m^3), using T,RH
                vapdens2int=((qint./(1-qint)).*pressint.*100./(Rd*(tempint+273.15))).*(1+(qint./(1-qint)).*Rv/Rd);           % Absolute humidity interpolation (kg/m^3), using T,Tdew
                Mwind=[];
                for iwind=1:size(M,1) % build Mwind matrix which is clear from wind vecloity anavliable data
                    if((M(iwind,5)>=0 && M(iwind,5)<=100) && (M(iwind,20)>=0 && M(iwind,20)<=360))
                        Mwind=[Mwind;M(iwind,2),M(iwind,5),-M(iwind,5)*sin(M(iwind,20)*pi/180),-M(iwind,5)*cos(M(iwind,20)*pi/180)];   % the matrix fileds are: height,total horizantal wind, u and v
                    end
                end 
                if size(Mwind,1)>=ml % if the number of lines (levels) in WIND matrix is equal of above "ml"                  
                     wsint=interp1(Mwind(:,1),Mwind(:,2),DH,'pchip');                     % wind speed interpolation by pressure that was already interpolted before
                     uint=interp1(Mwind(:,1),Mwind(:,3),DH,'pchip');
                     vint=interp1(Mwind(:,1),Mwind(:,4),DH,'pchip');
                else
                     wsint=zeros(size(DH,1),1)-999.9;
                     uint=zeros(size(DH,1),1)-999.9;
                     vint=zeros(size(DH,1),1)-999.9;
                end
               
                rsatu=(epsi.*esint)./(pressint.*100-esint);
                qint=(eint.*epsi)./(pressint.*100-eint*(1-epsi));
                rint=(eint.*epsi)./(pressint.*100-eint);
                gammawint=g.*(1+(Lint.*rsatu)./(Rd.*(tempint+273.15)))./(Cpd+(Lint.^2).*rsatu.*epsi./(Rd.*(tempint+273.15).^2)); % The wet lapse rate interpolation [K/m]
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % matrix Mint contains only the relvant station and relvant time with all levels data, after interpolation
                Mint=[pressint,DH,tempint,rhint,wsint,vapdens1int,vapdens2int,qint,rint,gammawint,Lint,eint,esint,uint,vint];
                
                nline1=size(Mint,1);    %number of lines
                tempgrad=zeros(nline1,1);
                for ir=1:nline1-1
                    tempgrad(ir)=(Mint(ir,3)-Mint(ir+1,3))./height_step;
                end
                tempgrad(nline1)=tempgrad(nline1-1);
                Mint=[Mint,tempgrad]; % adding temperature gradient to Mint !
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Mparcel=[];
                temp1=Mint(1,3);
                pressobs1=Mint(1,1);
                staheight=Mint(1,2);    % sounding first level
                tempmass=temp1;         % tempmass is the air parcel temperature [C]
                tempvir=(Mint(:,3)+273.15).*(Mint(:,9)+epsi)./(epsi*(1+Mint(:,9)))-273.15; % The virtual air profile temperature [C]
                %       Mparcel - the path of the air parcel during its rising with the follwing information:
                %       air pressure, height above msl, parcel temp., surrounding temp., parcel virtual temp., surrounding virtual temp.
                Mparcel=[pressobs1,staheight,temp1,temp1,tempvir(1),tempvir(1)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step1: forced dry lapse rate until saturation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                dryoff=0; % dryoff flag is 1 when the air parcel is in saturation
                hrun=Mint(1,2);
                istep=1;
                rfirst=(eint(1)*epsi)./(Mint(1,1).*100-eint(1));  % mixing ratio is constant during dry lapse
                erun=eint(1);                                   % vapor pressure
                esrun=611.3.*exp((L(1)./Rv)*(1/273.15-1./(tempmass+273.15))); % saturation vapor pressure
                Mdry=Mparcel;% The dry path of the air parcel
                while(esrun>=erun && istep<=(nline1-1)) % dry lapse, which occurs untill saturation or max hight
                    hrun=hrun+height_step;
                    tempmass=tempmass-height_step*gammadry;
                    istep=istep+1;
                    erun=(rfirst*Mint(istep,1)*100)/(epsi+rfirst); % running parcel vapor pressure, note: rfirst=const
                    Lrun=2500800-2360*tempmass+1.6*(tempmass.^2)-0.06*(tempmass.^3); % Lrun is the running Latent heat of the parcel
                    esrun=611.3*exp((Lrun/Rv)*(1/273.15-1/(tempmass+273.15))); % saturated running parcel vapor pressure
                    tempmassvir=(tempmass+273.15)*(rfirst+epsi)/(epsi*(1+rfirst))-273.15; % running virt. temp. of parcel
                    Mparcel=[Mparcel;Mint(istep,1),Mint(istep,2),tempmass,Mint(istep,3),tempmassvir,tempvir(istep)];
                    Mdry=[Mdry;Mint(istep,1),Mint(istep,2),tempmass,Mint(istep,3),tempmassvir,tempvir(istep)];
                end
                
                Heightf=[];             % important heights matrix
                if(istep<nline1)       % reached 100% (before max height)
                    dryoff=1;
                    Hlcl=Mint(istep,2); % Hlcl - height where air parcel RH starts to be 100%
                    Heightf=[Heightf,Hlcl];
                else                    % no saturation happened until the max height
                    Heightf=[Heightf,-999.9];
                    %Mdry=[-999.9,-999.9,-999.9,-999.9,-999.9,-999.9];
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step1 END: Finished forced dry lapse, starting wet lapse
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step2: Starting wet lapse.
                %        The air parcel contiunes to rise in 100% with wet lapse rate which depends on the mixing ratio
                %        and temp. (through L) of the air parcel.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if(dryoff==1)               % general case (happens always), when the parcel reaches saturation at some level
                    
                    Mwet=Mdry(istep,:);     % the first line of Mwet parcel data is the last Mdry parcel line
                    for ii=istep+1:nline1-1   % running with wet lapse untill max height (!!!)
                        Lrun=2500800-2360*tempmass+1.6*(tempmass.^2)-0.06*(tempmass.^3);  % Lrun is the running Latent heat of the parcel
                        esrun=611.3*exp((Lrun/Rv)*(1/273.15-1/(tempmass+273.15)));        % = erun = saturated running parcel vapor pressure
                        rrun=(epsi*esrun)/(Mint(ii-1,1)*100-esrun);                         % mixing ratio, which is no more constant, now it is decreasing with height because of condensation
                        gammawrun=g*(1+(Lrun*rrun)/(Rd*(tempmass+273.15)))/(Cpd+(Lrun^2)*rrun*epsi/(Rd*(tempmass+273.15).^2));  % wet lapse rate is not constant
                        tempmass=tempmass-gammawrun*height_step;                          % running temp. of parcel
                        tempmassvir=(tempmass+273.15)*(rrun+epsi)/(epsi*(1+rrun))-273.15; % running virt. temp. of parcel
                        Mparcel=[Mparcel;Mint(ii,1),Mint(ii,2),tempmass,Mint(ii,3),tempmassvir,tempvir(ii)];
                        Mwet=[Mwet;Mint(ii,1),Mint(ii,2),tempmass,Mint(ii,3),tempmassvir,tempvir(ii)];
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % In the following: Kwetold gets all the levels where parcel virt. temp > env. virt. temp.
                    %                   Kwet gets only the first consecutive levels where parcel virt. temp > env. virt. temp.
                    Kwet=find(Mwet(:,5)>Mwet(:,6));             % find all the levels where parcel virt. temp > env. virt. temp.
                    if(length(Kwet)>=2)                         % if found at least 2 levels where parcel virt. temp > env. virt. temp.
                        %firstflag=0;
                        Kwetold=Kwet;
                        for jj=1:(length(Kwet)-1)
                            %if(Kwetd>1.999 && firstflag==0)
                            if (Kwetold(jj+1)-Kwetold(jj))>1.999
                                Kwet=Kwetold(1:jj);
                                %firstflag=1;
                                break;
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % CAPE & CIN calculation :
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Mcape=[];
                    Mcin=[];
                    capeflag=0;
                    
                    if(length(Kwet)>=1)                         % case when there is CAPE
                        Hlfc=Mwet(Kwet(1),2);                   % Hlfc = height where the wet (staurated) air parcel starts to be warmer than the air around
                        Hel=Mwet(Kwet(length(Kwet)),2);         % Hel = TOP height where the wet (staurated) air parcel is still warmer than the air around
                        Heightf=[Heightf;Hlfc;Hel];
                        
                        if(length(Kwet)>=3)
                            Mcape=Mwet(Kwet,:);
                            capeflag=1;
                        else
                            capeflag=0;
                        end
                        
                        Kcinmax=Kwet(1)-1;
                        if(Kcinmax>=1)                          % Case when cin continues above dry lapse region
                            Mcin=[Mdry;Mwet(1:Kcinmax,:)];
                        else
                            Mcin=Mdry;                          % Case when cin is only at the dry lapse region
                        end
                    else                                        % Case when there is NO CAPE
                        Mcin=[Mdry;Mwet];                       % In this case, cin goes until the max height
                    end
                else
                    Mcin=Mdry;                                  % Unprobable case, when the parcel doesn't reach saturation at any level. Then cin is just all the dry lapse route
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%
                % CAPE CALCULATION :
                %%%%%%%%%%%%%%%%%%%%%
                cape=0;
                if(capeflag==1)
                    for i=1:(size(Mcape,1)-1)
                        i1=i;
                        i2=i+1;
                        cape=cape+g*(Mcape(i2,2)-Mcape(i1,2))*((Mcape(i1,5)+Mcape(i2,5))/2-(Mcape(i1,6)+Mcape(i2,6))/2)/((Mcape(i1,6)+Mcape(i2,6))/2+273.15);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%
                % CIN CALCULATION :
                %%%%%%%%%%%%%%%%%%%%%
                cin=0;
                for i=1:(size(Mcin,1)-1)
                    i1=i;
                    i2=i+1;
                    cinadd=g*(Mcin(i2,2)-Mcin(i1,2))*((Mcin(i1,5)+Mcin(i2,5))/2-(Mcin(i1,6)+Mcin(i2,6))/2)/((Mcin(i1,6)+Mcin(i2,6))/2+273.15);
                    if(cinadd<=0)           % add to the cin only the negative values (ignore superadiabatic region)
                        cin=cin+cinadd;
                    end
                end
                cin=-cin;
                
                %%%%%%%%%%%%%%%%%%%%%
                % Wind Shear Calculation :
                %%%%%%%%%%%%%%%%%%%%% 
                wshearout=zeros(length(windshear)/2,1)-999.9;% realoction of 3 windshears (The difinitation of the vector
                    for jjj=1:length(windshear)/2               % the loop for 3 windshears 
                        jjj1=jjj*2;                             % for the bottom level
                        jjj2=jjj*2-1;                           % for the upper level
                        if (windshear(jjj1)~=1100) 
                            [mdiff ku]=min(abs(Mint(:,1)-windshear(jjj2))); % find the upper level 
                            [mdiff kd]=min(abs(Mint(:,1)-windshear(jjj1))); % find the bottom level
                        else
                            [mdiff ku]=min(abs(Mint(:,1)-windshear(jjj2)));
                            kd=1; % the lowest level of Mint matrix (shoud be the surface level or at least the closest to the surface level if surface level does not exist)
                        end
                        if ((length(ku)==1 && length(kd)==1) && abs((Mint(ku,1)-windshear(jjj2))/windshear(jjj2))<0.01 ...
                        	&& (windshear(jjj1)==1100 || abs((Mint(kd,1)-windshear(jjj1))/windshear(jjj1))<0.01))   % if the pressure level difference between the closest interpolated level to the relevnt level is above 1% than don't calcute the windshear
                   % if the bottom windshear==1100 than there is no problem
                   % because it is the lowest or the surface level
                   % otherwise the difference should be less than 1% , the
                   % same condition that was done for the upper level is for the
                   % bottom level.
                            if ((Mint(ku,2)-Mint(kd,2))>0 && (Mint(ku,5)<100 && Mint(ku,5)>=0) && (Mint(kd,5)<100 && Mint(kd,5)>=0)) % no wind above 100m/s it is impossible and also the difference between the upper height to the lower height should be above zero
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wind shear scalaric calculation
                                if strcmp(windshearopt,'scalar')                % the windshear cacluation by the scalar options 
                                    wshearout(jjj)=(Mint(ku,5)-Mint(kd,5))/(Mint(ku,2)-Mint(kd,2)); % the windshear output
                                elseif strcmp(windshearopt,'vector')            % VECTORIC CACLUATION OF WIND SHEAR
                                    wshearout(jjj)=(sqrt((Mint(ku,14)-Mint(kd,14))^2+(Mint(ku,15)-Mint(kd,15))^2))/(Mint(ku,2)-Mint(kd,2)); % the windshear output
                                end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            end
                        end
                    end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % END of CAPE , CIN calculation & WindShear calculation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TCWC calculation :
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                TCWC1=0;
                TCWC2=0;
                for i=1:(size(Mint,1)-1)
                    TCWC1=TCWC1+height_step*(Mint(i,6)+Mint(i+1,6))/2;      % integral over the absolute humidity (kg/m^3), using Tdew
                    TCWC2=TCWC2+height_step*(Mint(i,7)+Mint(i+1,7))/2;      % integral over the absolute humidity (kg/m^3), using q
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                output(1:18,h_num,st_num)=[cape,cin,TCWC1,wshearout(1:3)',extrafields];
                % end
            end % if of data_exist==1 (if there are enough measured height levels (Above or equal "ml" levels)
        end % if of h_num loop
    end % if of st_num loop
else % else for if the file exist
    display(['The sounding file ' filestr ' does not exist']);
end % end for if the file exist

    