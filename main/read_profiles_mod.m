function output = read_profiles_mod(simdir,year,month,day,simtype,height_step,sound_exist,windshear,windshearopt)

% NAME
%   read_profiles_mod
% PURPOSE
%   Read and interpolate profiles data, calculate profiles characteristics  
% INPUTS
%   simdir - path to simulations files 
%   year
%   month
%   day
%   simtype - simulation name
%   height_step - The interpolation height step in meters
%   sound_exist - binary matrix [Day,Hour,Sounding location] with ones where the sounding data exist
%   windshear - [v1u,v1d,v2u,v2d,v3u,v3d] - pressure levels for calculating wind shears:
%       v1u - the upper pressure level for wshear1
%       v1d - the bottom pressure level for wshear 1
%       v2u - as mentioned above but for wshear 2
%       v2d - as mentioned abobe but for wshear 2
%       v3u - as mentioned above but for wshear 3
%       v3d - as mentioned above but for wshear 3, whereas 1100 is the surface level or the lowest level ("below surface")
%   windshearopt - windshear calculation method: 'scalar' or 'vector'
% OUTPUTS
%   output - Data matrix with dimensions [Field,Day,simulation,Hour,Sounding location]
% AUTHOR
%   Itsik Carmona (carmonai@ims.gov.il)

global maindir curdir extdir

soundingcord=load([extdir '/radiosondes_metadata.txt']);    % 'soundingcord' = stations number, lat and lon
station_list=soundingcord(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ml=20; % minimum number of levels in the COSMO atmospheric profile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dd=num2str(day); if day<10 dd=['0' num2str(day)]; end
mm=num2str(month); if month<10 mm=['0' num2str(month)]; end
display(['The forecasting date ' num2str(year) mm dd]); % the forecasting time
dayfile=datenum([num2str(year) mm dd '12'],'yyyymmddHH'); % The date forecasting time
daynext2=dayfile-1; % The file name is one day before !!! (or for +12h to +36h , therefore I subtracted -1 day)
monfile=str2num(datestr(daynext2,5)); % The month of in the file name
dayfile=str2num(datestr(daynext2,7)); % The day in the file name
yearfile=str2num(datestr(daynext2,10)); % The year in the file name
mmfile=num2str(monfile); if monfile<10 mmfile=['0',num2str(monfile)]; end % the string month in the file name
ddfile=num2str(dayfile); if dayfile<10 ddfile=['0',num2str(dayfile)]; end % the string day in the file name

Mread=[]; % The total file matrix with vertical profile data
for st_num=1:length(station_list)
    for h_num=1:8
        Mread{st_num,h_num}=[]; % The total file matrix with vertical profile data
    end
end

simtype=char(simtype);
%dirdate=[num2str(yearfile-2000),mmfile,ddfile,'12'];
dirdate=[num2str(yearfile-2000),mmfile,ddfile,'00'];
if strcmp(simtype,'REF')
    %fname=['LMvert_' dirdate '_DEF_20' dirdate];
    fname=['LMvert_DEF_20' dirdate];
else
    display(['reading LMvert_' simtype '_20' dirdate ' ...'])
    fname=['LMvert_' simtype '_20',dirdate];
end

fullfname=[simdir '/' fname];


fclose('all');
FID=fopen(fullfname,'r');
linestop=0;
%%%%%%%%%%%%%%%%%%%%%
%dry air gas constant (Rd) [J/(kg*K)]
Rd=287.058;
%%%%%%%%%%%%%%%%%%%%%
% water vapur gas constant (Rv) [J/(kg*K)]
Rv=461.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The ratio between Rd to Rv
epsi=Rd/Rv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DP is the layer pressure thick [mb]
DP=10;
DPS=num2str(DP);
%%%%%%%%%%%%w%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Accelration [m/s^2]
g=9.8076;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cpd is the specific heat for dry air [J/(kg*K)]
Cpd=1003.5;
% 'Lambada' is the temperature lapse rate [K/m]
gammadry=g/Cpd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=[]; Mint=[]; Mparcel=[]; Mdry=[]; Mwet=[]; Mcape=[]; Mcin=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

year2digits=year-floor(year/100)*100;
capeflag=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output=zeros(nvars,8,length(station_list))-999.9; % definiation of OUTPUT 3D MATRIX
output=zeros(18,8,length(station_list))-999.9; % definiation of OUTPUT 3D MATRIX

if exist(fullfname,'file') % if COSMO vertical file exist
    while(linestop==0) % go over the file until end of file 
        line1=fgets(FID);
        line1size=length(line1);
        if (line1size>1)    %read specific line
            typen=str2num(line1(16:17));    %typen=14,53,54 for wind speed in cm/s and =13 for m/s
            line2=fgets(FID); clear line2;
            line3=fgets(FID);   % the first needed line in file, which includes date,lat,lon, station number
            %Pavel% staheightchar=line3(9:12);
            %Pavel% staheight=str2num(staheightchar);
            %Pavel% Heightf=[staheight];
            yearn=str2num(line1(4:5));%year of the foreacst
            monthn=str2num(line1(6:7));%month of the foreacst
            dayn=str2num(line1(8:9));%day of the foreacst
            hourn=str2num(line1(10:11));%hour of the foreacst
            latsta=str2num(line3(13:17))/100;%lat of the g.p.
            lonsta=str2num(line3(18:22))/100;%lon of the g.p.
            stationnum=str2num(line3(4:8)); % Station number
            
            if ((~isempty(find(station_list==stationnum))) && (yearn==year2digits && monthn==month && dayn==day))
                
                sta_place=find(station_list==stationnum); % the place of the station in Mread matrix ;
                hour_place=round(hourn/3+1);
                
                if sound_exist(hour_place,sta_place)==1
                
                    linenum=3;
                    iend=0;
                    while(iend==0) % go over one specific profile
                        linenext=fgets(FID); %read specific line
                        
                        if ~isempty(strfind(linenext,'****'))
                            display(['Bad model profile, station num=' num2str(stationnum) ', forecast hour=' num2str(hourn) '. Exitting...']);
                            %break
                            % output=zeros(nvars,8,length(station_list))-999.9;
                            output=zeros(18,8,length(station_list))-999.9;
                            return
                        end
                        
                        linenum=linenum+1;
                        lineend=str2num(linenext(1:2));
                        if(lineend~=25) %flag 25 means line without data (3 lines afterwards will be header, or file ended)
                            if(lineend==70) %flag 70 means surface level
                                pressobs1=str2num(linenext(3:7))/10;    % pressure at surface
                                height1=str2num(linenext(8:12));        % height at surface
                                temp1=str2num(linenext(13:18))/10;      % temperature at surface
                                difftdew1=str2num(linenext(21:24))/10;  % difference T-Tdew
                                tdew1=temp1-difftdew1;                  % Tdew point
                                L1=2500800-2360*temp1+1.6*(temp1^2)-0.06*(temp1^3); % Latent heat near the surface
                                es1=611.3*exp((L1/Rv)*(1/273.15-1/(temp1+273.15))); % saturated water vapor pressure near the surface
                                e1=611.3*exp((L1/Rv)*(1/273.15-1/(tdew1+273.15)));  % water vapor pressure near the surface
                                rh1=e1/es1*100;                                     % calculated relative humdity near the surface
                                q1=(e1*epsi)./(pressobs1*100-e1*(1-epsi));          % calculated specific humidty near the surface
                                r1=q1/(1-q1);                                       % calculated mixing ration near the surface
                                ws1=str2num(linenext(30:33));                       % the wind speed (cm/s) near the surface
                                if(typen==14 || typen==53 || typen==54)
                                    ws_fact=100;
                                elseif(typen==13)
                                    ws_fact=1;
                                end
                                ws1=ws1/ws_fact;
                                wdir1=str2num(linenext(26:28));                     % the wind direction (in degree) near the surface
                                rsatau1=(es1/(pressobs1*100))*epsi;                 % is The Surface saturated Wet vapur mass devided by dry Air mass [kg/kg]L^2
                                gammaw1=g*(1+L1*rsatau1/(Rd*(temp1+273.15)))/(Cpd+(L1^2)*rsatau1*epsi/(Rd*(temp1+273.15)^2));  % the wet lapse rate near the surface
                                % 'vapdens11' is the density of the water vapour [kg/m^3] by using Tdew near the surface
                                vapdens11=((rh1/100)*es1)/(Rv*(temp1+273.15));
                                % 'vapdens21' is the density of the water vapour [kg/m^3] by using P,q (specific water) and T near the surface
                                vapdens21=((q1/(1-q1))*pressobs1*100/(Rd*(temp1+273.15)))*(1+(q1/(1-q1))*Rv/Rd);

                                Mread{sta_place,hour_place}=[Mread{sta_place,hour_place};pressobs1,height1,temp1,rh1,ws1,vapdens11,vapdens21,...
                                    q1,r1,gammaw1,latsta,lonsta,yearn,monthn,dayn,hourn,tdew1,wdir1];

                            end
                            if(lineend==0) %flag 0 means levels above surface
                                pressobs=str2num(linenext(3:7))/10;     %pressure
                                height=str2num(linenext(8:12));         %height
                                temp=str2num(linenext(13:18))/10;       %temperature
                                difftdew=str2num(linenext(21:24))/10;   %difference T-Tdew
                                tdew=temp-difftdew;                     % Tdew point
                                rh=str2num(linenext(61:70));            %rel.humidity
                                q=str2num(linenext(39:48));             %spec. humidity = (moist mass/tot mass) in kg/kg
                                % the latent heat ('L') is a function of temperature and the next line is an assumption of 'L' between -25 to +40 celcuis
                                L=2500800-2360*temp+1.6*(temp^2)-0.06*(temp^3);
                                % The water pressure at saturation 'es' is a function of temperautre the next 'es' is in unit [Pascal]
                                es=611.3*exp((L/Rv)*(1/273.15-1/(temp+273.15)));
                                % r is The wet ratio of Wet vapur mass devided by dry Air mass [kg/kg]L^2
                                r=q/(1-q);
                                rsatau=(es/(pressobs*100))*epsi; % is The saturated Wet vapur mass devided by dry Air mass [kg/kg]L^2
                                % The wet lapse rate (gamma W) [K/m]
                                gammaw=g*(1+L*rsatau/(Rd*(temp+273.15)))/(Cpd+(L^2)*rsatau*epsi/(Rd*(temp+273.15)^2));
                                % 'vapdens11' is the density of the water vapour [kg/m^3] by using Tdew near the surface
                                vapdens1=((rh/100)*es)/(Rv*(temp+273.15));
                                % 'vapdens21' is the density of the water vapour [kg/m^3] by using P,q (specific water) and T near the surface
                                vapdens2=((q/(1-q))*pressobs*100/(Rd*(temp+273.15)))*(1+(q/(1-q))*Rv/Rd);

                                % read wind data:
                                ws=str2num(linenext(30:33));
                                wdir=str2num(linenext(26:28));
                                if(typen==14 || typen==53 || typen==54)
                                    ws_fact=100;
                                elseif(typen==13)
                                    ws_fact=1;
                                end
                                ws=ws/ws_fact;
                                %---------------

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % read matrix Mread
                                Mread{sta_place,hour_place}=[Mread{sta_place,hour_place};pressobs,height,temp,rh,ws,vapdens1,vapdens2,...
                                    q,r,gammaw,latsta,lonsta,yearn,monthn,dayn,hourn,tdew,wdir];
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            iend=1;
                        end
                    end
                end
            else
                % if the comming data is not suitable for the intput station and date, just run over it quickly
                linenext=fgets(FID);
                lineend=str2num(linenext(1:2));
                while (lineend~=25) %flag 25 means line withoud data (3 lines afterwards will be header, or file ended)
                    linenext=fgets(FID);
                    lineend=str2num(linenext(1:2));
                end
            end
        else
            linestop=1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % At this stage matrix M contains all the levels data for the input station and time, before interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The next stage contatins 2 loops for the stations and hours in order to calcuate the sounding fileds if they exist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ww=0; % counter variable for test number of missing stations and time
    for st_num=1:length(station_list)
        for h_num=1:8
            clear M;
            M=Mread{st_num,h_num};
            nm=size(M,1);               % number of lines in M
            if (nm>=ml)                 % if the number of levels are equal or higher than "ml" (actually if the profile was really found in the file), calculate CAPE, CIN and TCWC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %display(['Found grid point close to station ' num2str(st_num) ' and hour ' num2str((h_num-1)*3) 'UTC no. of lines: ' num2str(nm)]);
%display(['Found grid point close to station ' num2str(st_num) ' and hour ' num2str((h_num-1)*3) 'UTC no. of lines: ' num2str(nm)]);


%%%            savename=['Mint2_data/Mint2_sim_',simtype,'_at_station_',num2str(st_num),'_on_',num2str(year),'_',num2str(month),'_',num2str(day),'_',num2str(h_num*3-3),'UTC.txt'];   % file savename
%%%            save('-ascii',savename,'M') % save the data in directory Mint_data
                pmin=50; % min pressure for interpolation
                pmax=1000; %max pressure for interpolation
                dp=50;% max pressure for interpolatio
               extrafields=get_press_levs(M,pmin,pmax,dp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                H11=ceil(M(1,2));       % round station height real number (to up)
                H22=floor(max(M(:,2))); % round upper model level height (to down)
                %height_step=10;        % The interpolation height step in meter unit. CHECK how bad the result becomes if we enlarge this value
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                DH=H11:height_step:H22;
                DH=DH';
                %display(['calculating CAPE for station ',num2str(st_num),' and hour ',num2str((h_num-1)*3),'UTC'])
                %display('INTREPOLATION STEP');
                pressint=interp1(M(:,2),M(:,1),DH,'pchip');     %pressure interpolation
                tempint=interp1(M(:,2),M(:,3),DH,'pchip');      %temperature interpolation
                tempdint=interp1(M(:,2),M(:,17),DH,'pchip');    %dew point interpolation
                vapdens1int=interp1(M(:,2),M(:,6),DH,'pchip');  %water vapor density by the method of using Tdew and Air temp
                vapdens2int=interp1(M(:,2),M(:,7),DH,'pchip');  %water vapor density by the method of using q, Air pressure and Air temp
                %wsint=interp1(M(:,2),M(:,5),DH,'pchip');        %wind speed interpolation
                %wdirint=interp1(M(:,2),M(:,18),DH,'pchip');     %wind direction interpolation
                Lint=2500800-2360.*tempint+1.6.*(tempint.^2)-0.06*(tempint.^3);   %latent heat released in condensation of kg (J/kg)
                eint=611.3.*exp((Lint./Rv).*(1/273.15-1./(tempdint+273.15)));	%water vapor pressure
                esint=611.3.*exp((Lint./Rv).*(1/273.15-1./(tempint+273.15)));	%water vapor pressure in saturation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Mwind=[];
                for iwind=1:size(M,1) % build Mwind matrix which is clear from wind vecloity anavailable data
                    if((M(iwind,5)>=0 && M(iwind,5)<=100) && (M(iwind,18)>=0 && M(iwind,18)<=360))
                        Mwind=[Mwind;M(iwind,2),M(iwind,5),-M(iwind,5)*sin(M(iwind,18)*pi/180),-M(iwind,5)*cos(M(iwind,18)*pi/180)];   % the matrix fileds are: height,total horizantal wind, u and v
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
                
                rhint=100*eint./esint;                          %rel.hum. interpolation
                rsatu=(epsi.*esint)./(pressint.*100-esint);       %spec. humidity if the rel.hum. was 100%
                qint=(eint.*epsi)./(pressint.*100-eint.*(1-epsi)); %spec. humidity interpolation = kg(vapor)/kg(dry air + vapor)
                rint=qint./(1-qint);          %mixing ratio interpolation = kg(vapor)/kg(dry air)
                gammawint=g.*(1+(Lint.*rsatu)./(Rd.*(tempint+273.15)))./(Cpd+(Lint.^2).*rsatu.*epsi./(Rd.*(tempint+273.15).^2)); %wet adiabat lapse rate
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % matrix Mint contains only the relvant station and relvant time with all levels data, after interpolation
                Mint=[pressint,DH,tempint,rhint,wsint,vapdens1int,vapdens2int,qint,rint,gammawint,Lint,eint,esint,uint,vint];
                
                if Mint(1,4)>=100
                   Mint(1,4)=99.9;
                end
                
                nline1=size(Mint,1);    %number of lines
                tempgrad=[];
                if height_step<0.00001
                    display('Too low height step. Exiting...');
                    return;
                end
                for ir=1:nline1-1;
                    ir1=ir;
                    ir2=ir+1;
                    tempgrad=[tempgrad;(Mint(ir1,3)-Mint(ir2,3))./height_step];  %temperature gradient
                end
                tempgrad(nline1)=tempgrad(nline1-1);
                Mint=[Mint,tempgrad]; % adding temperature gradient to Mint !
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                temp1=Mint(1,3);
                pressobs1=Mint(1,1);
                staheight=Mint(1,2);
                tempmass=temp1; % tempmass is the air parcel temperature [C]
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
                rfirst=(eint(1)*epsi)/(Mint(1,1)*100-eint(1));  % mixing ratio is constant during dry lapse
                erun=eint(1);                                   % vapor pressure
                esrun=611.3*exp((L(1)/Rv)*(1/273.15-1/(tempmass+273.15))); % saturation vapor pressure
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
                        % cape=g*dh*(Tvir_parcel-Tvir_sour)/Tvir_sour
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
                    cinadd=0;
                    if(Mcin(i1,2)<=10000 && Mcin(i2,1)<=10000) % Cin is caculated only for Heights below 10000m (10km)
                        cinadd=g*(Mcin(i2,2)-Mcin(i1,2))*((Mcin(i1,5)+Mcin(i2,5))/2-(Mcin(i1,6)+Mcin(i2,6))/2)/((Mcin(i1,6)+Mcin(i2,6))/2+273.15);   
                    end
                    if(cinadd<=0)           % add to the cin only the negative values (ignore superadiabatic region)
                        cin=cin+cinadd;
                    end
                end
                cin=-cin;
                
               %%%%%%%%%%%%%%%%%%%%%
                % Wind Shear Calculation :
                %%%%%%%%%%%%%%%%%%%%% 
                wshearout=zeros(length(windshear)/2,1)-999.9;% realoction of 3 windshears (The definitation of the vector
                %if strcmp(windshearopt,'scalar')                % the windshear cacluation by the scalar options
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
                %end
                %elseif strcmp(windshearopt,'vector')            % VECTORIC CACLUATION OF WIND SHEAR
                %    for jjj=1:length(windshear)/2               % the loop for 3 windshears 
                %        jjj1=jjj*2;                             % for the bottom level
                %        jjj2=jjj*2-1;                           % for the upper level
                %        if (windshear(jjj1)~=1100) 
                %            [mdiff ku]=min(abs(Mint(:,1)-windshear(jjj2))); % find the upper level 
                %            [mdiff kd]=min(abs(Mint(:,1)-windshear(jjj1))); % find the bottom level
                %        else
                %            [mdiff ku]=min(abs(Mint(:,1)-windshear(jjj2)));
                %            kd=1; % the lowest level of Mint matrix (shoud be the surface level or at least the closest to the surface level if surface level does not exist)
                %        end
                %        if ((length(ku)==1 && length(kd)==1) && abs((Mint(ku,1)-windshear(jjj2))/windshear(jjj2))<0.01 ...
                %        	&& (windshear(jjj1)==1100 || abs((Mint(kd,1)-windshear(jjj1))/windshear(jjj1))<0.01))   % if the pressure level difference between the closest interpolated level to the relevnt level is above 1% than don't calcute the windshear
                %   % if the bottom windshear==1100 than there is no problem
                %   % because it is the lowest or the surface level
                %   % otherwise the difference should be less than 1% , the
                %   % same condition that was done for the upper level is for the
                %   % bottom level.
                %            if ((Mint(ku,2)-Mint(kd,2))>0 && abs(Mint(ku,14)<100) && abs(Mint(kd,15)<100)) % no wind above 100m/s it is impossible and also the difference between the upper height to the lower height should be above zero
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wind shear scalaric calculation
                %                wshearout(jjj)=(sqrt((Mint(ku,14)-Mint(kd,14))^2+(Mint(ku,15)-Mint(kd,15))^2))/(Mint(ku,2)-Mint(kd,2)); % the windshear output
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %            end
                %        end
                %    end
                %end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % END of CAPE , CIN & WindShear calculation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                
                %%%%%%%%% Sanity check plot :
                %{
                if ~isempty(Mparcel) && ~isempty(Mdry) && ~isempty(Mwet) && ~isempty(Mcape)
                    close all;
                    %plot(Mint(:,3),Mint(:,2),'b')
                    %plot(Mint(:,3),Mint(:,2)/1000,'b')
                    %plot(Mint(:,3),Mint(:,2)/1000,'b'); hold on;
                    plot(Mparcel(:,6),Mparcel(:,2)/1000,'k'); hold on;
                    plot(Mparcel(:,5),Mparcel(:,2)/1000,'b'); hold on;
                    plot(Mdry(:,5),Mdry(:,2)/1000,'r'); hold on;
                    plot(Mwet(:,5),Mwet(:,2)/1000,'g'); hold on;
                    plot(Mcape(:,5),Mcape(:,2)/1000,'m'); hold on;
                    grid on;
                end
                %}
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
 %{                              
                for ii=1:length(varsound) % the number of output fileds
                    if(strcmp('CAPE',varsound(ii)))
                        output(ii,h_num,st_num)=cape;
                    elseif(strcmp('CIN',varsound(ii)))
                        output(ii,h_num,st_num)=cin;
                    elseif(strcmp('TCWC',varsound(ii)))
                        output(ii,h_num,st_num)=TCWC1;
                    elseif(strcmp('WSHEAR1',varsound(ii)))
                        output(ii,h_num,st_num)=wshearout(1);
                    elseif(strcmp('WSHEAR2',varsound(ii)))
                       output(ii,h_num,st_num)=wshearout(2);
                    elseif(strcmp('WSHEAR3',varsound(ii))) 
                       output(ii,h_num,st_num)=wshearout(3);
                    end
                end
%}
            output(1:18,h_num,st_num)=[cape,cin,TCWC1,wshearout(1:3)',extrafields];
            else
                output(:,h_num,st_num)=-999.9;
            end
        end % end of the hour loop ('h_num' variable)
    end % end of stations place loop ('st_num' variable)
else
    display(['The COSMO vertical FILE ' fullfname '  does not exist']);
end
