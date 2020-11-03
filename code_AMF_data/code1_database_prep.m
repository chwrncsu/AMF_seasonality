clear all
close all
clc

% This code prepares the data for regression analysis.
% The code calls findmonth.m to estimate what date belongs to which month.

Tyear=1970;

% The INPUT file contains:
% 1. stream: Annual Maximum Flow (AMF) data for each of the 975 stations in
% every cell. Each cell has 5 columns, the 1st column is the USGS station
% ID, 2nd column is AMF value, 3rd to 5th column are the month, day and
% year when the AMF was recorded
% 2. latlon: Latitude and Longitude of each station
% 3. wrr: Water Resources Region number of each station
% 4. attr: Drainage area(DA), Slope, Length, Elevation of each of the 975
% stations obtained from Bureau of Reclamation
% 5: PRECIPdaily: Daily maximum precipitation for a given month
% corresponding to each year. Each cell represent precipitation data for
% each station. Inside each cell, the number of rows corresponds to every
% year's data and each column (1-12) represent every month, with column 1
% is January and column 12 is December. The 13th column is the year.
% 6: MAXtemp: Daily maximum temperature for a chosen month. The format of
% temperature data is the same as the precipitation data.

load output0.mat %stream latlon wrr attr PRECIPdaily MAXtemp

N=size(stream,1);
for i=1:N
    clear dummy0 dummy1 dummydays
    dummy0=stream{i,1}; dummy1(:,[1 3])=dummy0(:,[2 5]);
    for j=1:size(dummy0,1)
        dummydays(j,1)=datenum(dummy0(j,5),12,31)-datenum(dummy0(j,5),1,1)+1;
    end
    mbar(i,1)=mean(dummydays);
    for j=1:size(dummy0,1)
        dummy1(j,4)=datenum(dummy0(j,[5 4 3]))-datenum(dummy0(j,5),1,1)+1;
        dummy1(j,2)=dummy1(j,4)*2*pi/mbar(i,1);
    end
    stream(i,2)={dummy1};
end


% 1: SlNo; 2: ID; 3: WRR; 4: Region
% 5: LAT; 6: LON
% 7: SI1; 8: SI2; 9: Theta1; 10: Theta2
% 11: SI2-SI1; 12: Theta2-Theta1 13: Precip2-Precip1; 14: Temp2-Temp1;
% 15: DA; 16: Slope1; 17: Length; 18: Elev1;

DATA(:,1)=1:N; % Sl No. %
for s=1:N
    DATA(s,2)=stream{s,1}(1,1); % ID %
end
DATA(:,3)=wrr; % WRR %
for s=1:N
    if(DATA(s,3)==1||DATA(s,3)==2) DATA(s,4)=1; % Region %
    elseif(DATA(s,3)==3||DATA(s,3)==6) DATA(s,4)=2; % Region %
    elseif(DATA(s,3)==4||DATA(s,3)==5) DATA(s,4)=3; % Region %
    elseif(DATA(s,3)==7||DATA(s,3)==9||DATA(s,3)==10) DATA(s,4)=4; % Region %
    elseif(DATA(s,3)==8||DATA(s,3)==11||DATA(s,3)==12) DATA(s,4)=5; % Region %
    elseif(DATA(s,3)==13||DATA(s,3)==14||DATA(s,3)==15||DATA(s,3)==16) DATA(s,4)=6; % Region %
    elseif(DATA(s,3)==17) DATA(s,4)=7; % Region %
    elseif(DATA(s,3)==18) DATA(s,4)=8; % Region %
    end
end
% WRR1+WRR2=Region(1);
% WRR3+WRR6=Region(2);
% WRR4+WRR5=Region(3);
% WRR7+WRR9+WRR10=Region(4);
% WRR8+WRR11+WRR12=Region(5);
% WRR13+WRR14+WRR15+WRR16=Region(6);
% WRR17=Region(7);
% WRR18=Region(8);
DATA(:,5)=latlon(:,1); % LAT %
DATA(:,6)=-latlon(:,2); % LON %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:N
    clear dummy dummy1 dummy2
    dummy=stream{s,2};
    n=size(dummy,1);
    dummy1=dummy(:,3);
    dummy2=find(dummy1<=Tyear);
    record(s,1)=length(dummy2); record(s,2)=n-length(dummy2);
end
fprintf('record before:%0.2f\nrecord after:%0.2f\n',mean(record));
clear datasort
for s=1:N
    if(record(s,1)>=1 && record(s,2)>=1)
        clear dummy dummy1 dummy2 dummy3
        dummy=stream{s,2};
        n=size(dummy,1);
        dummy1=dummy(:,3);
        dummy2=find(dummy1<=Tyear);
        dummy3=find(dummy1>Tyear);
        datasort{s,1}=dummy(dummy2,:);
        datasort{s,2}=dummy(dummy3,:);
    end
end
clear SI1
for s=1:N
    clear amf
    amf=datasort{s,1}(:,1:2);
    sumx=0; sumy=0;
    for i=1:size(amf,1)
        sumx=sumx+cos(amf(i,2));
        sumy=sumy+sin(amf(i,2));
    end
    xbar=sumx/size(amf,1); ybar=sumy/size(amf,1);
    r=sqrt(xbar^2+ybar^2);
    if(xbar<=0)
        theta=(atan(ybar/xbar)+pi)*mbar(s,1)/(2*pi);
    elseif(xbar>0 && ybar>=0)
        theta=atan(ybar/xbar)*mbar(s,1)/(2*pi);
    elseif(xbar>0 && ybar<0)
        theta=(atan(ybar/xbar)+2*pi)*mbar(s,1)/(2*pi);
    end
    SI1(s,1)=r; SI1(s,2)=theta; %SI1(s,3)=r/sum(amf(:,1));
end
clear SI2
for s=1:N
    clear amf
    amf=datasort{s,2}(:,1:2);
    sumx=0; sumy=0;
    for i=1:size(amf,1)
        sumx=sumx+cos(amf(i,2));
        sumy=sumy+sin(amf(i,2));
    end
    xbar=sumx/size(amf,1); ybar=sumy/size(amf,1);
    r=sqrt(xbar^2+ybar^2);
    if(xbar<=0)
        theta=(atan(ybar/xbar)+pi)*mbar(s,1)/(2*pi);
    elseif(xbar>0 && ybar>=0)
        theta=atan(ybar/xbar)*mbar(s,1)/(2*pi);
    elseif(xbar>0 && ybar<0)
        theta=(atan(ybar/xbar)+2*pi)*mbar(s,1)/(2*pi);
    end
    SI2(s,1)=r; SI2(s,2)=theta; %SI2(s,3)=r/sum(amf(:,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA(:,7)=SI1(:,1); % SI1 %
DATA(:,8)=SI2(:,1); % SI2 %
DATA(:,9)=SI1(:,2); % Theta1 %
DATA(:,10)=SI2(:,2); % Theta2 %
    DATA(:,11)=DATA(:,8)-DATA(:,7); % SI2-SI1 %
DATA(:,12)=DATA(:,10)-DATA(:,9); % Theta2-Theta1 %
for i=1:N
    if(abs(DATA(i,12))>365/2)
        if(DATA(i,12)<0) DATA(i,12)=DATA(i,12)+mbar(i,1);
        elseif(DATA(i,12)>0) DATA(i,12)=DATA(i,12)-mbar(i,1);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:N
    clear dummy dummy1 dummy2
    m1=findmonth(SI1(s,2)); m2=findmonth(SI2(s,2));
    months(s,1)=m1; months(s,2)=m2;
    dummy=PRECIPdaily{s};
    if(length(dummy>0))
        clear d; d=find(dummy(:,13)<=1970); dummy1=dummy(d,:);
        clear d; d=find(dummy(:,13)>=1971); dummy2=dummy(d,:);
        precip_temp(s,1)=mean(dummy1(:,m1));
        precip_temp(s,2)=mean(dummy2(:,m1));
        precip_temp(s,3)=mean(dummy1(:,m2));
        precip_temp(s,4)=mean(dummy2(:,m2));
    else
        precip_temp(s,1:4)=NaN;
    end
    clear dummy dummy1 dummy2
    dummy=MAXtemp{s};
    if(length(dummy>0))
        clear d; d=find(dummy(:,13)<=1970); dummy1=dummy(d,:);
        clear d; d=find(dummy(:,13)>=1971); dummy2=dummy(d,:);
        precip_temp(s,5)=mean(dummy1(:,m1));
        precip_temp(s,6)=mean(dummy2(:,m1));
        precip_temp(s,7)=mean(dummy1(:,m2));
        precip_temp(s,8)=mean(dummy2(:,m2));
    else
        precip_temp(s,5:8)=NaN;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA(:,13)=precip_temp(:,2)-precip_temp(:,1); % Percip2-Percip1 %
DATA(:,14)=precip_temp(:,6)-precip_temp(:,5); % Temp2-Temp1 %
DATA(:,15:18)=attr(:,1:4); % DA, Slope, Length, Elev %
clear results; results=DATA;
save output1.mat results Tyear