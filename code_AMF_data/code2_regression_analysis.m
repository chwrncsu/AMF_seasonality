clear all
close all
clc

load output1.mat results
vari=[11 13 14 15 18];
% 1: SlNo; 2: ID; 3: WRR; 4: Region
% 5: LAT; 6: LON
% 7: SI1; 8: SI2; 9: Theta1; 10: Theta2
% 11: SI2-SI1; 12: Theta2-Theta1 13: Precip2-Precip1; 14: Temp2-Temp1;
% 15: DA; 16: Slope1; 17: Length; 18: Elev1;

regs=unique(results(:,4));
R=length(regs);

clear a; a=find(abs(results(:,11))<0.15); dummySI0=results(a,:);
clear dummy
for s=1:size(dummySI0,1)
    if(abs(dummySI0(s,12))<=30)
        dummy(s,1)=0;
    else
        dummy(s,1)=1;
    end
end
clear a; a=find(dummy==0); SI_Theta_0_0=dummySI0(a,:);
clear a; a=SI_Theta_0_0(:,1); SI_Theta_other=results; SI_Theta_other(a,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SI_change(11) vs Precip_change(13); Temp_change(14) %
for r=1:R
    clear dummy dummya dummyb
    clear a; a=find(SI_Theta_other(:,4)==regs(r));
    dummy=SI_Theta_other(a,:);
    clear a; a=find(dummy(:,11)<0); dummya=dummy(a,:);
    clear a; a=find(dummy(:,11)>=0); dummyb=dummy(a,:);
    Regions(r,1)={dummya}; Regions(r,2)={dummyb};
end
for r=1:R
    for i=1:2
        clear lm; lm={};
        clear dummy; dummy=[];
        if(size(Regions{r,i},1)>=5)
            dummy=Regions{r,i}(:,vari);
            dummy(:,4:5)=log10(dummy(:,4:5));
            clear dummyx dummyy dummycomb K V ;
            dummyy=dummy(:,1); dummyx=dummy(:,2:end); K=size(dummyx,2); V=1:K;
            count=0;
            for j=1:K
                clear dummy1; dummy1=nchoosek(V,j);
                for k=1:size(dummy1,1)
                    count=count+1;
                    dummycomb(count,1)={dummy1(k,:)};
                end
            end
            for j=1:length(dummycomb)
                clear mdl
                mdl=fitlm(dummyx(:,dummycomb{j,1}),dummyy,'linear');
                clear aic no k SSE; no=mdl.NumObservations; k=mdl.NumVariables; SSE=mdl.SSE;
                aic=no*log(SSE/no)+2*k;
                lm{j,1}=mdl;
                lm{j,2}=aic;
            end
            lm(:,3)=dummycomb;
        end
        LM(r,i)={lm};
    end
end

clear SLOPE; SLOPE={};
for r=1:R
    for i=1:2
        clear dummy aic slope
        dummy=LM{r,i};
        slope=ones(2,K)*NaN;
        if(length(dummy)>0)
            for j=1:size(dummy,1)
                aic(j,1)=dummy{j,2};
            end
            AIC(r,i)={aic};
            clear a; a=find(aic==min(aic)); BB(r,i)=a; AICmin(r,i)=aic(a,1);
            clear dummy1 dummy2
            dummy1=table2array(LM{r,i}{a,1}.Coefficients); dummy2=LM{r,i}{a,3};
            for j=1:length(dummy2)
                slope(1,dummy2(j))=dummy1(j+1,1);
                slope(2,dummy2(j))=dummy1(j+1,4);
            end
        end
        SLOPE(r,i)={slope};
    end
end
open SLOPE

for r=1:R
    for i=1:2
        clear dummy a1 a2 p1 p2 p3 p4; dummy=[];
        if(size(Regions{r,i},1)>=5)
            dummy=Regions{r,i}(:,[vari 12]);
            dummy(:,4:5)=log10(dummy(:,4:5));
        end
        for j=1:size(dummy,1)
            if(abs(dummy(j,length(vari)+1))<30)
                dummy(j,length(vari)+2)=0;
            else
                dummy(j,length(vari)+2)=1;
            end
        end
        if(length(dummy)==0)
        else
            p1=polyfit(dummy(:,2),dummy(:,1),1);
            p2=polyfit(dummy(:,3),dummy(:,1),1);
            p3=polyfit(dummy(:,4),dummy(:,1),1);
            p4=polyfit(dummy(:,5),dummy(:,1),1);
            %dummy=sortrows(dummy,length(vari)+2);
            clear b; b=find(dummy(:,end)==0); a1=dummy(b,:);
            clear b; b=find(dummy(:,end)==1); a2=dummy(b,:);
            data(r,i)={dummy};
            clear a; a=r*100+i*10;
            figure(a+1)
            hold on
            plot(a1(:,2),a1(:,1),'o','color',[0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5])
            plot(a2(:,2),a2(:,1),'o','color',[1 0 0],'MarkerFaceColor', [1 0 0])
            clear bb1 bb2; bb1=axis; bb2=linspace(bb1(1),bb1(2),50);
            plot(bb2,polyval(p1,bb2),'b--')
            xlabel('Change in Precipitation')
            ylabel('Change in SI')
            set(gca,'fontsize',15)
            figure(a+2)
            hold on
            plot(a1(:,3),a1(:,1),'o','color',[0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5])
            plot(a2(:,3),a2(:,1),'o','color',[1 0 0],'MarkerFaceColor', [1 0 0])
            clear bb1 bb2; bb1=axis; bb2=linspace(bb1(1),bb1(2),50);
            plot(bb2,polyval(p2,bb2),'b--')
            xlabel('Change in Temperature')
            ylabel('Change in SI')
            set(gca,'fontsize',15)
            figure(a+3)
            hold on
            plot(a1(:,4),a1(:,1),'o','color',[0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5])
            plot(a2(:,4),a2(:,1),'o','color',[1 0 0],'MarkerFaceColor', [1 0 0])
            clear bb1 bb2; bb1=axis; bb2=linspace(bb1(1),bb1(2),50);
            plot(bb2,polyval(p3,bb2),'b--')
            xlabel('log_1_0(DA)')
            ylabel('Change in SI')
            set(gca,'fontsize',15)
            figure(a+4)
            hold on
            plot(a1(:,5),a1(:,1),'o','color',[0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5])
            plot(a2(:,5),a2(:,1),'o','color',[1 0 0],'MarkerFaceColor', [1 0 0])
            clear bb1 bb2; bb1=axis; bb2=linspace(bb1(1),bb1(2),50);
            plot(bb2,polyval(p4,bb2),'b--')
            xlabel('log_1_0(Elev)')
            ylabel('Change in SI')
            set(gca,'fontsize',15)
        end
    end
end
% open data