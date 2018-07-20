distanceShow=1;

FS=6; % Fontsize
FSL=14; % Fontsize labels
FSL2=14; % tickmarks
markerObs='b'; % Label
markerMW='r';

%% Make posteriors
HNposteriorStack=[];

maxposteriors=5;
leg={};
saveConsistent=zeros(N,maxposteriors);

for pp=1:maxposteriors
    colonyTime=40e6;
    %    colonyTime=1;
    switch pp
        case 1
            posterior=0; % Base case
            leg=[leg; 'No update'];
        case 2
            posterior=3;
            leg=[leg; 'Random sampling'];
        case 3
            posterior=2;
            leg=[leg; 'Spatial Poisson'];
        case 4
            % G-hat survey K3 civs
            PK3=0.5; % Probability that becoming K3 is a good thing
            PK3success = 0.01; % Probability that achieves K3
            KK3=1e5; % number of galaxies not seen as having K3
            posterior=11;
            leg=[leg; 'No K3 civilization observed'];
        case 5
            posterior=13; % No past colonization
            leg=[leg; 'Settlement update'];
    end
    generatePosterior; % script selecting the consistent cases
    Nconsistent=sum(consistent>0);
    saveConsistent(:,pp)=(consistent>0);
    makeHistograms;
    HNposteriorStack=[HNposteriorStack; HNposterior];
    
    
    %% Make table
    Nconsistent=sum(consistent>0);
    ffc=find(consistent);
    meanN = mean(10.^log10N(ffc));
    medianN = median(10.^log10N(ffc));
    PrN1 = sum(log10N(ffc)<0)/Nconsistent;
    PrN10 = sum(log10N(ffc)<log10(1/150e9))/Nconsistent;
    meanFl = mean(10.^fl(ffc));
    medianFl = median(10.^fl(ffc));
    meanL = mean(10.^L(ffc));
    medianL = median(10.^L(ffc));
    
    logshiftMeanFl = (mean(fl(ffc))-mean(fl));
    logshiftMedianFl = (median(fl(ffc))-median(fl));
    logshiftMeanL = (mean(L(ffc))-mean(L));
    logshiftMedianL =(median(L(ffc))-median(L));
    
    PaloneGalaxy =  mean(exp(-10.^log10N));
    
    %  ind=find(100*mean(dcdfP)>=50,1);
    %   if (isempty(ind))
    %       medD=NaN;
    %   else
    %       medD=10^d(ind); % median distance in pc to nearest
    %   end
    medD=NaN; % placeholder!
    
  %  disp( sprintf('%s&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}',char(leg(pp)),num2str(meanN,2),num2str(medianN,2),num2str(PrN1,2),num2str(PrN10,2),num2str(medD,2),num2str(meanFl,2),num2str(medianFl,2),num2str(meanL,2),num2str(medianL,2),num2str(logshiftMeanFl,2),num2str(logshiftMedianFl,2),num2str(logshiftMeanL,2),num2str(logshiftMedianL,2)))
  
    
 %   disp( sprintf('%s&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}',char(leg(pp)),num2str(meanN,2),num2str(medianN,2),num2str(PrN1,2),num2str(PrN10,2),num2str(medianFl,2),num2str(medianL,2),num2str(logshiftMedianFl,2),num2str(logshiftMedianL,2)))
    disp( sprintf('%s&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}&\\num{%s}',char(leg(pp)),num2str(meanN,2),num2str(medianN,2),num2str(PrN1,2),num2str(PrN10,2),num2str(medianFl,2),num2str(medianL,2)))
  
    
    
    %disp(sum(L(consistent)>3)/N), % Probability of L>1000
end


%% PDF
clf

subplot(3,1,1)
%semilogx(1,1,'w');
%fill(10.^[x(1) x x(end)],[1e-6 HN 1e-6],[0.7 0.7 1])
%semilogx(10.^x,HN)
for i=1:maxposteriors
    %fill(10.^[x(1) x x(end)],[1e-6 HNposteriorStack(i,:) 1e-6],[0.7 1 0.7],'FaceAlpha',0.5)
    semilogx(10.^x,HNposteriorStack(i,:))
    hold on
end
grid on

height=max(HNposteriorStack(:))*1.05;

% % lines
Nobs=log10(1/150e9); Ngal=log10(1);
l1=plot(10.^[Nobs Nobs],[0 height],markerObs);
% hh=max(HN)/10; %text(10.^(Nobs-.3),hh,'Alone in observable universe','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
l2=plot(10.^[Ngal Ngal],[0 height],markerMW);
% %text(10.^(Ngal-.3),hh,'Alone in Milky Way','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
%


set(gca,'FontSize',FSL2);

axis([10^-40 10^15 0 height])
%xlabel('log_{10}(N)')
xlabel('N','FontSize',FSL)
ylabel('Frequency','FontSize',FSL)

legend(leg,'Location','west')


% %% CDF
subplot(3,1,2)

for i=1:maxposteriors
    semilogx(10.^x,100*cumsum(HNposteriorStack(i,:)))
    hold on
end
grid on

% % lines
Nobs=log10(1/150e9); Ngal=log10(1);
plot(10.^[Nobs Nobs],[0 100],markerObs);
%iNobs=find(x>Nobs,1);   plot(10^Nobs, interp1([x(iNobs-1) x(iNobs)],[y(iNobs-1) y(iNobs)],Nobs),'o')
hh=100/10; %  text(10^(Nobs-.3),hh,'Alone in observable universe','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
plot(10.^[Ngal Ngal],[0 100],markerMW)
%iNgal=find(x>Ngal,1);   plot(10^Ngal, interp1([x(iNgal-1) x(iNgal)],[y(iNgal-1) y(iNgal)],Ngal),'o')
hh=100/10; %   text(10^(Ngal-.3),hh,'Alone in Milky Way','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')

axis([10^-40 10^15 0 100])
set(gca,'YTick',0:25:100)
set(gca,'FontSize',FSL2);

xlabel('N','FontSize',FSL)
ylabel('P(N<x) (%)','FontSize',FSL)



%% Distances
subplot(3,1,3)

d=0:0.1:15; % distance in log10 pc
V=2.3e11; % pc^3
MWStars=3e11;
sigfunc = @(A, x)(A(4)+A(1)./(1+exp(A(2)*(x-A(3)))));
A=[6.5915    1.5716   21.9808  -26.6385];
rho=(3.08567758e16^3/2e30)*10.^(sigfunc(A,d+16.4894)); % add to calculate in parsec

subsample=.001; % fraction of runs to use

for ii=1:maxposteriors
    
    ffc=find(saveConsistent(:,ii)); % but if there are fewer ones for the posterior, use that instead
    runs=min(length(ffc), round(N*subsample))
    dcdf=zeros(runs,length(d));
    dcdfP=zeros(runs,length(d));
    j=1;
    for i=1:round(subsample*N)
      %  if (rem(i,10000)==0) disp(i); end
        dcdf(i,:)=1-exp(-(4*pi/3)*rho.*(10.^(3*d + log10N(i) - log10(MWStars))));
        dcdfP(i,:)=1-exp(-(4*pi/3)*rho.*(10.^(3*d + log10N(ffc(j)) - log10(MWStars))));
        j=j+1;
        if (j>length(ffc)) break; end
    end
    
    semilogx(10.^(d),100*mean(dcdfP));
    
    q=100*mean(dcdfP);
    i2=find(q>=50,1);
    if (~isempty(i2))
 loc=interp1([q(i2-1) q(i2)],10.^[d(i2-1) d(i2)],50);
    else
        loc=NaN;
    end
    disp(loc)
 
    
hold on
drawnow
end

l2=plot(10.^([5 5]),[0 100],markerMW);
hh=100/10; %  text(10^(5-.13),hh,'Milky way','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
l1=plot(10.^([1 1]*log10(93e9)),[0 100],markerObs);
hh=100/10;  % text(10.^(log10(93e9)-.13),hh,'Observable universe','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
axis([10^0 10^15 0 100])
set(gca,'YTick',0:25:100)
set(gca,'FontSize',FSL2);

xlabel('d (pc)','FontSize',FSL)
ylabel('P(D<d) (%)','FontSize',FSL)
grid on


legend([l1,l2],'Alone in observable universe','Alone in Milky Way','Location','southoutside','Orientation','horizontal')
