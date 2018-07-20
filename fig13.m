% make multiple posterior plots
distanceShow=0;

FS=6; % Fontsize
FSL=14; % Fontsize labels
FSL2=14; % tickmarks
markerObs='b'; % Label
markerMW='r';


% Posterior type
% 0 = none
% 1 = square cut-off in log-space
% 2 = spatial Poisson
% 3 = individual well-mixed
% 4 = simple colonization model
% 5 = advanced colonization model
% 6 = dark biosphere
% 7 = independent biospheres
% 8 = prehistoric intelligence
% 9 = extant aliens
% 10 = colonization, no extinction during expansion
% 11 = No K3 supercivs in 10^5 galaxies
% 12 = no past MW civ
plist=[1 2 3 4 5 11];

clf

for pindex=1:length(plist)
    posterior = plist(pindex);
    
    generatePosterior; % script selecting the consistent cases
    makeHistograms;
    

    subplot(2,1,1)
%    fill(10.^[x(1) x x(end)],[1e-7 HNposterior 1e-7],[0.7 1 0.7],'FaceAlpha',0.5)
if (pindex==1) semilogx(1,1,'w'); end

    semilogx(10.^x,HNposterior)
    hold on
    drawnow
    
    
    subplot(2,1,2)
    if (pindex==1) semilogx(1,1,'w'); end
    y=100*cumsum(HNposterior);
    y=y/y(end);
%plot(x,y)
semilogx(10.^x,y)
hold on
drawnow 

end


% lines
YTop=0.03;
subplot(2,1,1)
Nobs=log10(1/150e9); Ngal=log10(1);
l1=plot(10.^[Nobs Nobs],[0 YTop],markerObs);
l2=plot(10.^[Ngal Ngal],[0 YTop],markerMW)
legend([l1,l2],'Alone in observable universe','Alone in Milky Way','Location','northoutside','Orientation','horizontal')

subplot(2,1,2)
l1=plot(10.^[Nobs Nobs],[0 1],markerObs);
l2=plot(10.^[Ngal Ngal],[0 1],markerMW);

grid on
ylabel('P(N<x) (%)','FontSize',FSL)

set(gca,'FontSize',FSL2);
xlabel('N','FontSize',FSL)

axis([1e-40 1e15 0 1])
subplot(2,1,1)
ylabel('Frequency','FontSize',FSL)
axis([1e-40 1e15 0 0.04])



if (distanceShow>0)
    subplot(2+distanceShow,1,3)
    
    d=0:0.1:15; % distance in log10 pc
    
    %             V=1000*pi*50000^2; % ly
    V=2.3e11; % pc^3
    MWStars=3e11;
    sigfunc = @(A, x)(A(4)+A(1)./(1+exp(A(2)*(x-A(3)))));
    A=[6.5915    1.5716   21.9808  -26.6385];
    rho=(3.08567758e16^3/2e30)*10.^(sigfunc(A,d+16.4894)); % add to calculate in parsec
    
    subsample=.01; % fraction of runs to use
    ffc=find(consistent); % but if there are fewer ones for the posterior, use that instead
    runs=min(length(ffc), round(N*subsample))
    dcdf=zeros(runs,length(d));
    if (plotType==12) dcdfP=zeros(runs,length(d));  end
    j=1;
    for i=1:round(subsample*N)
        if (rem(i,10000)==0) disp(i); end
        
        dcdf(i,:)=1-exp(-(4*pi/3)*rho.*(10.^(3*d + log10N(i) - log10(MWStars))));
        %             dcdf(i,:)=1-exp(-(4*pi/3)*(10.^(log10N(i)-log10(V)+3*d)));
        
        if (plotType==12)   % Make plot for posterior
            dcdfP(i,:)=1-exp(-(4*pi/3)*rho.*(10.^(3*d + log10N(ffc(j)) - log10(MWStars))));
            j=j+1;
        end
        if (j>length(ffc)) break; end
    end
    
    %% convert to pc!!!
    semilogx(10.^(d),100*mean(dcdf))
    hold on
    if (plotType==12) semilogx(10.^(d),100*mean(dcdfP)); end
    
    plot(10.^([5 5]),[0 100],markerMW)
    hh=100/10; %  text(10^(5-.13),hh,'Milky way','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
    
    plot(10.^([1 1]*log10(93e9)),[0 100],markerObs)
    hh=100/10;  % text(10.^(log10(93e9)-.13),hh,'Observable universe','FontSize',FS,'Rotation',90,'VerticalAlignment','bottom')
    
    %axis([-40 15 0 1])
    axis([10^0 10^15 0 100])
    
    set(gca,'YTick',0:25:100)
    
    set(gca,'FontSize',FSL2);
    
    xlabel('d (pc)','FontSize',FSL)
    ylabel('P(D<d) (%)','FontSize',FSL)
    grid on
    
    
    %ax=gca;
    %ax.XTick = [0:15];
end
drawnow
