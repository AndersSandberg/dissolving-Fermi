%% Data source
% 1: Data from our guesstimates
% 2: Data from literature
% 3: Data with just life, no uncertainty other factors
dataSource=1;

% For lit data, splice in a life rate-probability model?
spliceInLife=0;

% For uniform distributions, use lin-uniform
linuniform=0;

% Rate estimate intelligence
rateintell=0;

% Remove the k highest or lowest literature claims
OutlierRemoval=0;

% Perform a bootstrap resample of the literature before resampling it
BootStrap=0;

%% Posterior type
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
posterior = 0;

% Parameters in case 1
% par1: log10(cutoffN), par2 remaining probl
par1=2; par2=0.1;

% Detection limit case 2
%dDetect=60; % Detection limit in lightyears
dDetect=18; % Detection limit in parsec
probfaildetect=0.0; % Probability of missing a civ

% Number of checked stars for case 3
NChecked=1e3;

% Time for colonizing the galaxy in posterior 5
colonyTime=40e6;
%colonyTime=50000; % fastest possible
%colonyTime=250e6; % one galactic year, panspermia like

% Beta parameters for external biospheres
nSearched=100;
nFound=0;

% G-hat survey K3 civs
PK3=0.5; % Probability that becoming K3 is a good thing
PK3success = 0.01; % Probability that achieves K3
KK3=1e5; % number of galaxies not seen as having K3
    

%% Plot type
% 1=PDF
% 2=CDF
% 3=rate;
% 4=distance to nearest civ
% 5=Scatterplot of fl vs rest.
% 6=Show factor distributions
% 10 = fig1 paper.
% 11 = fig2 paper
% 12 = fig showing posterior
% 13 = plot of L vs fl. Dependent on logRange below
% 14 = plot several posteriors 
% 15 = plot posteriro for paper
plotType=11;

% Show just the number
onlyN=0;

% Do not erase previous plot
noClf=1;

% Semilog plots for pdf, cdf
logPlot=0;

histBins=100;

% Range in logspace to plot
%logRange=[-100 .1 50];
logRange=[-100 .3 50];


showstats=0; % For plottype 6, show statistical fits

%% Data printout
printout=0;

%% Number of Monte Carlo realizations
N=1e7;

%% Reset random number generator

rng('default')

%% Create data
switch(dataSource)
    case 1
        % Astrophysics
        %Ns=log10(unifrnd(1*ones(N,1),100));
        Ns=unifrnd(0*ones(N,1),2);
        fp=unifrnd(-1*ones(N,1),0);
        ne=unifrnd(-1*ones(N,1),0);
        
        if (linuniform==1)
        Ns=log10(unifrnd(1*ones(N,1),100));
        fp=log10(unifrnd(0.1*ones(N,1),1));
        ne=log10(unifrnd(0.1*ones(N,1),1));
        end
        
        % Life
        %fl=unifrnd(-30*ones(N,1),0); % fraction estimate
        
        % From rate estimate
        % rate per planet up till now
        loglambdalife=normrnd(-0*ones(N,1),50);
        
        %loglambdalife=unifrnd(-50*ones(N,1),50); %%% Log-uniform
        %   loglambdalife=normrnd(-70*ones(N,1),50);
        
        fl=log10(1-exp(-10.^loglambdalife));
        fl(loglambdalife<-10)=loglambdalife(loglambdalife<-10); % Fix the underflow numbers
        
        % fl=unifrnd(-100*ones(N,1),0); $uniform fl distribution
        %if (linuniform==1)
        %    fl=log10(unifrnd(0*ones(N,1),1));
        %end
        
        % Intelligence
        fi=unifrnd(-3*ones(N,1),0); % fraction estimate
        % From rate estimate
        if (rateintell==1)
            loglambda=normrnd(-0*ones(N,1),50); % 0,50
            fi=log10(1-exp(-10.^loglambda));
            fi(loglambda<-10)=loglambda(loglambda<-10); % Fix the underflow numbers
        end
        if (linuniform==1)
            fi=log10(unifrnd(0.001*ones(N,1),1));
        end
        
        % Communication
        fc=unifrnd(-2*ones(N,1),0);
        if (linuniform==1)
            fc=log10(unifrnd(0.01*ones(N,1),1));
        end
        
        % Longevity
        L=unifrnd(2*ones(N,1),10);
        if (linuniform==1)
            L=log10(unifrnd(100*ones(N,1),1e10));
        end
        
    case 2  %% Data from literature
        litResample;
        
        if (spliceInLife)
            % Life
            % From rate estimate
            % rate per planet up till now
            loglambdalife=normrnd(-0*ones(N,1),50);
            fl=log10(1-exp(-10.^loglambdalife));
            fl(loglambdalife<-10)=loglambdalife(loglambdalife<-10); % Fix the underflow numbers
            
        end
        
        %% No uncertainty except fl, all values max optimistic
    case 3
        Ns=unifrnd(log10(50)*ones(N,1),log10(50));
        fp=unifrnd(0*ones(N,1),0);
        ne=unifrnd(1*ones(N,1),1);
        
        % Life
        % From rate estimate
        % rate per planet up till now
        loglambdalife=normrnd(-0*ones(N,1),50);
        
        fl=log10(1-exp(-10.^loglambdalife));
        fl(loglambdalife<-10)=loglambdalife(loglambdalife<-10); % Fix the underflow numbers
        
        % Intelligence
        fi=unifrnd(0*ones(N,1),0); % fraction estimate
        
        % Communication
        fc=unifrnd(0*ones(N,1),0);
        
        % Longevity
        L=unifrnd(9*ones(N,1),9);
end

%% Calc
log10N=(Ns+fp+ne+fl+fi+fc+L);

%% Posterior

generatePosterior; % script selecting the consistent cases

Nconsistent=sum(consistent>0);

%% Plot
makeHistograms;

% Inverse cube root of density, for distance estimation in ly
%lamb=(10.^log10N/(pi*1000*50000^2)).^(-1/3);


%% Actual plotting

switch plotType
    case 1
        %% PDF plots
        if (noClf==0) clf; end
        if (logPlot==1)
            if (onlyN)
                semilogy(x,HN,x,HNposterior)
            else
                subplot(8,1,1);
                semilogy(x,HN,x,HNposterior)
                subplot(8,1,2)
                semilogy(x,HNs,x,HNsposterior)
                subplot(8,1,3)
                semilogy(x,Hfp,x,Hfpposterior)
                subplot(8,1,4)
                semilogy(x,Hne,x,Hneposterior)
                subplot(8,1,5)
                semilogy(x,Hfl,x,Hflposterior)
                subplot(8,1,6)
                semilogy(x,Hfi,x,Hfiposterior)
                subplot(8,1,7)
                semilogy(x,Hfc,x,Hfcposterior)
                subplot(8,1,8)
                semilogy(x,HL,x,HLposterior)
            end
        else
            if (onlyN)
                plot(x,HN,x,HNposterior)
            else
                subplot(8,1,1);
                plot(x,HN,x,HNposterior)
                subplot(8,1,2)
                plot(x,HNs,x,HNsposterior)
                subplot(8,1,3)
                plot(x,Hfp,x,Hfpposterior)
                subplot(8,1,4)
                plot(x,Hne,x,Hneposterior)
                subplot(8,1,5)
                plot(x,Hfl,x,Hflposterior)
                subplot(8,1,6)
                plot(x,Hfi,x,Hfiposterior)
                subplot(8,1,7)
                plot(x,Hfc,x,Hfcposterior)
                subplot(8,1,8)
                plot(x,HL,x,HLposterior)
            end
        end
        set(gca,'XTick',-100:10:50)
    case 3
        %% Rate plots
        if (noClf==0) clf; end
        [Hllife,rate]=hist(loglambdalife,1000);
        [Hllifeposterior,rate]=hist(loglambdalife(consistent),rate);
        plot(rate,Hllife/N,rate,Hllifeposterior/Nconsistent)
        
    case 2
        %% CDF graphs
        if (noClf==0) clf; end
        if (logPlot==1)
            if (onlyN)
                semilogy(x,cumsum(HN),x,cumsum(HNposterior))
                grid on
            else
                subplot(8,1,1)
                semilogy(x,cumsum(HN),x,cumsum(HNposterior))
                grid on
                subplot(8,1,2)
                semilogy(x,cumsum(HNs),x,cumsum(HNsposterior))
                grid on
                subplot(8,1,3)
                semilogy(x,cumsum(Hfp),x,cumsum(Hfpposterior))
                grid on
                subplot(8,1,4)
                semilogy(x,cumsum(Hne),x,cumsum(Hneposterior))
                grid on
                
                subplot(8,1,5)
                semilogy(x,cumsum(Hfl),x,cumsum(Hflposterior))
                grid on
                subplot(8,1,6)
                semilogy(x,cumsum(Hfi),x,cumsum(Hfiposterior))
                grid on
                subplot(8,1,7)
                semilogy(x,cumsum(Hfc),x,cumsum(Hfcposterior))
                grid on
                subplot(8,1,8)
                semilogy(x,cumsum(HL),x,cumsum(HLposterior))
                grid on
            end
        else
            if (onlyN)
                plot(x,cumsum(HN),x,cumsum(HNposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
            else
                subplot(8,1,1)
                plot(x,cumsum(HN),x,cumsum(HNposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,2)
                plot(x,cumsum(HNs),x,cumsum(HNsposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,3)
                plot(x,cumsum(Hfp),x,cumsum(Hfpposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,4)
                plot(x,cumsum(Hne),x,cumsum(Hneposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,5)
                plot(x,cumsum(Hfl),x,cumsum(Hflposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,6)
                plot(x,cumsum(Hfi),x,cumsum(Hfiposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,7)
                plot(x,cumsum(Hfc),x,cumsum(Hfcposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
                subplot(8,1,8)
                plot(x,cumsum(HL),x,cumsum(HLposterior))
                grid on
                axis([logRange(1) logRange(3) 0 1]);
            end
        end
    case 4
        % Distance to nearest
        
        [nm,x]=hist(log10(0.54901*lamb),1000);
        nmposterior=hist(log10(0.54901*lamb(consistent)),x);
        nm=nm/sum(nm);
        nmposterior=nmposterior/sum(nmposterior);
        
        if (noClf==0) clf; end
        subplot(2,1,1)
        plot(x,nm,x,nmposterior)
        ytop=max(max(nm),max(nmposterior))*1.4
        hold on
        plot(5*[1 1],[0,ytop],':'); text(5,ytop*.95,'Milky Way','FontSize',9,'Rotation',-90)
        plot(10.667*[1 1],[0,ytop],':'); text(10.667,ytop*.95,'Observable universe','FontSize',9,'Rotation',-90)
        xlabel('Log10 median distance to nearest civ (lightyears)')
        axis([0 max(x) 0 ytop])
        subplot(2,1,2)
        plot(x,cumsum(nm),x,cumsum(nmposterior))
        hold on
        plot(x(find(cumsum(nm)>=0.5,1)),0.5,'b+', x(find(cumsum(nmposterior)>=0.5,1)),0.5,'+r')
        xlabel('Log10 median distance to nearest civ (lightyears)')
        grid on
        axis([0 max(x) 0 1])
        
        %         subplot(3,1,3)
        %          f=find(consistent);
        %          distCDF=sum((1-exp(-4*pi*10.^(3*x)./lamb(1:100))))/100;
        %          distCDFposterior=sum((1-exp(-4*pi*10.^(3*x)./lamb(f(1:100)))))/100;
        %          plot(x,distCDF,x,distCDFposterior)
        %          grid on
        % axis([0 max(x) 0 1])
        % xlabel('Log10 median distance to nearest civ (lightyears)')
        %
        
    case 5
        if (noClf==0) clf; end
        plot(fl(consistent),log10N(consistent)-fl(consistent),'b.','MarkerSize',1)
        hold on
        plot(fl(~consistent),log10N(~consistent)-fl(~consistent),'r.','MarkerSize',1)
        axis equal
        xlabel('f_l');
        ylabel('N/f_l')
        
    case {10,11,12} % Make fig 1/2 paper
        fig101112
    case 13
        subplot(2,1,1);
        loglog(1,1,'w')
        hold on
        fill(10.^[x(1) x x(end)],10.^[-6 log10(Hfl+1e-6) -6],[0.7 0.7 1])
        fill(10.^[x(1) x x(end)],10.^[-6 log10(Hflposterior+1e-6) -6],[0.7 1 0.7],'FaceAlpha',0.5)
        axis([1e-40 10 1e-5 1]);
        
        subplot(2,1,2)
        loglog(1,1,'w')
        hold on
        fill(10.^[x(1) x x(end)],10.^[-6 log10(HL+1e-7) -6],[0.7 0.7 1])
        fill(10.^[x(1) x x(end)],10.^[-6 log10(HLposterior+1e-7) -6],[0.7 1 0.7],'FaceAlpha',0.5)
        axis([1 1e10 1e-6 1]);
    case 14
        fig14;
    case 15
        fig14P;
end


%% Results
if (printout>0)
    disp('log10N:')
    disp(sprintf('\tUnconditioned median log10N:%f\n\tUnconditioned mean log10N:%f',median(log10N),mean(log10N)))
    disp(sprintf('\tConditioned median log10N:%f\n\tConditioned mean log10N:%f',median(log10N(consistent)),mean(log10N(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(log10N(consistent))-median(log10N),100*10^median(log10N(consistent))/10^median(log10N)))
    
    disp('Ns:')
    disp(sprintf('\tUnconditioned median Ns:%f\n\tUnconditioned mean Ns:%f',median(Ns),mean(Ns)))
    disp(sprintf('\tConditioned median Ns:%f\n\tConditioned mean Ns:%f',median(Ns(consistent)),mean(Ns(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(Ns(consistent))-median(Ns),100*10^median(Ns(consistent))/10^median(Ns)))
    
    disp('fp:')
    disp(sprintf('\tUnconditioned median fp:%f\n\tUnconditioned mean fp:%f',median(fp),mean(fp)))
    disp(sprintf('\tConditioned median fp:%f\n\tConditioned mean fp:%f',median(fp(consistent)),mean(fp(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(fp(consistent))-median(fp),100*10^median(fp(consistent))/10^median(fp)))
    
    disp('ne:')
    disp(sprintf('\tUnconditioned median ne:%f\n\tUnconditioned mean ne:%f',median(ne),mean(ne)))
    disp(sprintf('\tConditioned median ne:%f\n\tConditioned mean ne:%f',median(ne(consistent)),mean(ne(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(ne(consistent))-median(ne),100*10^median(ne(consistent))/10^median(ne)))
    
    disp('fl')
    disp(sprintf('\tUnconditioned median fl:%f\n\tUnconditioned mean fl:%f',median(fl),mean(fl)))
    disp(sprintf('\tConditioned median fl:%f\n\tConditioned mean fl:%f',median(fl(consistent)),mean(fl(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(fl(consistent))-median(fl),100*10^median(fl(consistent))/10^median(fl)))
    
    disp('fi')
    disp(sprintf('\tUnconditioned median fi:%f\n\tUnconditioned mean fi:%f',median(fi),mean(fi)))
    disp(sprintf('\tConditioned median fi:%f\n\tConditioned mean fi:%f',median(fi(consistent)),mean(fi(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(fi(consistent))-median(fi),100*10^median(fi(consistent))/10^median(fi)))
    
    disp('fc')
    disp(sprintf('\tUnconditioned median fc:%f\n\tUnconditioned mean fc:%f',median(fc),mean(fc)))
    disp(sprintf('\tConditioned median fc:%f\n\tConditioned mean fc:%f',median(fc(consistent)),mean(fc(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(fc(consistent))-median(fc),100*10^median(fc(consistent))/10^median(fc)))
    
    disp('L')
    disp(sprintf('\tUnconditioned median L:%f\n\tUnconditioned mean L:%f',median(L),mean(L)))
    disp(sprintf('\tConditioned median L:%f\n\tConditioned mean L:%f',median(L(consistent)),mean(L(consistent))))
    disp(sprintf('\tMedian shift: %f\n\tPercentage: %f', median(L(consistent))-median(L),100*10^median(L(consistent))/10^median(L)))
    
    disp(sprintf('Probability we are alone in the galaxy\n\tUnconditioned: %f\n\tConditioned:%f', sum(log10N<=0)/N, sum(log10N(consistent)<=0)/Nconsistent))
    disp(sprintf('Probability we are alone in the observable universe\n\tUnconditioned: %f\n\tConditioned:%f', sum(log10N<=-12)/N, sum(log10N(consistent)<=-12)/Nconsistent))
    
    disp(sprintf('Distance to nearest other civ (ly):\n\tMedian unconditioned:%f\n\tMedian conditioned:%f', 0.54901*median(lamb),0.54901*median(lamb(consistent))))
    
end

% c=clock; save(sprintf('Runs/fermi_%d_%d_%d%d%d%d%d',dataSource,posterior,c(1:5)))