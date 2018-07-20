showcond=1;
LW=2;
FS=18;

mu1=0;
sigma1=10;
mu2=0;
sigma2=2;

theta=4;

air=1.1;

Sx=20;
Sy=10;
s=0.1;
[x,y]=meshgrid(-Sx:s:Sx,-Sy:s:Sy);
marginal1=normpdf(x,mu1,sigma1);
marginal2=normpdf(y,mu2,sigma2);

marginal1post = marginal1.*normcdf(theta-x);
marginal1post = marginal1post/sum(s*marginal1post(1,:));
marginal2post = marginal2.*normcdf(theta-y);
marginal2post = marginal2post/sum(s*marginal2post(:,1));

clf
axes('Position',[0.1,0.1,0.7,0.7])
contourf(x,y,marginal1.*marginal2)
hold on
if (showcond) plot(-Sx:s:Sx,theta-(-Sx:s:Sx),'w','LineWidth',LW); end
set(gca,'XLim',[-Sx Sx]);
set(gca,'YLim',[-Sy Sy]);
xlabel('log_{10}(x_1)','FontSize',FS)
ylabel('log_{10}(x_2)','FontSize',FS)

axes('Position',[0.1,0.85,0.7,0.1])
plot(-Sx:s:Sx,marginal1(1,:),'LineWidth',LW)
hold on
if (showcond) plot(-Sx:s:Sx,marginal1post(1,:),'LineWidth',LW); end
plot(s*sum((-Sx:s:Sx).*marginal1(1,:))*[1 1],[0 max(marginal1post(1,:))*air],'b:','LineWidth',LW)
if (showcond) plot(s*sum((-Sx:s:Sx).*marginal1post(1,:))*[1 1],[0 max(marginal1post(1,:))*air],'r:','LineWidth',LW); end
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'YLim',[0 max(max(marginal1post(1,:)))*air])

axes('Position',[0.85,0.1,0.1,0.7])
plot(marginal2(:,1),-Sy:s:Sy,'LineWidth',LW)
hold on
if (showcond) plot(marginal2post(:,1),-Sy:s:Sy,'LineWidth',LW); end;
plot([0 max(marginal2post(:,1))*air],s*sum((-Sy:s:Sy)'.*marginal2(:,1))*[1 1]','b:','LineWidth',LW)
plot([0 max(marginal2post(:,1))*air],s*sum((-Sy:s:Sy)'.*marginal2post(:,1))*[1 1]','r:','LineWidth',LW)
set(gca,'XLim',[0 max(marginal2post(:,1))*air])

set(gca,'YTick',[])
set(gca,'XTick',[])