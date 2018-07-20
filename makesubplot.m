function makesubplot(X,showstats)
% Make subplots of logs of factors
lx=log10(X);
xmin=min(floor(lx))-1;
xmax=max(ceil(lx))+1;

x=linspace(xmin,xmax,1000);

hold on

if (showstats)
[MUHAT,SIGMAHAT] = normfit(lx)
fill([x(1) x x(length(x))],[0 normpdf(x,MUHAT,SIGMAHAT) 0],'r','FaceAlpha',.2)

[ahat,bhat] = unifit(lx)
fill(x,unifpdf(x,ahat,bhat),'g','FaceAlpha',.2)

text(0.1,0.85,sprintf('\\mu=%2.2f, \\sigma=%2.2f\n|b-a|=%2.2f',MUHAT,SIGMAHAT,(bhat-ahat)),'FontSize',8,'Units','normalized')
end

%n=hist(lx,x); n=n/sum(n);
%bar(x,n)

histogram(lx,'Normalization','pdf','BinWidth',0.25)


%histogram(X,xx,'Normalization','pdf')
%set(gca,'xscale','log'); % scale the x-axis


%axis([xmin xmax 0 1])