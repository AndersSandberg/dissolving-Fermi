N=1e7;
x=rand(N,9)*.2;

clf 
subplot(2,1,1)

N=1e11*prod(x');
[n,x]=hist(N,logspace(-8,8,100));
n=n/sum(n);
fill(log10(x),n,[.7 .7 1])
xlabel('log_{10}(N)');
ylabel('Frequency');
hold on
plot(log10(mean(N))*[1 1],[0 max(n)],'r')
plot(log10(median(N))*[1 1],[0 max(n)],'r:')

% 
% % probability of being alone in galaxy
% Palone = (1-(N/1e11)).^(1e11);
% [n2,x2]=hist(Palone,linspace(0,1,150));
% n2=n2/sum(n2);
% subplot(2,1,2)
%  n2=[0 n2 0];
%  x2=[0 x2 1];
%  fill(x2,n2,[.7 .7 1])
% %bar(x2,(n2))
% xlabel('Pr(alone)');
% ylabel('Frequency');
% hold on
% plot(mean(Palone)*[1 1],[0 max(n2)],'r')
% plot(median(Palone)*[1 1],[0 max(n2)],'r:')
% axis([0 1 0 max(n2)*1.1])




% probability of being alone in galaxy
Pone = exp(log(N) + (1e11-1)*log(1-(N/1e11)));
[n2,x2]=hist(Pone,linspace(0,1,150));
n2=n2/sum(n2);
subplot(2,1,2)
 n2=[0 n2 0];
 x2=[0 x2 1];
 fill(x2,n2,[.7 .7 1])
%bar(x2,(n2))
xlabel('Pr(N=1)');
ylabel('Frequency');
hold on
plot(mean(Pone)*[1 1],[0 max(n2)],'r')
plot(median(Pone)*[1 1],[0 max(n2)],'r:')
axis([0 1 0 max(n2)*1.1])