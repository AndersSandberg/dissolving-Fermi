[n,x]=hist(log10N,logRange(1):logRange(2):logRange(3));
%hh=bar(x,n);

HN=n/N;
HNposterior = hist(log10N(consistent),x)/Nconsistent;

HNs = hist(Ns,x)/N;
HNsposterior = hist(Ns(consistent),x)/Nconsistent;
Hfp= hist(fp,x)/N;
Hfpposterior = hist(fp(consistent),x)/Nconsistent;
Hne= hist(ne,x)/N;
Hneposterior = hist(ne(consistent),x)/Nconsistent;
Hfl= hist(fl,x)/N;
Hflposterior = hist(fl(consistent),x)/Nconsistent;
Hfi= hist(fi,x)/N;
Hfiposterior = hist(fi(consistent),x)/Nconsistent;
Hfc= hist(fc,x)/N;
Hfcposterior = hist(fc(consistent),x)/Nconsistent;
HL = hist(L,x)/N;
HLposterior = hist(L(consistent),x)/Nconsistent;