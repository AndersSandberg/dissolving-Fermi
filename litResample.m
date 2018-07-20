%% Data
xR=[
    1
10
20
20
20
10
1
10
4/0.7
2.7/0.7
2
10
10
10
5
1
7
10
7
10
40
50
7
10
5/0.7
2
16
    ];

xfp=[0.2
0.5
0.06
1
1
0.25
0.5
0.25
0.05
1
0.05*0.5
0.1
0.04
0.1
0.333
0.333
0.025
0.1
0.25
0.01
1
0.01
1
0.45
0.1
0.1
0.2
0.05
0.2
0.5
0.5
0.4
0.5
0.01
0.1
0.5
0.1
0.5
1
1
0.1
1
0.55
1
1
1
1
    ];


xne=[ 0.05
1
5
1
0.3
0.055
1
1
1
0.1
1
0.01
2
1
0.1
0.05
0.01
2
2
8.60E-04
0.02
0.05
10
3
1
0.012
0.003
0.012
1
1.00E-03
2.80E-06
0.5
3
1
2
1
0.5
2
1
1
1
0.1
0.1
0.3
0.1
0.5
0.14
0.4
0.1
1
0.65
0.74
0.1
0.2
0.2
0.064
    ];

xfl=[1
1
1
1
1.00E-18
2.00E-01
0.2
1
1
1
1.00E-10
0.2
0.01
1.00E-30
0.3333
0.3333
0.5
1
0.01
1
0.01
0.1
0.2
1
1
1
0.5
2.00E-07
1
0.5
0.5
0.5
1
2.00E-17
1
1
0.0005
0.05
1];

xfi=[1
0.1
1
1
1
0.1
1
0.01
1
0.5
0.5
0.1
0.1
6.60E-03
0.75
0.3
1
1
0.05
1
0.01
0.1
0.01
1
0.2
0.5
0.2
1
1.00E-02
0.1
0.05
0.25
1.00E-05
1.00E-30
1
    ];

xfc=[0.1
0.2
0.1
0.5
0.5
0.1
0.1
0.1
1
0.1
1
0.1
0.25
0.1
0.1
1
1
0.5
1
0.5
0.1
0.01
0.2
0.8
0.2
0.5
0.2
0.1
0.5
1
1
1
1
    ];

xL=[6500
1.00E+03
1.00E+08
1.00E+05
100000
1.00E+06
1.00E+07
1.00E+06
100
1.00E+08
1.00E+08
3.00E+06
1.00E+08
1.00E+02
100
1.00E+09
100
1.00E+09
1.00E+05
1.00E+09
1.00E+04
45
9.00E+07
1.00E+08
1.00E+06
1.00E+06
300000
420.6
1.50E+02
1.00E+07
1.00E+07
500
10000
200
10000
1000
1.00E+04
1.00E+09
1.00E+06
10000
65
500
550
5.50E+03
1.00E+08
    ];

Nestimates=[78000
1000
1.00E+08
1.00E+06
1.00E+06
1.00E+08
1.00E+06
100
150
0.003
1.00E+06
27000
4.00E+03
1.00E+06
1.13
30000
1
400
0.0125
2
40000
500
50
0.0003
85
];

if (OutlierRemoval>0)
    xR = sort(xR); xR=xR((OutlierRemoval+1):(length(xR)-OutlierRemoval));
    xfp = sort(xfp); xfp=xfp((OutlierRemoval+1):(length(xfp)-OutlierRemoval));
    xne = sort(xne); xne=xne((OutlierRemoval+1):(length(xne)-OutlierRemoval));
    xfl = sort(xfl); xfl=xfl((OutlierRemoval+1):(length(xfl)-OutlierRemoval));
    xfi = sort(xfi); xfi=xfi((OutlierRemoval+1):(length(xfi)-OutlierRemoval));
    xfc = sort(xfc); xfc=xfc((OutlierRemoval+1):(length(xfc)-OutlierRemoval));
    xL = sort(xL); xL=xL((OutlierRemoval+1):(length(xL)-OutlierRemoval));
end

% Bootstrap resampling
if (BootStrap)
    xR=xR(ceil(rand(length(xR),1)*length(xR)));
    xfp=xfp(ceil(rand(length(xfp),1)*length(xfp)));
    xne=xne(ceil(rand(length(xne),1)*length(xne)));
    xfl=xfl(ceil(rand(length(xfl),1)*length(xfl)));
    xfi=xfi(ceil(rand(length(xfi),1)*length(xfi)));
    xfc=xfc(ceil(rand(length(xfc),1)*length(xfc)));
    xL=xL(ceil(rand(length(xL),1)*length(xL)));
end

Ns=log10(xR(ceil(rand(N,1)*length(xR))));
fp=log10(xfp(ceil(rand(N,1)*length(xfp))));
ne=log10(xne(ceil(rand(N,1)*length(xne))));
fl=log10(xfl(ceil(rand(N,1)*length(xfl))));
fi=log10(xfi(ceil(rand(N,1)*length(xfi))));
fc=log10(xfc(ceil(rand(N,1)*length(xfc))));
L=log10(xL(ceil(rand(N,1)*length(xL))));

if (plotType==6)
    clf
    subplot(4,2,1)
    makesubplot(xR,showstats)
    xlabel('log(N_s)')
   %  axis([-1 3 0 .5])
    
    subplot(4,2,2)
    makesubplot(xfp,showstats)
    xlabel('log(f_p)')
   % axis([-20 0 0 .5])
    
    subplot(4,2,3)
    makesubplot(xne,showstats)
    xlabel('log(n_e)')
    % axis([-6 2 0 .5])
     
    subplot(4,2,4)
    makesubplot(xfl,showstats)
    xlabel('log(f_l)')
    %    axis([-20 0 0 .5])
        
    subplot(4,2,5)
    makesubplot(xfi,showstats)
    xlabel('log(f_i)')
     %   axis([-20 0 0 .5])
        
    subplot(4,2,6)
    makesubplot(xfc,showstats)
    xlabel('log(f_c)')
     %   axis([-20 0 0 .5])
        
    subplot(4,2,7)
    makesubplot(xL,showstats)
    xlabel('log(L)')
   %  axis([0 10 0 .5])
     
    stop
end
