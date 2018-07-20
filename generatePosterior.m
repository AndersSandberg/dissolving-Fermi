%% Posterior for Drake equation given non-observation

switch posterior
    case 1
        % Simplest model: perfect detection below certain density, prob detection
        % above
        consistent = (log10N < par1) | (log10N >= par1).*(rand(N,1)<par2);
    case 2
        % Version recalculating based on distance and random site
        %V=pi*1000*50000^2; % Milky way volume cubic light years.
        % Wolfram alpha: 8e12 ly^3, 2.3e11 pc^3
        V=2.3e11; % pc^3
        consistent = exp(-4*pi*dDetect^3*(1-probfaildetect)*(10.^log10N)/V)>rand(N,1);
    case 3
        % Version where NChecked stars have been checked and found empty
        % this is pretty weak!
        logNMW=log10(300e9);
        Pconsistent = (1-10.^min(0,log10N-logNMW)).^NChecked;
        consistent = rand(N,1)<Pconsistent;
        
    case 4
        % colonization model basic: P(consistent)=exp(-N)
        Pconsistent = exp(-10.^log10N);
        consistent = rand(N,1)<Pconsistent;
        
    case 5
%         % Advanced colonization model
%         Tconsistent=(10.^L)<colonyTime;
%         Pconsistent = exp(-10.^log10N); %time galaxy entirely empty
%         Pconsistent = Pconsistent +Tconsistent.*(1-exp(-10.^log10N)).*(1-(10.^(2.5*L)/(3.5*colonyTime^2.5)));
%         %Pconsistent = Pconsistent + (1-Tconsistent).*(colonyTime.*10.^(-L))*(1-1/3.5);
%         Pconsistent = Pconsistent + (1-Tconsistent).*(1-exp(-10.^log10N)).*(1-1/3.5); % UPD
%         consistent = rand(N,1)<Pconsistent;
        
        
        Pconsistent = exp(-10.^log10N); %time galaxy entirely empty
        Tconsistent=(10.^L)<colonyTime; % case when dies out during colonization
       PconsLT = Pconsistent + (1-Pconsistent).*((10.^L).^2.5./(3.5*colonyTime^2.5)); % case when L<T
       Pconsistent = Pconsistent + (1-Pconsistent).*((colonyTime./(10.^L))*(1-1/3.5));
       Pconsistent(Tconsistent)=PconsLT(Tconsistent);
        consistent = rand(N,1)<Pconsistent;
    case 6
        % 1 dark biosphere
        consistent = rand(N,1)<10.^fl;
    case 7
        % searched biospheres
        p=10.^fl;
        consistent = rand(N,1) < nchoosek(nSearched,nFound)*p.^nFound.*(1-p).^(nSearched-nFound);

    case 8
        % Prehistoric intelligence
        consistent = rand(N,1)<10.^fi; 
    case 9
        % Extant aliens
        consistent = log10N>log10(2);
        
    case 10
        % Colonization, but extinction cannot happen before colonytime
%         Tconsistent=(10.^L)<colonyTime;
%         Pconsistent = exp(-10.^log10N); %time galaxy entirely empty
%         Pconsistent = Pconsistent +Tconsistent.*(1-exp(-10.^log10N)).*(1-(10.^(2.5*L)/(3.5*colonyTime^2.5)));
%         %Pconsistent = Pconsistent + (1-Tconsistent).*(colonyTime.*10.^(-L))*(1-1/3.5);
%         Pconsistent = Pconsistent + (1-Tconsistent).*(1-exp(-10.^log10N)).*(1-1/3.5);
%         consistent = (rand(N,1)<Pconsistent) & (Tconsistent==1);
        
        Pconsistent = exp(-10.^log10N); %time galaxy entirely empty
        Tconsistent=(10.^L)<colonyTime; % case when dies out during colonization
       PconsLT = Pconsistent + (1-Pconsistent).*((10.^L).^2.5./(3.5*colonyTime^2.5)); % case when L<T
       Pconsistent = Pconsistent + (1-Pconsistent).*((colonyTime./(10.^L))*(1-1/3.5));
       Pconsistent(Tconsistent)=PconsLT(Tconsistent);
       consistent = (rand(N,1)<Pconsistent) & (Tconsistent==1);
      
        
        
    case 11
        % G-hat cutoff    
        Pconsistent = 1-PK3*(1-PK3success.^(KK3*10.^log10N));
        consistent = rand(N,1)<Pconsistent;
        
    case 12
        % No past civ colonized galaxy
        
        % Version using past stars: this one is suspect!
        %npaststars= (1e10*Ns)/300e9;
        %consistent = rand(N,1)>(npaststars.*10.^log10N) ;
        
        % version using formula in supplement
        consistent = rand(N,1)<exp(-(10.^(log10N-L))*13e9);
        
    case 13 
           % No past civ colonized galaxy, or one is on its way here now
        Pconsistent = exp(-(10.^(log10N-L))*(13e9-colonyTime));
        Pconsistent = Pconsistent+exp(-(10.^(log10N-L))*colonyTime)*(1-1/3.5);
             consistent = rand(N,1)<Pconsistent;
        
    otherwise
        consistent = logical(ones(N,1));
end

% Sigmoid etc.