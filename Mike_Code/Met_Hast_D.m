%%%%Met_Hast.m

%%% Mike's implementation of Metroplis Hasting
%%% Goal: Output S_hat (the spike train or observations) for coupling matrix J and field h

%% Inputs:
%			T: Number of observations or datapoints
%			N: number of nodes or neurons
% 			J: input coupling field
%			h: external field
% 		Other inputs may be applied in the future
%
%% Outputs: 
%			S_hat: The entire spike train, known as S_hat in Nicola's version
%			mfinal: mfinal in Nicola code, final magnetizations
%			Cfinal: final correlation matrix
% 		Need to create structs which allow recording of in between outputs for debugging purposes
%		
%% NOTES.
%	- No functionality for 0/1 values applied in this
%
%

function SstructDimer = Met_Hast_D(T,N,jn,JHstruct,sparsity,time)

	
% initial settings, copied wholesale from Nicola
Nequil = 1e4*N;				%that means 1e4 upgrades of the whole system for equilibrium (1e4 swips)
Int = 10*N; 				%steps interval between two different measurements (10 swips) 				
Navg = T*Int;			%| end state number of iterations post equilibrium reached
Energy = zeros(1,Nequil+Navg); 	%|

SstructDimer = struct('S_hat',{},'mfinal',{},'Cfinal',{},'mequil',{},'Cequil',{});


for i = 1:jn

	% Set J and H from JHstruct
	J =JHstruct(i).Jsparse;
	h = JHstruct(i).Hsparse;
		
	% initial configuration & representation 
	S = double(rand(1,N) > 0.5);
	S = 2*S-1;

	% initial magnetization and correlations
	m = S';							%| Complex transpose the vector of random charge
	C = S'*S;						%| Conjugate Transpose matrix
	Npoints = 1e2;					%| setting 100 points? why an arbitrar number? 
	Step = floor(Nequil/Npoints);	%| 
	err_m = zeros(1,Npoints);		%| Zero vector size of npoints
	err_C = zeros(1,Npoints);		%| Zero vector size of npoints
	m_old = m;						%| Store previous complex transpose
	C_old = C;						%| Store previous Conjugate transpose matrix	
	l = 0;							%| Not sure what this is #?		
	
			%%%%%%%%%%%%%%%%%%% INITIAL LOOP FOR EQUILIBRIUM %%%%%%%%%%%%%%%%%%%%%%%
		
		S_eq = zeros(T,N);			      		%| Zeros matrix
		mu_eq = 0;						      		%| #?
      		
		for j = 1:Nequil 				      		%| ah so just creating a high amount of loops to go through?
		     
		    % store last points at intervals of Int
		    if j > Nequil-Navg			%| what does the view function do here. Notsure what the whole thing here does or what Navg looks like, in my tests this is a high numbe, way above J , I thought it would be a negative number
		 
		        if mod(j-(Nequil-Navg),Int) == 0	%| 
		            mu_eq = mu_eq + 1;				%|
		            S_eq(mu_eq,:) = S;				%|
		        end
		    end
		    
		    k = randi(N);          					% Select a spin at random uniform distribution
		    E = -0.5*S*J*S' - S*h';       			% Energy of starting configuration | What does J come in as, wtf does this look like?
		    Energy(j) = E;							% |Store energy value 
		    St = S;									% |Store S, to what purpose #?
		
		    
		    St(k) = -S(k); 						%|Flipping spin for -1/1 mode
		    
		
		    Et = -0.5*St*J*St' - St*h';   %Energy of the trial configuration | simila
		    DE = Et - E;							%| Energy is a scalar, this takes the difference of Energy 
		    
		    %METROPOLIS ALGORITHM
		    % Comparing the energyconfiguration of the two different states 
		
		    if DE<0                   				%| This tests the value of DE to see if it's an automatic accept                           
		        S = St;               				%| If it's not accepted, it will move to random selection        
		    elseif DE>=0              				%| Random Selection acceptance is very low        
		        X = rand;             				%|         
		        if X < exp(-DE)       				%|               
		            S = St;           				%|           
		        end                   				%|   
		    end                       
		    
		    % updating magnetizations and correlations
		    m = m + (S'-m)/(j+1);
		    C = C + (S'*S-C)/(j+1);
		        
		    % storing updates differences at some steps to check convergence
		    if mod(j,Step) == 0
		        l = l + 1;
		        err_m(l) = sqrt((m - m_old)'*(m - m_old)/N);
		        err_C(l) = sqrt(2*sum(sum(triu((C - C_old).^2,1)))/(N*(N-1)));
		        m_old = m;
		        C_old = C;
		    end
		    
		end


		%%%%%%%%%%%%%%%%%%%%%% LOOP AFTER EQUILIBRIUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%| "At the bottom of the valley"
		%|
		%|This is for collecting data once at the area of highest probability
		%| 
		
		S_hat = zeros(T,N);
		mu = 0;
		for j = Nequil+1: Nequil + Navg
		    
		    % collect data at intervals of Int
		    if mod(j,Int) == 0
		        mu = mu + 1;
		        S_hat(mu,:) = S;
		    end
		    
		    k = randi(N);          %Select a spin at random, uniform distribution
		    E = -0.5*S*J*S' - S*h';       %Energy of starting configuration
		    Energy(j) = E;
		    St = S;
		    
		    % flipping spin k
		    St(k) = -S(k);
		    
		    Et = -0.5*St*J*St' - St*h';   %Energy of the trial configuration
		    DE = Et - E;
		    
		    if DE<0                        %METROPOLIS ALGORITHM
		        S = St;
		    elseif DE>=0
		        X = rand;
		        if X < exp(-DE)
		            S = St;
		        end
		    end
		    
		end

	%%%%%%%%%%%%%%%%%%% END OF LOOP AFTER EQUILIBRIUM %%%%%%%%%%%%%%%%%%%%%%%%
	
	% final estimates of magnetisation and correlation
	SstructDimer(i).S_hat = S_hat;
	SstructDimer(i).mfinal = mean(S_hat);
	SstructDimer(i).Cfinal = S_hat'*S_hat/T;
	SstructDimer(i).mequil = mean(S_eq);
	SstructDimer(i).Cequil = S_eq'*S_eq/T;
    
    disp(['Saving ',num2str(i)])
	save(['SstructDimer_N',num2str(N),'_T',num2str(T),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time,'.mat'],'SstructDimer','-v7.3');
	disp(['End S run ',num2str(i)])
end


%%save(['SstructDimer_N',num2str(N),'_T',num2str(T),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time,'.mat'],'SstructDimer');

end