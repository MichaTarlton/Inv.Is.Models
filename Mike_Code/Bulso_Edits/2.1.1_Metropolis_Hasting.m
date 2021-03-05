%%% GLAUBER DYNAMICS SPIN 0/1 EVOLUTION PERFORMED WITH METROPOLIS MONTECARLO ALGORITHM %%%
%INPUT : NDATA, NSPINS, COUPLINGS, FIELDS, representation (optionally view, dafault is off) 
%OUTPUT : SPIN CONFIGURATION (matrix Ndata X N)
%25/01/2018 w.r.t the previous code for "Metropolis Evolution", I have
%cleaned the code, improved it a bit and make it work for 0/1 data and -1/1
%data
% 06/07/2018 improvement w.r.t previous code "Metropolis" for visualisation
% of large sistem, checking saturation and convergence


function S_big = Metropolis_Hasting(Ndata,N,J,h,representation,vw)
%				 Metropolis_Hasting(3*M(numel(M)),Nx,beta*W,beta*H,rep)
%	  			 M = vector of sample sizes , Nx = number of inputs , beta = vbeta(i); vbeta = distribution of input featurs according to Ising with parameters , 
%       		J is connections, and h is externaøl field I'm guessing, represenation sets to -1/1 mode, and vw have no idea, doesn't appearused though.	




% plot results
try view = vw; catch
    view = false;
end





% initial settings
Nequil = 1e4*N;%that means 1e4 upgrades of the whole system for equilibrium (1e4 swips)		| Not certain #? 10k updates of the whole system, expecting the biggest lost of the system. The number of steps expected to reach the "bottom" or equilibrium
Int = 10*N; %steps interval between two different measurements (10 swips) 					| Ah ok only for a single run of N. The seperation of samples bt Metroplis polls, we don't want them to be too correlated

Navg = Ndata*Int;																			%|Number of Data points times #?
Energy = zeros(1,Nequil+Navg);




% initial configuration & representation (0/1 false; -1/1 true)
S = double(rand(1,N) > 0.5); 	%|Creates a random vector of 0,1s  the size of N (Nx), Double converts a logical array to numerical


if ~representation				%|Converts whole of S vector from 0,1s to -1,1s, unless representation = true
    mini = 0; maxi = 1;
elseif representation
    S = 2*S-1;
    mini = -1; 
    maxi = 1;
end




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

S_eq = zeros(Ndata,N);			      		%| Zeros matrix
mu_eq = 0;						      		%| #?
      		
for j = 1:Nequil 				      		%| ah so just creating a high amount of loops to go through?
     
    % store last points at intervals of Int
    if j > Nequil-Navg && view 				%| what does the view function do here. Notsure what the whole thing here does or what Navg looks like, in my tests this is a high numbe, way above J , I thought it would be a negative number
 
        if mod(j-(Nequil-Navg),Int) == 0	%| 
            mu_eq = mu_eq + 1;				%|
            S_eq(mu_eq,:) = S;				%|
        end
    end
    
    k = randi(N);          					% Select a spin at random uniform distribution
    E = -0.5*S*J*S' - S*h';       			% Energy of starting configuration | What does J come in as, wtf does this look like?
    Energy(j) = E;							% |Store energy value 
    St = S;									% |Store S, to what purpose #?

    if representation   % flipping spin k
        St(k) = -S(k); 						%|Flipping spin for -1/1 mode
    else
        St(k) = 1-S(k); 					%|Flipping spin for 0/1 mode
    end

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

S_big = zeros(Ndata,N);
mu = 0;
for j = Nequil+1: Nequil + Navg
    
    % collect data at intervals of Int
    if mod(j,Int) == 0
        mu = mu + 1;
        S_big(mu,:) = S;
    end
    
    k = randi(N);          %Select a spin at random, uniform distribution
    E = -0.5*S*J*S' - S*h';       %Energy of starting configuration
    Energy(j) = E;
    St = S;
    if representation   % flipping spin k
        St(k) = -S(k);
    else
        St(k) = 1-S(k); 
    end
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
mfinal = mean(S_big);
Cfinal = S_big'*S_big/Ndata;
mequil = mean(S_eq);
Cequil = S_eq'*S_eq/Ndata;

%NOTE: since I want to calculate the averages using indipendent
%measurements, I let the system evolve Int=10 steps between each measurament

if view
    
    %energy
    figure;
    semilogx(1:1e3:Nequil+Navg,Energy(1:1e3:Nequil+Navg))
    hold on;line([Nequil,Nequil],[0,max(Energy)]);
    title('Energy')
    xlabel('single spin update')
    ylabel('Energy')
    
    %magnetizations
    figure;
    subplot(1,2,1)
    semilogy(err_m)
    hold on;fplot(@(x) sqrt(1/Ndata),'k:'); % order of magnitude of std 1/sqrt(T)
    xlim([0,numel(err_m)]);
    ylabel('magnetization')
    xlabel('updates')
    subplot(1,2,2)
    scatter(mequil,mfinal);
    hold on;fplot(@(x) x); % order of magnitude of std 1/sqrt(T)
    hold on;fplot(@(x) x + sqrt(1/Ndata),'k:');    
    hold on;fplot(@(x) x - sqrt(1/Ndata),'k:');
    axis([mini,maxi,mini,maxi]);
    ylabel('magnetization recorded')
    xlabel('magnetization accumulated equilibrium')
    
    %correlations
    figure;
    subplot(1,2,1)
    semilogy(err_C);
    xlim([0,numel(err_C)]);
    hold on;fplot(@(x) sqrt(1/Ndata),'k:'); % order of magnitude of std 1/sqrt(T)
    ylabel('correlation')
    xlabel('updates')
    subplot(1,2,2)
    scatter(Cequil(:),Cfinal(:));
    hold on;fplot(@(x) x); % order of magnitude of std 1/sqrt(T)
    hold on;fplot(@(x) x + sqrt(1/Ndata),'k:');    
    hold on;fplot(@(x) x - sqrt(1/Ndata),'k:');
    axis([mini,maxi,mini,maxi]);
    ylabel('correlations recorded')
    xlabel('correlations last point equilibrium')
    
    % magnetizations
    width = 0.5;
    figure;
    bar(1:N,mfinal,width)
    title('Magnetization')
    ylabel('m(i)')
    xlabel('i')
    
    % correlations
    figure;
    bar3(Cfinal,width)
    title('Correlations')
    zlabel('C(i,j)')
    xlabel('i')
    ylabel('j')
end
