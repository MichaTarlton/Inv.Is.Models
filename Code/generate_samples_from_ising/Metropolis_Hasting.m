%%% SAMPLING NETWORK CONFIGURATIONS WITH METROPOLIS ALGORITHM %%%
% sample configurations of the network using Metropolis_algorithm. This
% is useful both if you simply want to sample from an Ising model and also
% for "boltzmann learning" where the gradient is calculated by sampling

%INPUT : number of data points, parameters (fields and couplings)
%OUTPUT : partition function, magnetisation and correlations

% Note: it works for both -1/1 and 0/1 representations;

function S_big = Metropolis_Hasting(Ndata,J,h,representation,vw)

% plot results
try view = vw; catch
    view = false;
end

% select -1/1 representation if rep = true; 0/1 otherwise
try rep = representation; catch
    rep = true;
end

% settings, initializations and initial configuration
N = numel(h);
Int = 10*N; %steps interval between two different measurements (10 swips)
Navg = Ndata*Int;
Nequil = max(1e4*N,Navg);%that means 1e4 upgrades of the whole system for equilibrium (1e4 swips)
Energy = zeros(1,Nequil+Navg);
S = double(rand(1,N) > 0.5);
if ~rep
    mini = 0; maxi = 1;
elseif rep
    S = 2*S-1;
    mini = -1; maxi = 1;
end

% initial magnetization and correlations
m = S';
C = S'*S;
Npoints = 1e2;
Step = floor(Nequil/Npoints);
err_m = zeros(1,Npoints);
err_C = zeros(1,Npoints);
m_old = m;
C_old = C;
l = 0;

%%%%%%%%%%%%%%%%%%% INITIAL LOOP FOR EQUILIBRIUM %%%%%%%%%%%%%%%%%%%%%%%
S_eq = zeros(Ndata,N);
mu_eq = 0;
for j = 1:Nequil
    
    % store last points at intervals of Int
    if j > Nequil-Navg && view
        if mod(j-(Nequil-Navg),Int) == 0
            mu_eq = mu_eq + 1;
            S_eq(mu_eq,:) = S;
        end
    end
    
    k = randi(N);          %Select a spin at random uniform distribution
    E = -0.5*S*J*S' - S*h';       %Energy of starting configuration
    Energy(j) = E;
    St = S;
    if rep   % flipping spin k
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
    if rep   % flipping spin k
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
