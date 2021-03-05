%%%ForIs2.m
%% Forward ising
%% 07.08.20: key changes are:
%% - No longer using Overstruct as an amalgamation of allstructs, instead it's now a more table like version of our data (perhaps should be stored that way in the future)
%% - trying to set it to work with trials
%% - Removing sanity checks and Ezaki PLLH (see Mike_code_3 if you need to restore them)
%% - making the OverStruct the main storage device
%% - Changing up how things are stored so I can tell the difference bt trial runs and legit runs
%% - fix name is organized/works, actually make it more complex
%% - Making make container for topology thing, which is just a straightforward calculator of the J and topo array. Any loops should be made outside of it
%% - Not adding a sparseness loop (at least not one that can coexist with c value), but it can be done in the future
%% - Option for adjust the "couplings" or the coupling strength distribution, not adding a loop for this but possibility to do so later
%% -
%% -
%% -




%% Input:
%	Nx: number of nodes
%	T: number of observations or "columns"
% 	h_on: external field, boolean, [topdir,'/',betadir,'/',time(1:5),'JGraphs','_beta',num2str(beta),'_',method()jj,'_',time(6:12),'.png']
%      Should be generated but likely we'll have it set to 0 mostly
%		This input is an on/of deciding whether or not one will be generated
%
%% Inputs/Gens left out (so far):
%	alpha:	no idea what this is for, seems used to create sample array M
%   M: 		number of sample to be taken in an experiment, I think.
% 			In Bcode it's  M = alpha*Nx_v
%  			"Number of observations: notice that K (for CV with l1) has to be an exact divisor for M (usually I take K = 5)""

%			Sparsity of the model? In which case I'm renaming this bitch
%	vbeta:	"Inverse Temp" Represents correlation between nodes (I think)
%			Used multiple betas in code
%	theta: 	coupling strength range


clear all;


%%% File Structure and Storage
%cd(cd)
addpath(genpath('/lustre1/home/michaeta/Mike_Code_4'));
cd('/lustre1/home/michaeta/Mike_Code_4/Warm Storage');

time = datestr(now,'HHMM-ddmmmyy');
disp(time)

mkdir(cd,time);
cd(time);
disp(cd)

topdir = cd;    


%%%% Parameters
    
    %% Trials
    %jn = 100; %| number of Trials
    jn = 200; 

    %test params
    % Nvec = [10];     
    % Tvec = [1e2];   
    % betavec = [0.1];
    
    Nvec =[30,50,70,80,100,120]		%Nvec = [50,100,150]; %Nvec = [50,100,300]; % Nvec = [100,200,300,400,500];
    
Tvec = [1e3,1e4,1e5,1e6] % Tvec = [1e3,1e4,1e5,1e6]; %| originally preset to be a multiple of the node number:T.calculation: 3*M(numel(M))
    
    %betavec = [0.01]
    %betavec = [0.05]
    %betavec = [0.1]
    %betavec = [0.2]  % Beta Set 1     		
    %betavec = [0.3]
    betavec = [0.4]
    %betavec = [0.5]  % Beta Set 2    
    %betavec = [0.6]
    %betavec = [0.7]
    %betavec = [0.8]  % Beta Set 3    
    %betavec = [1.0]
    %betavec = [1.2]      % Beta Set 4
    %betavec = [1.4]
    %betavec = [1.6]      % Beta Set 5
    
    %% Sparsity measure, used in old "sk" distribution method
    sprsvec = [0];						% I need to convert this to a vector as well
    sprs = 0; % only here as temp measure
    
    %% Coordination number, what I'm doing here is replacing the previous sparsity measure. This might break some stuff
    coordvec = [1,2,4,8,12,16]                      
    %coordvec = [2,4,6];
    
    %h_on = 1; 							%% h field genereation
    h_on = 0; 

    topovec = {1,2,3,4,5,6,7}
    %topovec = [2,5]; %need tomake not a cell
    %topovec = {1,3,4,5,6};		%topologies = {'sk',1,2,3,4,5,6}	%set sk to 0 if you don't want to use it(for now) % make topology the outer looop in the future
                            % if first place is not occupied the ADJSET gets fucked up causing problems further down

%%% for distribtions, see TCS.m for details
    %   I'd like over different or random distribution's later (lol basically monte carlo?)
    couplings = 4; % "SK" or Mike's Gaussian
    %couplings = 5;  % double mean gauss
    %couplings = 1; %---Gaussian
    %couplings = 2; %---Delta Function
    %couplings = 3; %---Double Delta Function. on average have the same amount of +/-1 values and will make sure all the same weight values
    J0 = 1;         %---"The Mean" but not exactly clear what it actually is. Keep at 1 for now, changes internally if different distribution is used
    %sigJ = 1/sqrt(N);%---stddev | I used the SK model deviation: beta./sqrt(N) . probably not right, but Nicola was alright with it.
                    % using the std dev of 1/sqrt(N), the beta is already in the system.
                    %changes internally depending on what's chosen

%% For displaying and monitoring the number of trials that are being run
jta = 1;                % Random name. Our measure of how many trials are run so far
jtatot = length(coordvec)*length(betavec)*length(Tvec)*length(Nvec)*length(topovec)*jn

runs = 1;   % for indexing the trials ran for sprs, beta, T , N
            % keep out of the trials loop

lowdir = [time(1:12)];
mkdir(cd,lowdir);
    
paramdir = [time(1:12),'-parameters'];
mkdir(cd,paramdir);
save([paramdir,'\',time(1:12),'-parameters.mat'],'Tvec','Nvec','betavec','topovec','jn','sprsvec','h_on','time','coordvec','topdir','jtatot','-v7.3');                        



overdir = [time(1:12),'-OverStructs'];
mkdir(cd,overdir);

imgdir = [time(1:12),'-Figures'];
mkdir(cd,imgdir);

OverStruct = struct;

OverStruct.Nvec 	   = Nvec ;
OverStruct.Tvec 	   = Tvec ;
OverStruct.betavec 	   = betavec ;
OverStruct.topologies   = topovec;
OverStruct.jn 		    = jn 	;
OverStruct.sprsvec      = sprsvec;
OverStruct.h_on 	   = h_on ;
OverStruct.topdir 	 = topdir ;
OverStruct.time 	   = time ;

%AllStruct = struct;

for Ti = 1:length(Tvec)
    
    T = Tvec(Ti);
    
    tic
    for Ni = 1:length(Nvec)
        N = Nvec(Ni);
        
        %for  Si = 1:length(sprsvec)
        for  Ci = 1:length(coordvec)
         
            %sprs = sprsvec(Si);
            c = coordvec(Ci);
        
        
            for Bi = 1:length(betavec)
                
                beta = betavec(Bi);
                
                cd(topdir)


               %name = ['St',num2str(Si),'Bt',num2str(Bi),'N',num2str(N),'T1E',num2str(log10(T))];
               name = ['T',num2str(Ti),'N',num2str(N),'St',num2str(Ci),'Bt',num2str(Bi)]; % have to keep in St for now even though it reps the coord no
                   
                %% Non-changing variable (at least not with trials)
                %AllStruct.(name).name        = name;
                %AllStruct.(name).T           = T;
                %AllStruct.(name).N           = N;
                %AllStruct.(name).beta        = beta;
                %AllStruct.(name).sprsvec     = sprsvec;
                %AllStruct.(name).coordvec   = coordvec;
                %AllStruct.(name).topology    = topovec;
    
                %AllStruct() = {};
                OverStruct.list(runs).name        = name;
                OverStruct.list(runs).T           = T;
                OverStruct.list(runs).N           = N;
                OverStruct.list(runs).beta        = beta;
                OverStruct.list(runs).sprsvec     = sprsvec;
                OverStruct.list(runs).coordvec    = coordvec;
                OverStruct.list(runs).topology    = topovec;
                    
                %AllStruct.(name).topology    = topovec; % starting to duplicate data everywhere
                %AllStruct.(name).couplings   = couplings;
                %AllStruct.(name).c           = c;

                OverStruct.list(runs).topology               = topovec; % starting to duplicate data everywhere
                OverStruct.list(runs).couplings              = couplings;
                OverStruct.list(runs).c                      = c;
               
                for tp = 1:length(topovec) %topology = 1:6 % This doesn't work if you don't have a first place value of any sort in addition to the
                                    
                    topo = topovec{tp};
                   
                    %%%Forward Ising Topologies and Distributions
                    JHnorm = struct;
                    %Adjset = {}; %make sure these aren't duplicated in the topology creator below
                    %Jtoposet = {};
                    %hfieldset = {};
                    sigJ = 1/sqrt(N);   %---stddev | I used the SK model deviation: beta./sqrt(N) . probably not right, but Nicola was alright with it.
                                        % using the std dev of 1/sqrt(N), the beta is already in the system.
                                        %changes internally depending on what's chosen
                   
                    parfor trn = 1:jn
                        % call rng for reproducibility
                        rng(trn);

                        [Adj,J,hfield] = TCS(tp,N,c,couplings,beta,J0,sigJ,sprs,h_on);
                        JHnorm(trn).Adjset= Adj;
                        JHnorm(trn).Jtopo = J;
                        JHnorm(trn).Htopo = hfield;
                        
                        
                        %Adjset = [Adjset,Adj];
                        %Jtoposet = [Jtoposet,J];
                        %hfieldset = [hfieldset,hfield];
                    end
                    
                    OverStruct.list(runs).Jcontru(tp).topo = {JHnorm.Adjset};
                    OverStruct.list(runs).Jtru(tp).topo    = {JHnorm.Jtopo};
                    OverStruct.list(runs).htru(tp).topo    = {JHnorm.Htopo};
                                        
                    %JHnorm.Jtopo = Jtoposet; 	% replaces the used J graph if we want to use this topology
                    %JHnorm.Htopo = hfieldset;		% Nicola recommends using fixed h for these. iirc setting it too high fucked it up
                    
                    %AllStruct.(name).Jcontru(tp).topo     = Adjset;
    %                %AllStruct.(name).Jtru(tp).topo        = Jtoposet;
    %                %AllStruct.(name).htru(tp).topo        = hfieldset;
    %
                    %OverStruct.list(runs).Jcontru(tp).topo     = Adjset;
                    %OverStruct.list(runs).Jtru(tp).topo        = Jtoposet;
                    %OverStruct.list(runs).htru(tp).topo        = hfieldset;

                    %AllStruct.(name).Jcontru(tp).topo          = JHnorm.Adjset;
                    %AllStruct.(name).Jtru(tp).topo             = JHnorm.Jtopo;
                    %AllStruct.(name).htru(tp).topo             = JHnorm.Htopo;
                    %OverStruct.list(runs).Jcontru(tp).topo     = JHnorm.Adjset;
                    %OverStruct.list(runs).Jtru(tp).topo        = JHnorm.Jtopo;
                    %OverStruct.list(runs).htru(tp).topo        = JHnorm.Htopo;
    
                    
                    %%% Part 2, generate samples (or spike train) S_hat
                    
                    %SstructNorm = Met_Hast_norm(T,N,jn,JHnorm,sprsvec,time,lowdir,beta);
                    %SStruct = Met_Hast_norm(T,N,size(Jtoposet,2),JHnorm,sprs,time,lowdir,beta); % Replacing the jn with the number of topologies, so each S is a different topology
                    SStruct = Met_Hast_norm(T,N,jn,JHnorm,sprs,time,lowdir,beta);
                    
                    %AllStruct.(name).S(tp).topo = SStruct; % for working over multiple topologies and multiple output S
                    OverStruct.list(runs).S(tp).topo = SStruct;
                    
                                    
                    %%% Part 3, inference and model select
                    %% Bulso Likelihood Estimator
                    %[LLH, jta] = BLLH2(T,N,beta,sprs,h_on,AllStruct.(name).S,Adjset,jta,jtatot); % more or less depreciated since I can parfor the lower part
                    [LLH,statvecs,stats,jta] = PBLLH2(T,N,beta,sprs,h_on,SStruct,JHnorm,jta,jtatot);
                    
                    %AllStruct.(name).BLLH(tp).topo = LLH;
                    %AllStruct.(name).BLLH(tp).statvecs = statvecs;
                    %AllStruct.(name).BLLH(tp).stats = stats;

                    OverStruct.list(runs).BLLH(tp).topo = LLH;
                    OverStruct.list(runs).BLLH(tp).statvecs = statvecs;
                    OverStruct.list(runs).BLLH(tp).stats  = stats;
                    save([overdir,'\',time(1:12),'-OverStruct.mat'],'OverStruct','-v7.3');
                    disp(saving) 
                end
                runs = runs + 1; %stays out of trial loop, for measuring the runs per other parameters
            end

            %save([cd,'/',time(1:5),'AllStruct_sprs',num2str(sprs),'_beta',num2str(beta),'_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',time(6:12),'.mat'],'AllStruct','-v7.3');
            %save([alldir,'/',time(1:5),'AllStruct_C',num2str(c),'_beta',num2str(beta),'_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',time(6:12),'.mat'],'AllStruct','-v7.3');
        end
    end
    toc
end
    
save([overdir,'\',time(1:12),'-OverStruct_final.mat'],'OverStruct','-v7.3');  %In theory no different and should be deleted

% comp will set what tpe of figures we want
% comp = 1;   % 1. Perconerr v beta
% 
% figstor = multiallgraph3(OverStruct,comp,topdir,time);
% save([overdir,'/',time(1:12),'-figstor.mat'],'figstor','-v7.3'); 

