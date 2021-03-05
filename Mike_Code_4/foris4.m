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
% 	h_on: external field, boolean, [topdir,'\',betadir,'\',time(1:5),'JGraphs','_beta',num2str(beta),'_',method()jj,'_',time(6:12),'.png']
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
%cd('/home/michaeta/Documents/MATLAB/Mike_Code_3.1');
%cd('/lustre1/home/michaeta/Mike_Code_3');
%addpath(genpath('/lustre1/home/michaeta/Mike_Code_3'));
addpath(genpath('E:\GitHub\Inv.Is.Models\Mike_Code_4'));
savepath

%%% Storage

%% Is this a test run with disposable files? For deciding where and how to store stuff mostly
testrun = 0 %% Yes test run
%testrun = 0 %% No this is a real run of useful stuff

if testrun == 0
    cd('E:\GitHub\Inv.Is.Models\Mike_Code_4\Warm Storage');  %Changes the working directory to where you want to store your output
elseif testrun == 1
    cd('E:\GitHub\Inv.Is.Models\Disposable Storage');
end

time = datestr(now,'HHMM-ddmmmyy');
disp(time)

mkdir(cd,time);
cd(time);
disp(cd)

topdir = cd;    


%%%% Parameters
    
    %% Trials
    %jn = 100; %| number of Trials
    jn = 5; 

    %test params
    % Nvec = [10];     
    % Tvec = [1e2];   
    % betavec = [0.1];
    
    Nvec =[5,10];					%Nvec = [50,100,150]; %Nvec = [50,100,300]; % Nvec = [100,200,300,400,500];
    
    Tvec = [10,1E2];					%Tvec = [1e3,1e4,1e5]; % Tvec = [1e3,1e4,1e5,1e6]; %| originally preset to be a multiple of the node number:T.calculation: 3*M(numel(M))
    
    betavec = [0.1,0.2]; 			%betavec = [0.4,0.9,1.4];
    
    %% Sparsity measure, used in old "sk" distribution method
    sprsvec = [0.2];						% I need to convert this to a vector as well
    sprs = 0.2; % only here as temp measure
    
    %% Coordination number, what I'm doing here is replacing the previous sparsity measure. This might break some stuff
    %coordvec = [1,2,4,8,12,16];                      
    coordvec = [2,4,8,16];
    
    h_on = 1; 							%% h field genereation
    
    topovec = {2,5,6};
    %topovec = {1,3,4,5,6};		%topologies = {'sk',1,2,3,4,5,6}	%set sk to 0 if you don't want to use it(for now) % make topology the outer looop in the future
                            % if first place is not occupied the ADJSET gets fucked up causing problems further down

%%% for distribtions, see TCS.m for details
    %   I'd like over different or random distribution's later (lol basically monte carlo?)
    %couplings = 4; % "SK" or Mike's Gaussian
    couplings = 5;  % double mean gauss
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

alldir = [time(1:12),'-AllStructs'];
mkdir(cd,alldir);

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

AllStruct = struct;

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
               
               
               for trn = 1:jn
                    
                    % call rng for reproducibility
                    rng(trn);
                    
                    %% Non-changing variable (at least not with trials)
                    AllStruct.(name).name        = name;
                    AllStruct.(name).T           = T;
                    AllStruct.(name).N           = N;
                    AllStruct.(name).beta        = beta;
                    AllStruct.(name).sprsvec     = sprsvec;
                    AllStruct.(name).coordvec   = coordvec;
                    AllStruct.(name).topology    = topovec;
    
                    %AllStruct() = {};
                   OverStruct.list(runs).name        = name;
                   OverStruct.list(runs).T           = T;
                   OverStruct.list(runs).N           = N;
                   OverStruct.list(runs).beta        = beta;
                   OverStruct.list(runs).sprsvec     = sprsvec;
                   OverStruct.list(runs).coordvec    = coordvec;
                   OverStruct.list(runs).topology    = topovec;
                    
                    %%%Forward Ising Topologies and Distributions
                                      
                    Adjset = {}; %make sure these aren't duplicated in the topology creator below
                    Jtoposet = {};
                    hfieldset = {};
                   sigJ = 1/sqrt(N);%---stddev | I used the SK model deviation: beta./sqrt(N) . probably not right, but Nicola was alright with it.
                    % using the std dev of 1/sqrt(N), the beta is already in the system.
                    %changes internally depending on what's chosen
                   for tp = 1:length(topovec) %topology = 1:6 % This doesn't work if you don't have a first place value of any sort in addition to the
                        topo = topovec{tp};
                        [Adj,J,hfield] = TCS(tp,N,c,couplings,beta,J0,sigJ,sprs,h_on)
                        Adjset = [Adjset,Adj];
                        Jtoposet = [Jtoposet,J];
                        hfieldset = [hfieldset,hfield];
                        
                    end
                    
                    
                    JHnorm.Jtopo = Jtoposet; 	% replaces the used J graph if we want to use this topology
                    JHnorm.Htopo = hfieldset;		% Nicola recommends using fixed h for these. iirc setting it too high fucked it up
                    
                    %% These don't change and can probably be rolled up into the non changing variables above
                    AllStruct.(name).topology    = topovec; % starting to duplicate data everywhere
                    AllStruct.(name).couplings   = couplings;
                    AllStruct.(name).c           = c;

                    OverStruct.list(runs).topology               = topovec; % starting to duplicate data everywhere
                    OverStruct.list(runs).couplings              = couplings;
                    OverStruct.list(runs).c                      = c;


                    %% this however does change on each trial
                    AllStruct.(name).Jcontru(trn).trial     = Adjset;
                    AllStruct.(name).Jtru(trn).trial        = Jtoposet;
                    AllStruct.(name).htru(trn).trial        = JHnorm.Htopo;
    
                    OverStruct.list(runs).Jcontru(trn).trial     = Adjset;
                    OverStruct.list(runs).Jtru(trn).trial        = Jtoposet;
                    OverStruct.list(runs).htru(trn).trial        = JHnorm.Htopo;
    
                    
                    %%% Part 2, generate samples (or spike train) S_hat
                    
                    %SstructNorm = Met_Hast_norm(T,N,jn,JHnorm,sprsvec,time,lowdir,beta);
                    SstructNorm = Met_Hast_norm(T,N,size(Jtoposet,2),JHnorm,sprs,time,lowdir,beta); % Replacing the jn with the number of topologies, so each S is a different topology
                    
                    % AllStruct.(name).S = SstructNorm.S_hat; %  %| Adding a fucking ' here so sublime doesn'tlose it's shit %Removed the transpose so it might not work with regular sanity checks, fix that down the line
                    AllStruct.(name).S(trn).trial = SstructNorm; % for working over multiple topologies and multiple output S
                    OverStruct.list(runs).S(trn).trial = SstructNorm;
                    
                                    
                    %%% Part 3, inference and model select
                    %% Bulso Likelihood Estimator
                    
                    %[LLH, jta] = BLLH2(T,N,beta,sprs,h_on,AllStruct.(name).S,Adjset,jta,jtatot); % more or less depreciated since I can parfor the lower part
                    [LLH, jta] = PBLLH(T,N,beta,sprs,h_on,SstructNorm,Adjset,Jtoposet,jta,jtatot);
                    
                                        
                    AllStruct.(name).BLLH(trn).trial = LLH;
                    OverStruct.list(runs).BLLH(trn).trial = LLH;
                    save([overdir,'\',time(1:12),'-OverStruct.mat'],'OverStruct','-v7.3'); 
                end
                runs = runs + 1; %stays out of trial loop, for measuring the runs per other parameters
            end

            %save([cd,'\',time(1:5),'AllStruct_sprs',num2str(sprs),'_beta',num2str(beta),'_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',time(6:12),'.mat'],'AllStruct','-v7.3');
            save([alldir,'\',time(1:5),'AllStruct_C',num2str(c),'_beta',num2str(beta),'_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',time(6:12),'.mat'],'AllStruct','-v7.3');
        end
    end
    toc
    %Jstor = JGraphs4(AllStruct,time,beta,topdir,betadir,Tvec,Nvec);
    
end
    

    %Jstor = modelgraphs2(OverStruct,sprsvec,betavec,Nvec,Tvec,topovec);
    %OverStruct(1).Jstor = Jstor;
    %modelgraphs3(Jstor,sprsvec,betavec,Nvec,Tvec,topologies);
    %Jstor2 = modelgraphs4(Jstor,sprsvec,betavec,Nvec,Tvec,topologies);
    %figstor = modelgraphs6(OverStruct);
    %figstorAUROC = modelgraphs7(OverStruct);
    %figstorTFRatio = modelgraphs8(OverStruct);
save([overdir,'\',time(1:12),'-OverStruct_final.mat'],'OverStruct','-v7.3');  %In theory no different and should be deleted



