%%%fi0201.m
%% Full Ising problem

function fi0201(jobname,intbeta)
%load('ERmissparam.mat')
% Total time of  run, measured at the end in hours
totaltime = tic;

% Displays the HPC array number input, originally used here to select beta value
disp(['Array number: ',num2str(intbeta),' Job Name: ', num2str(jobname)])

%%% File Structure and Storage
%cd(cd)
%addpath(genpath('E:\GitHub\Inv.Is.Models\Mike_Code_4.1'));
%cd('E:\GitHub\Inv.Is.Models\Mike_Code_4.1\Warm Storage'); 
%savepath

addpath(genpath('/lustre1/home/michaeta/Mike_Code_4'));
cd(['/lustre1/home/michaeta/',num2str(jobname)]);
%jobdir=convertCharsToStrings(intjobdir)
%cd(jobdir)

%%% Storage
%% Is this a test run with disposable files? For deciding where and how to store stuff mostly
%testrun = 1; %% Yes test run
%testrun = 0 %% No this is a real run of useful stuff

% if testrun == 0
%     cd('E:\GitHub\Inv.Is.Models\Mike_Code_4.1\Warm Storage');  %Changes the working directory to where you want to store your output
% elseif testrun == 1
%     cd('E:\GitHub\Inv.Is.Models\Mike_Code_4.1\Disposable Storage');
% end

%%% Current Time
time = datestr(now,'HHMM-ddmmmyy');
disp(time)

%% for rng
hpct=clock();

seed=hpct(6) * 1000; % Seed with the second part of the clock array.

rng(seed);

%rng('shuffle','philox')
%rndy1 = num2str(randi([1 99],1))
%rng('shuffle','philox')
%rndy2 = num2str(randi([1 99],1))
%rng('shuffle','philox')
%rndy3 = num2str(randi([1 99],1))
%rng('shuffle','philox')
%rndy4 = num2str(randi([1 999],1))

rndy = num2str(randi([100 999],1))

%dirname = [time,'-',rndy]
dirname = [num2str(jobname),'_',num2str(intbeta),'-',time]

mkdir(cd,dirname);
cd(dirname);
disp(cd)

topdir = cd;    


%%%% Parameters
    
    %% Trials
    jn = 100; %| number of Trials
    %jn = 5; 

    h_on = 0 							%% h field genereation
    %% Ok so important to note in the decimation logistical regression code,
    %% having this set to off will fuck things up as in it will not attempt to
    %% infer our h values and return an array which doesn't work with the rest
    %% of my code. Here my h_on represents setting the h values to all zeros,
    %% or having the same distribution as is being used for J.
    
    %test params
    % Nvec = [10];     
    % Tvec = [1e2];   
    % betavec = [0.1];
    
    %Tvec = [10,30]
    Tvec = [100,200]				%Tvec = [1e3,1e4,1e5]; % Tvec = [1e3,1e4,1e5,1e6]; %| originally preset to be a multiple of the node number:T.calculation: 3*M(numel(M))
    %Tvec = [1]

    %Nvec =[120]
    Nvec =[100] 	% N val 1	   %Nvec = [50,100,150]; %Nvec = [50,100,300]; % Nvec = [100,200,300,400,500];
    %Nvec =[50]  % N val 2
    %Nvec =[70]  % N val 3
    %Nvec =[80]  % N val 4
    %Nvec =[100] % N val 5
    %Nvec =[120] % N val 6
   
    betavec = [0.3,1,1.6]
    %betavecint = [1]
    %betavec = betavecint(1)
    %betavec = [0.01]   %betavec value 1
    %betavec = [0.1]    %betavec value 2
    %betavec = [0.2]    %betavec value 3
    %betavec = [0.3]    %betavec value 4
    %betavec = [0.4]    %betavec value 5
    %betavec = [0.5]    %betavec value 6
    %betavec = [0.6]    %betavec value 7 
    %betavec = [0.7]    %betavec value 8
    %betavec = [0.8]    %betavec value 9
    %betavec = [1.0]    %betavec value 10
    %betavec = [1.2]    %betavec value 11
    %betavec = [1.4]    %betavec value 12
    %betavec = [1.6]    %betavec value 13
        
    
    %% Sparsity measure, used in old "sk" distribution method
    sprsvec = [0];						% I need to convert this to a vector as well
    sprs = 0; % only here as temp measure
    
    %% Coordination number, what I'm doing here is replacing the previous sparsity measure. This might break some stuff
    %coordvec = [1,2,4,8,12,16];                      
    coordvec = [50,70,90]
    %coordvecint = [28]
    %coordvec = coordvecint(intbeta)
    
    topo    = 7
    topovec = [7]; %no longer using cell for this, I'll fix this down the line
    %topovec = {5}
    %topovec = {1,2,5};
    %topovec = [2,5]; %need tomake not a cell
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

lnc = length(coordvec)
lnb = length(betavec)
lnt = length(Tvec)
lnn = length(Nvec)
lntop   = length(topovec)

%% For displaying and monitoring the number of trials that are being run
jta = 1;                % Random name. Our measure of how many trials are run so far
jtatot = length(coordvec)*length(betavec)*length(Tvec)*length(Nvec)*length(topovec)*jn

runs = 1;   % for indexing the trials ran for sprs, beta, T , N
            % keep out of the trials loop

% lowdir = [time(1:12)];
% mkdir(cd,lowdir);
    
% paramdir = [time(1:12),'-parameters'];
% mkdir(cd,paramdir);
% save([paramdir,'/',time(1:12),'-parameters.mat'],'Tvec','Nvec','betavec','topovec','jn','sprsvec','h_on','time','coordvec','topdir','jtatot','-v7.3');                        
save([topdir,'/',num2str(jobname),'_',num2str(intbeta),'-',time(1:12),'-parameters.mat'],'Tvec','Nvec','betavec','topovec','jn','sprsvec','h_on','time','coordvec','topdir','jtatot','-v7.3');                        
% 
% alldir = [time(1:12),'-AllStructs'];
% mkdir(cd,alldir);
% 
% overdir = [time(1:12),'-OverStructs_2'];
% mkdir(cd,overdir);
% 
% imgdir = [time(1:12),'-Figures'];
% mkdir(cd,imgdir);

% Making paramstruct---------------------------------------
paramstruct = struct;
paramstruct.t     = [];
paramstruct.n     = [];
paramstruct.b     = [];
paramstruct.c     = [];
paramstruct.top   = [];
prun = 1;
for Ti = 1:lnt
    for Ni = 1:lnn
        for  Ci = 1:lnc
            for Bi = 1:lnb
                 for tp = 1:lntop
                    paramstruct(prun).T     = Tvec(Ti) ;
                    paramstruct(prun).N     = Nvec(Ni) ;
                    paramstruct(prun).beta     = betavec(Bi) ;
                    paramstruct(prun).c     = coordvec(Ci) ;
                    paramstruct(prun).top   = topovec(tp);
prun = prun + 1;
end
end
end                    
end                
end

%----------------------------------------------------------



OverStruct = struct;

OverStruct.Nvec 	   = Nvec ;
OverStruct.Tvec 	   = Tvec ;
OverStruct.betavec 	   = betavec ;
OverStruct.coordvec     = coordvec ;
OverStruct.topologies   = topo;
OverStruct.jn 		    = jn 	;
OverStruct.sprsvec      = sprsvec;
OverStruct.h_on 	   = h_on ;
OverStruct.topdir 	 = topdir ;
OverStruct.time 	   = time ;
% save([overdir,'/',num2str(jobname),'_',num2str(intbeta),'-OverStruct.mat'],'OverStruct','-v7.3');
save([topdir,'/',num2str(jobname),'_',num2str(intbeta),'-OverStruct.mat'],'OverStruct','-v7.3');

%AllStruct = struct;

%create a local cluster object
%distcomp.feature( 'LocalUseMpiexec', false ) % highly experimental here
pc = parcluster('local')
%pc = parcluster('threads')
% explicitly set the Job Storage Location to the temp directory that was created in your sbatch script
mkdir(cd,'scratch')
parscratch = [topdir,'/scratch']
pc.JobStorageLocation = parscratch
parpool(pc,20)

%parpool('local',18)
tic
    N       = paramstruct(intbeta).N  ;
    T       = paramstruct(intbeta).T;
    beta    = paramstruct(intbeta).beta  ;
    c       = paramstruct(intbeta).c  ;
    topo    = paramstruct(intbeta).top ;

    tp = 1; %because we are only doing a single topology, probably need to fix down the line
    % origiinally for iterating over multiple topologies, but we are only doingone at a time currently



    %N       = paramstruct(intbeta).n  ;
    %T       = paramstruct(intbeta).t;
    %beta    = paramstruct(intbeta).b  ;
    %c       = paramstruct(intbeta).c  ;
    %topo    = paramstruct(intbeta).top;

    %beta = beta./sqrt(c);
            
    cd(topdir)


    %name = ['St',num2str(Si),'Bt',num2str(Bi),'N',num2str(N),'T1E',num2str(log10(T))];
    name = ['T',num2str(T),'N',num2str(N),'C',num2str(c),'Bt',num2str(beta)]; % have to keep in St for now even though it reps the coord no

    
               
    OverStruct.list(runs).name        = name;
    OverStruct.list(runs).T           = T;
    OverStruct.list(runs).N           = N;
    OverStruct.list(runs).beta        = beta;
    OverStruct.list(runs).sprsvec     = sprsvec;
    OverStruct.list(runs).coordvec    = coordvec;
    %OverStruct.list(runs).topology    = topovec;
    OverStruct.list(runs).topology    = topo;
    OverStruct.list(runs).couplings   = couplings;
    OverStruct.list(runs).c           = c;
                    
                   
    %%%Forward Ising Topologies and Distributions
    JHnorm = struct;
                   
    parfor trn = 1:jn
        % call rng for reproducibility
        rng(trn);

        [Adj,J,hfield] = TCS4(topo,N,c,couplings,beta,J0,sprs,h_on);
        JHnorm(trn).Adjset= Adj;
        JHnorm(trn).Jtopo = J;
        JHnorm(trn).Htopo = hfield;
    end
                    
    OverStruct.list(runs).Jcontru(tp).topo = {JHnorm.Adjset};
    OverStruct.list(runs).Jtru(tp).topo    = {JHnorm.Jtopo};
    OverStruct.list(runs).htru(tp).topo    = {JHnorm.Htopo};
                        
    %JHnorm.Jtopo = Jtoposet; 	% replaces the used J graph if we want to use this topology
    %JHnorm.Htopo = hfieldset;		% Nicola recommends using fixed h for these. iirc setting it too high fucked it up
        
    %OverStruct.list(runs).Jcontru(tp).topo     = Adjset;
    %OverStruct.list(runs).Jtru(tp).topo        = Jtoposet;
    %OverStruct.list(runs).htru(tp).topo        = hfieldset;

    %OverStruct.list(runs).Jcontru(tp).topo     = JHnorm.Adjset;
    %OverStruct.list(runs).Jtru(tp).topo        = JHnorm.Jtopo;
    %OverStruct.list(runs).htru(tp).topo        = JHnorm.Htopo;
    
    
    %%% Part 2, generate samples (or spike train) S_hat
    
    %SstructNorm = Met_Hast_norm(T,N,jn,JHnorm,sprsvec,time,lowdir,beta);
    %SStruct = Met_Hast_norm(T,N,size(Jtoposet,2),JHnorm,sprs,time,lowdir,beta); % Replacing the jn with the number of topologies, so each S is a different topology
    SStruct = Met_Hast_norm(T,N,jn,JHnorm,sprs,time,beta);
    
    %OverStruct.list(runs).S(tp).topo = SStruct; % not going to keep these for future runs because the memory is taking a beating
    
                    
    %%% Part 3, inference and model select

    %% Bulso Likelihood Estimator
    %[LLH, jta] = BLLH2(T,N,beta,sprs,h_on,AllStruct.(name).S,Adjset,jta,jtatot); % more or less depreciated since I can parfor the lower part
    %[LLH,statvecs,stats,jta] = PBLLH2(T,N,beta,sprs,h_on,SStruct,JHnorm,jta,jtatot);
    [LLH,statvecs,stats,jta] = PBLLH6(T,N,tp,beta,c,h_on,SStruct,JHnorm,jta,jtatot,topdir,intbeta);
    

    %OverStruct.list(runs).BLLH(tp).topo = LLH;
    OverStruct.list(runs).BLLH(tp).statvecs = statvecs;
    OverStruct.list(runs).BLLH(tp).stats  = stats;


    %save([overdir,'/',time(1:12),'-OverStruct.mat'],'OverStruct','-v7.3'); 
    save([topdir,'/',num2str(jobname),'_',num2str(intbeta),'-OverStruct.mat'],'OverStruct','-v7.3'); 
    disp('saving overstruct')

     runs = runs + 1; %stays out of trial loop, for measuring the runs per other parameters
    

           
toc

delete(gcp('nocreate'))  
%save([overdir,'/',time(1:12),'-OverStruct_final.mat'],'OverStruct','-v7.3');  %In theory no different and should be deleted
eval(['OverStruct_' num2str(intbeta) '= OverStruct']);
save(['/lustre1/home/michaeta/overstructs/',num2str(jobname),'_',num2str(intbeta),'-OverStruct.mat'],'OverStruct','-v7.3');

endtime = toc(totaltime);
disp(['Total Time: ', num2str(endtime./3600)])
end
% comp will set what tpe of figures we want
% comp = 1;   % 1. Perconerr v beta
% 
% figstor = multiallgraph3(OverStruct,comp,topdir,time);
% save([overdir,'\',time(1:12),'-figstor.mat'],'figstor','-v7.3'); 

