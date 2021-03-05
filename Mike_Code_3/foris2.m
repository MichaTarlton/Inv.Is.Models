%%%ForIs2.m
%% Forward ising
%% 01.06.20 Trying to rebuild this from scratch to streamline our parameter input and review process

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
%cd(cd)
cd('E:\GitHub\Inv.Is.Models\Mike_Code_3');  %Changes the working directory to where you want to store your output
time = datestr(now,'HHMM-ddmmmyy');
mkdir(cd,time);
cd(time);
topdir = cd;

% Parameters
jn = 1; %| number of Trials

%test params
% Nvec = [10];     
% Tvec = [1e2];   
% betavec = [0.1];

Nvec = [10,20,30];					%Nvec = [50,100,150]; %Nvec = [50,100,300]; % Nvec = [100,200,300,400,500];

Tvec = [1e3,1e4];					%Tvec = [1e3,1e4,1e5]; % Tvec = [1e3,1e4,1e5,1e6]; %| originally preset to be a multiple of the node number:T.calculation: 3*M(numel(M))

betavec = [0.3,0.4]; 			%betavec = [0.4,0.9,1.4];

sprsvec = [0];						% I need to convert this to a vector as well

h_on = 1; 							%% h field genereation

topovec = {1,3,4,5,6};		%topologies = {'sk',1,2,3,4,5,6}	%set sk to 0 if you don't want to use it(for now) % make topology the outer looop in the future
                        % if first place is not occupied the ADJSET gets fucked up causing problems further down


jta = 1;                % Random name. Our measure of how many trials are run so far
jtatot = length(sprsvec)*length(betavec)*length(Tvec)*length(Nvec)*length(topovec)*jn

runs = 1;   % for indexing the trials ran for sprs, beta, T , N
    




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

for  Si = 1:length(sprsvec)
 
    sprs = sprsvec(Si);

    AllStruct = struct;

    for Bi = 1:length(betavec)
        
        beta = betavec(Bi);
        
        cd(topdir)
        betadir = ['Beta_',num2str(beta)];
        mkdir(cd,betadir);
        cd(betadir)
        
               
        for Ti = 1:length(Tvec)
            
            T = Tvec(Ti);
            
            for Ni = 1:length(Nvec)
                
                N = Nvec(Ni);
                
                %name = ['N',num2str(N),'T1E',num2str(log10(T))];
                name = ['St',num2str(Si),'Bt',num2str(Bi),'N',num2str(N),'T1E',num2str(log10(T))];
                
                
                %lowdir = [time(1:5),'parameters_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_beta',num2str(beta),'_',time(6:12)];
                lowdir = [time(1:5),'parameters_',name,'_',time(6:12)];

                mkdir(cd,lowdir);
                %save([lowdir,'\',time(1:5),'parameters_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_beta',num2str(beta),'_',time(6:12),'.mat'],'T','N','beta','topologies','jn','sprsvec','h_on','time','name','-v7.3');
                save([lowdir,'\',time(1:5),'parameters_',name,'_',time(6:12),'.mat'],'T','N','beta','topovec','jn','sprsvec','h_on','time','name','-v7.3');
                
                
                AllStruct.(name) = {};
                AllStruct.(name).name = name;
                AllStruct.(name).T = T;
                AllStruct.(name).N = N;
                AllStruct.(name).beta = beta;
                AllStruct.(name).sprsvec = sprsvec;
                AllStruct.(name).topology = topovec;

                %AllStruct() = {};
                AllStruct.list(runs).name = name;
                AllStruct.list(runs).T = T;
                AllStruct.list(runs).N = N;
                AllStruct.list(runs).beta = beta;
                AllStruct.list(runs).sprsvec = sprsvec;
                AllStruct.list(runs).topology = topovec;
                
                
                
                %%%Forward Ising Topologies
                % This is a mess rn, we want to pull this to the outside tbh. ALso need to better streamline my SK model with the other topologies
                
                Adjset = {}; %make sure these aren't duplicated in the topology creator below
                Jtoposet = {};
                
                if topovec{1} == 'sk'
                    %% Mike's implementation of SK model (I think)
                    JHnorm = JH(N,jn,h_on,sprs,time,T,lowdir,beta); % SK model
                    % AllStruct.(name).Jtru = JHnorm.Jsparse;
                    % AllStruct.(name).htru = JHnorm.Hsparse;
                    
                    Jtoposet = [Jtoposet,JHnorm.Jsparse2]; % also cribbed from below, but should be in duplicate, this add the (single for now) topology to the front of the Jtoposet thingy.
                    % Will need to be fixed for multiple trials due to how the index, you can
                    % see the difference in methastnorm.m line 38 onward
                    Adjset = [Adjset,JHnorm.Jcon]; % this is for only the connectivity matricies
                end
                
                %% Topologies for sanity check purposes
                %JHdiscon = JHs(N,jn,h_on,sprsvec,time,T,lowdir); %for disconnected J
                %JHdimer = JD(N,jn,h_on,sprsvec,time,T,lowdir); 	%for dimers
                
                
                
                %% Additional topologies via Nicola
                %% Maybe put this in it's own container
                %% I'm making a for loop here, but it really should be put into it's own container.
                % this fails to support trials make into it's own container
                %% Consider doing so on another revision
                
                
                % topo_on = 1;			% Turns on this module for parts below. This was supposed to be a switch for Nicola's topologies, but isn't doing much atm
                % 		topology = 1;  	%---cayley tree with coordination number c
                % 		topology = 2;  	%---fully connected topology
                % 		topology = 3;  	%---independent pair topology
                % 		topology = 4;  	%---2D Ising lattice
                %	 	topology = 5;  	%---Erdos Reyni random graph
                % 		topology = 6;  	%--- Star network
                %       Add +1 if you include other topologies from above
                c = 2; 		  		%--- Coordination Number: Average number of conenctions each node has. c = 2,3,4,6,8 . How many children each node generates.
                couplings = 1;	        %---Gaussian
                %couplings = 2;     	%---Delta Function
                % 		couplings = 3;	%---Double Delta Function. on average have the same amount of +/-1 values and will make sure all the same weight values
                J0 = 1; 				%---"The Mean" but not exactly clear what it actually is. Keep at 1 for now
                sigJ = beta./sqrt(N);%---stddev | I used the SK model deviation: beta./sqrt(N) . probably not right, but Nicola was alright with it.
                
                %Adjset = {};
                %Jtoposet = {};
                
               if topovec{1} == 'sk'
                skon = 2;
                else
                    skon =1;
                end
                    

                for tp = skon:length(topovec) %topology = 1:6 % This doesn't work if you don't have a fist place value of any sort in addition to the
                    topo = topovec{tp};
                    Adj = set_topology(topo,N,c);
                    Jtopo = set_couplings(couplings,beta,J0,sigJ,Adj);
                    
                    Adjset = [Adjset,Adj];
                    Jtoposet = [Jtoposet,Jtopo];
                    %TopStruct.(['topo',num2str(topology)]) = {Adj,Jtopo}; %| I think I overcomplicated it
                end
                
                
                JHnorm.Jtopo = Jtoposet; 	% replaces the used J graph if we want to use this topology
                JHnorm.Htopo = 0 .* ones(1,N);		% Nicola recommends using fixed h for these. iirc setting it too high fucked it up
                
                AllStruct.(name).topology = topovec; % starting to duplicate data everywhere
                AllStruct.(name).Jcontru = Adjset;
                AllStruct.(name).Jtru = Jtoposet;
                AllStruct.(name).htru = JHnorm.Htopo;
                AllStruct.(name).couplings = couplings;
                AllStruct.(name).c = c;

                AllStruct.list(runs).topology = topovec; % starting to duplicate data everywhere
                AllStruct.list(runs).Jcontru = Adjset;
                AllStruct.list(runs).Jtru = Jtoposet;
                AllStruct.list(runs).htru = JHnorm.Htopo;
                AllStruct.list(runs).couplings = couplings;
                AllStruct.list(runs).c = c;



                
                
                
                
                %%% Part 2, generate samples (or spike train) S_hat, first using Met_Hast, then using Mean_Field
                %% 2.1 Generate courrelation and magnetization from field
                %% 2.2 Generate S_hat(s_big in bulso) one S vector at a time, for some length based on M
                
                
                %SstructNorm = Met_Hast_norm(T,N,jn,JHnorm,sprsvec,time,lowdir,beta);
                SstructNorm = Met_Hast_norm(T,N,size(Jtoposet,2),JHnorm,sprs,time,lowdir,beta); % Replacing the jn with the number of topologies, so each S is a different topology
                
                % AllStruct.(name).S = SstructNorm.S_hat; %  %| Adding a fucking ' here so sublime doesn'tlose it's shit %Removed the transpose so it might not work with regular sanity checks, fix that down the line
                AllStruct.(name).S = SstructNorm; % for working over multiple topologies and multiple output S
                AllStruct.list(runs).S = SstructNorm;
                %% Sanity check spike trains
                %SstructDisc =  Met_Hast_Disc(T,N,jn,JHdiscon,sprs,time,lowdir,beta);
                %SstructDimer = Met_Hast_D(T,N,jn,JHdimer,sprs,time,lowdir,beta);
               
                
                
                
                %%% Part 3 SANITY CHECK graphs and inference methods
                % Variational inerence methods are built in here
                
                %sanitynorm = sanitychknorm(jn,SstructNorm,JHnorm,sprs,time,T,N,lowdir,beta);
                %AllStruct.(name).Jmf = sanitynorm.mfJ;
                %AllStruct.(name).hmf = sanitynorm.mfh;
                %AllStruct.(name).Jtap	=	sanitynorm.tapJ;
                %AllStruct.(name).htap	=	sanitynorm.taph;
                %AllStruct.(name).Jplmf	=	sanitynorm.plJmf;
                %AllStruct.(name).hplmf	=	sanitynorm.plhmf;
                
                %sanitydisc = sanitychkdiscon(jn,SstructDisc,JHdiscon,sprs,time,T,lowdir);
                %sanitydimer = sanitychkdimer(jn,SstructDimer,JHdimer,sprs,time,T,lowdir);
               
                
                
                
                %% PLLH from Ezaki
                % woven into main loop instead of being an add-on like below
                %Spllh = SstructNorm.S_hat;
                %PLLHout = pfunc_02_Inferrer_PL(Spllh',name,time,beta);
                %AllStruct.(name).Jpllh = PLLHout.J;
                %AllStruct.(name).hpllh = PLLHout.h'; %|'
                
                
                %% Bulso Likelihood Estimator
                % create For loop for each observation step against all others
                % no actually put it in it's own conatiner
                disp('bllh2')
                tic
                %[LLH, jta] =
                %BLLH2(T,N,beta,sprs,h_on,AllStruct.(name).S,Adjset,jta,jtatot);
                %%Regular bllh, testing PBLLH rn
                [LLH, jta] = PBLLH(T,N,beta,sprs,h_on,AllStruct.(name).S,Adjset,jta,jtatot);
                toc
                
                %Testing Ncolas original code to see if I somehow fucked
                %something up
                % Checked it and nodifference as far as I can tell
%                 disp('Testing Ncolas original code to see if I somehow fucked something up')
%                 
%                 tic 
%                 LLHold = decimation_MS_log_reg_LocIsing_old(N,h_on,T,SstructNorm);
%                 toc
                
                AllStruct.(name).BLLH = LLH;
                AllStruct.list(runs).BLLH = LLH;
                %AllStruct.(name).BLLHold = LLHold;
                %AllStruct.(name).avgconerr = mean([LLH.totconerr]); %vestigial since made individual errors for each penalty method
             runs = runs + 1;   
                
            end
            
            
            
            
            %% Plot
            % Probably could break this out into a side thing
            
            %Graphs(SstructDisc,sanitydimer,sanitydisc,N,T,lowdir,time);
            %
            %Graphs(SstructNorm,sanitynorm,N,T,'Normal Distribution')
            %Graphs(SstructDisc,sanitydisc,N,T,'Disconnected')
            %Graphs(SstructDimer,sanitydimer,N,T,'Independent Pairs')
            %
            %
            %Jgraphs(JHnorm,sanitynorm,N,T,lowdir,time);
            %Jgraphs(JHdiscon,sanitydisc,N,T,lowdir,time);
            %Jgraphs(JHdimer,sanitydimer,N,T,lowdir,time);
            
            
            
            
            
        end
        
        save([cd,'\',time(1:5),'AllStruct_sprs',num2str(sprs),'_beta',num2str(beta),'_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',time(6:12),'.mat'],'AllStruct','-v7.3');
    end
    

    %Jstor = JGraphs4(AllStruct,time,beta,topdir,betadir,Tvec,Nvec);
    OverStruct(Si).AllStruct = AllStruct;
    
end    


    Jstor = modelgraphs2(OverStruct,sprsvec,betavec,Nvec,Tvec,topovec);
    OverStruct(1).Jstor = Jstor;
    %modelgraphs3(Jstor,sprsvec,betavec,Nvec,Tvec,topologies);
    %Jstor2 = modelgraphs4(Jstor,sprsvec,betavec,Nvec,Tvec,topologies);
    %figstor = modelgraphs6(OverStruct);
    %figstorAUROC = modelgraphs7(OverStruct);
    %figstorTFRatio = modelgraphs8(OverStruct);
    save([cd,'\',time(1:5),'OverStruct_',time(6:12),'.mat'],'OverStruct','-v7.3');



