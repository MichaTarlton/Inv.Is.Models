%%%ForIs.m
%% Forward ising

%% Creating forward Ising solution from scratch

%% We need to create S, the array that is the outcome of our forward Ising solution and meant to fed into the reverse Ising
%% Also need magnetizations and Complexity


%%% Create an array size: nodes(Nx) x observations(T)

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

%% Iteration tracking variables:
% 	jn: number of coupling matricies or configuration we want to make

%% Storage structures
% 	Jstruct


%% Might break these out into different scripts/modules

%% Generated (inputs to next step):
%	J: array of connection valules for array
%	h: external field, decides how likely node is to fire in abcense of other influence. vector of values which decide chance of node firing in lack of of other influence

%%% Part 0, So first take input parameters for overall function
% For this, N is number of nodes, 
% 			jn is the number of coupling matrixes we want to create, 
%			hon is whether or not we want to create exernal field
% Possibly add option for distribution chosen
clear all;

Tval = [1e3,1e4];
Nval = [30,40];

%Tval = [1e3,1e4,1e5]; %| Presetting the T.calculation: 3*M(numel(M))
%Nval = [50,100,150];

%Tval = [1e3,1e4,1e5]; %| Presetting the T.calculation: 3*M(numel(M))
%Nval = [50,100,300];

% Tval = [1e3,1e4,1e5,1e6]; %| Presetting the T.calculation: 3*M(numel(M))
% Nval = [100,200,300,400,500];

betavec = [0.1]; 
%betavec = [0.4,0.9,1.4]; 
jn = 1; %| Trials
sparsity = 0;
h_on = 1; %% h field genereation

cd('E:\GitHub\Inv.Is.Models\Mike_Code_3'); %Only for my particular instance

time = datestr(now,'HHMM-ddmmmyy');
mkdir(cd,time);
cd(time);
topdir = cd;

for Bi = 1:length(betavec)
	cd(topdir)

	beta = betavec(Bi);

	betadir = ['Beta_',num2str(beta)];
	mkdir(cd,betadir);
	cd(betadir)


Snamevec = {};

AllStruct = struct;

for Ti = 1:length(Tval)
    T = Tval(Ti);

for Ni = 1:length(Nval)
    N = Nval(Ni);
    
lowdir = [time(1:5),'parameters_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_beta',num2str(beta),'_',time(6:12)];
mkdir(cd,lowdir);
save([lowdir,'\',time(1:5),'parameters_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_beta',num2str(beta),'_',time(6:12),'.mat'],'T','N','jn','sparsity','h_on','time','-v7.3');

name = ['N',num2str(N),'T1E',num2str(log10(T))];

AllStruct.(name) = {};
AllStruct.(name).name = name;
AllStruct.(name).T = T;
AllStruct.(name).N = N;


%%%  Part 1, generate coupling matricies and h field
%% Gaussian Dist


%% Decimation algo, decimate random connections in coupling matrix
% Sparsity setting is set from inside JH fnc currently


%JHnorm = JH(N,jn,h_on,sparsity,time,T,lowdir,beta);
% AllStruct.(name).Jtru = JHnorm.Jsparse;
% AllStruct.(name).htru = JHnorm.Hsparse;

%JHdiscon = JHs(N,jn,h_on,sparsity,time,T,lowdir); %for disconnected J
%JHdimer = JD(N,jn,h_on,sparsity,time,T,lowdir); 	%for dimers
%JHferr = JF(N,jn,h_on,sparsity,time,T,lowdir); 	%for ferromagnetic lattice


%% Additional topologies via Nicola
%% Maybe put this in it's own container
% topo_on = 1;			% Turns on this module for parts below
% 		topology = 1;  	%---cayley tree with coordination number c 
% 		topology = 2;  	%---fully connected topology 
% 		topology = 3;  	%---indipendent pair topology
% 		topology = 4;  	%---2D Ising lattice 
 	topology = 5;  		%---Erdos Reyni random graph 
% 		topology = 6;  	%--- Star network
	c = 3; 		  		%--- Coordination Number: Average number of conenctions each node has. c = 2,3,4,6,8 . How many children each node generates.
 	couplings = 1;		%---Gaussian  
% 		couplings = 2	%---Delta Function 
% 		couplings = 3	%---Double Delta Function. on average have the same amount of +/-1 values and will make sure all the same weight values
	J0 = 1; 				%---"The Mean" but not exactly clear what it actually is. Keep at 1 for now
	sigJ = beta./sqrt(N);%---stddev | I used the SK model deviation: beta./sqrt(N) . probably not right, but Nicola was alright with it.

	Adj = set_topology(topology,N,c);
   	Jtopo = set_couplings(couplings,beta,J0,sigJ,Adj);

JHnorm.Jsparse = Jtopo; 	% replaces the used J graph if we want to use this topology
JHnorm.Hsparse = 0.4 .* ones(1,N);		% Nicola recommends using fixed h for these. iirc setting it too high fucked it up

AllStruct.(name).Jcontru = Adj;
AllStruct.(name).Jtru = Jtopo;
AllStruct.(name).htru = JHnorm.Hsparse;

AllStruct.(name).couplings = couplings;
AllStruct.(name).topology =  topology;
AllStruct.(name).c = c;






%%% Part 2, generate samples (or spike train) S_hat, first using Met_Hast, then using Mean_Field
%% 2.1 Generate courrelation and magnetization from field
%% 2.2 Generate S_hat(s_big in bulso) one S vector at a time, for some length based on M

SstructNorm = Met_Hast_norm(T,N,jn,JHnorm,sparsity,time,lowdir,beta);
AllStruct.(name).S = SstructNorm.S_hat'; % ' %| Adding a fucking ' here so sublime doesn'tlose it's shit % 

%SstructDisc =  Met_Hast_Disc(T,N,jn,JHdiscon,sparsity,time,lowdir,beta);
%SstructDimer = Met_Hast_D(T,N,jn,JHdimer,sparsity,time,lowdir,beta);
%SstructFerr = Met_Hast_F(T,N,jn,JHferr,sparsity,time,lowdir);



%%% Part ??? SANITY CHECK

%sanitynorm = sanitychknorm(jn,SstructNorm,JHnorm,sparsity,time,T,N,lowdir,beta);
%AllStruct.(name).Jmf = sanitynorm.mfJ;
%AllStruct.(name).hmf = sanitynorm.mfh;
%AllStruct.(name).Jtap	=	sanitynorm.tapJ;
%AllStruct.(name).htap	=	sanitynorm.taph;
%AllStruct.(name).Jplmf	=	sanitynorm.plJmf;
%AllStruct.(name).hplmf	=	sanitynorm.plhmf;
%
%sanitydisc = sanitychkdiscon(jn,SstructDisc,JHdiscon,sparsity,time,T,lowdir);
%sanitydimer = sanitychkdimer(jn,SstructDimer,JHdimer,sparsity,time,T,lowdir);
%sanityferr = sanitychkferr(jn,SstructFerr,JHferr,sparsity,time,T,lowdir);

%% PLLH from Ezaki
% woven into main loop instead of being an add-on like below
%PLLHout = pfunc_02_Inferrer_PL(SstructNorm.S_hat',name,time,beta);
%AllStruct.(name).Jpllh = PLLHout.J;
%AllStruct.(name).hpllh = PLLHout.h';


%% Bulso Likelihood Estimator 
% create For loop for each observation step against all others
% no actually put it in it's own conatiner
tic
LLH = BLLH(T,N,h_on,AllStruct.(name).S);
toc

AllStruct.(name).BLLH = LLH;

AllStruct.(name).conerr = sum(double(not((LLH(1).Jcon - Adj) == 0)),'all');






%% Plot

% Check to see if either of these two funcs are needed

% Probably could break this out into a side thing

%Graphs(SstructDisc,sanitydimer,sanitydisc,N,T,lowdir,time);

%Graphs(SstructNorm,sanitynorm,N,T,'Normal Distribution')
%Graphs(SstructDisc,sanitydisc,N,T,'Disconnected')
%Graphs(SstructDimer,sanitydimer,N,T,'Independent Pairs')
%Graphs(SstructFerr,sanityferr,N,T,'Ferromagnetic')

%Jgraphs(JHnorm,sanitynorm,N,T,lowdir,time);
% Jgraphs(JHdiscon,sanitydisc,N,T,lowdir,time);
% Jgraphs(JHdimer,sanitydimer,N,T,lowdir,time);
%Jgraphs(JHferr,sanityferr,N,T,lowdir,time);



%%%Part 3 (This is actuall part inference)
%% Create array X to regress on Y out of sampled s vectors from S_hat
%% Pull Y vector, Do we generate this as a probablistic output layer as Bulso does (eq. 1 in his paper and pg.7 in my notes) 
%%		He mentions something along these lines in his videochat

%% Actually built some analtic infferrence methods into the sanit checks

Snamevec = {Snamevec,name};
end
end
save([lowdir,'\',time(1:5),'AllStruct_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_beta',num2str(beta),'_',time(6:12),'.mat'],'AllStruct','-v7.3');


%Jstor = JGraphs3(AllStruct,time,Snamevec,beta,topdir,betadir,Tval,Nval);

%% PLLH from Ezaki
% This has been incorporated above, this can probably be removed
% this is temporary and should be done away with into it's own script
% Sstructvec = [SN50T3,SN50T4,SN50T5,SN100T3,SN100T4,SN100T5,SN300T3,SN300T4,SN300T5];
% Snamevec = ["SN50T3","SN50T4","SN50T5","SN100T3","SN100T4","SN100T5","SN300T3","SN300T4","SN300T5"];
% 
% PLLHout = Max_edits_pfunc_02_Inferrer_PL(Sstructvec,Snamevec,lowdir,time,beta);

% Jstor = JGraphs2(JHN50T3,JHN50T4,JHN50T5,JHN100T3,JHN100T4,JHN100T5,JHN300T3,JHN300T4,JHN300T5,...
%     sanityN50T3,sanityN50T4,sanityN50T5,sanityN100T3,sanityN100T4,sanityN100T5,sanityN300T3,sanityN300T4,sanityN300T5,...
%     PLLH,...
%     time,name,beta)


end
