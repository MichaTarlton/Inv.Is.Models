%%%ForIs.m
%% Forward ising

%% Creating forward Ising solution from scratch

%% We need to create S, the array that is the outcome of our forward Ising solution and meant to fed into the reverse Ising
%% Also need magnetizations and Complexity


%%% Create an array size: nodes(Nx) x observations(T)

%% Input:
%	Nx: number of nodes
%	T: number of observations or "columns"
% 	h_on: external field, boolean, 
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

T = 1e5; %| Presetting the T.calculation: 3*M(numel(M))
N = 50;
jn = 10; %| Trials
sparsity = 0;

time = datestr(now,'HHMM-ddmmmyy');

%%%  Part 1, generate coupling matricies and h field
%% Gaussian Dist

%% h field genereation
h_on = 1;

%% Decimation algo, decimate random connections in coupling matrix
% Sparsity setting is set from inside JH fnc currently

%JHstruct = JH(N,jn,h_on,sparsity);

JHstruct = JHs(N,jn,h_on,sparsity,time,T); %Turn off if runing the normal one
JHDstruct = JD(N,jn,h_on,sparsity,time,T); %for dimers


%%% Part 2, generate samples (or spike train) S_hat, first using Met_Hast, then using Mean_Field
%% 2.1 Generate courrelation and magnetization from field
%% 2.2 Generate S_hat(s_big in bulso) one S vector at a time, for some length based on M


Sstruct = Met_Hast(T,N,jn,JHstruct,sparsity,time);
SstructDimer = Met_Hast_D(T,N,jn,JHDstruct,sparsity,time);
%%% Part ??? SANITY CHECK

sanity = sanitychk(jn,Sstruct,JHstruct,sparsity,time,T);
sanitydimer = sanitychkdimer(jn,SstructDimer,JHDstruct,sparsity,time,T);

%Plot
figure
scatter(sort([Sstruct.mfinal]),sort([SstructDimer.mfinal]))
hold on 
title({'Magnetizations: mi',['N = ',num2str(N)],['T = ',num2str(T)]})
xlabel('Disconnected')
ylabel('Dimer')

figure
scatter(sort([sanity.tchk]),sort([sanitydimer.tchk]))
hold on 
title({'Magnetization Check: th(h) - mi',['N = ',num2str(N)],['T = ',num2str(T)]})
xlabel('Disconnected')
ylabel('Dimer')


figure
scatter(sort(mean([sanity.chi])),sort(mean([sanitydimer.chi])))
hold on 
title({'Correlation Check: chi = sisj - mimj',['N = ',num2str(N)],['T = ',num2str(T)]})
xlabel('Disconnected')
ylabel('Dimer')

%%%Part 3 (This is actuall part inference)
%% Create array X to regress on Y out of sampled s vectors from S_hat
%% Pull Y vector, Do we generate this as a probablistic output layer as Bulso does (eq. 1 in his paper and pg.7 in my notes) 
%%		He mentions something along these lines in his videochat


