
 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
%% 10/06/20 - sketching the code for simulations
% 
% Goal: network reconstruction in Ising models by implementing the 
%       approach of Bulso et al. 2019 within the pseudo-likelihood 
%       approximation. Testing on and comparing with other approaches
%       on synthetic data. Studying networks found with this approach
%       from real (or more realistic) neural data
%
%% 1st: create a function that does network reconstruction, that is a
% function that 
%   - takes as input: the matrix S of size M times N where M is the number
%   of samples and N is the size of the network. Each element in S is a
%   binary variable coded as (-1,1) which can represent for instance the
%   activity of a neuron in a certain time bin;
%   - return as output: the network that best describe the observed data
%   according to the criterion proposed in Bulso et al. 2019 (as well as
%   according to other criteria in order to compare reconstructions)
%
% Basically the implementation requires the following steps
%   - for each neuron i in the network (for loop)
%       - take the relative row and call it Y = S(i,:)
%       - take all the other rows and call them X = S(not i,:)
%       - run the routine "[w_ML,l_ML,posterior,cost,BestModel,IMAX] =
%       decimation_logistic_Model_Selection(X,Y,field)" with the option
%       field = 1 which indicates that we want to model also a field parameter for
%       each node.
%       - the function returns (among other things) the "BestModel" structure
%       according to different criteria which is a vector of logical values
%       indicating the presence or the absence of a coupling connecting the 
%       "not i" neurons to the "i-th" neuron. Basically this step in the
%       for loop tells us which elements in the row J(i,:) in the connectivity 
%       matrix is different from zero. Remember that J(i,i) = 0 and what
%       BestModel.method tells you is which elements of the J(i,not i) are
%       different from zero (pay attention to the order of the indices!).
%       The actual maximum likelihood estimates of J(i,not i) are contained 
%       in w_ML(best model,:)

%% 2nd: create the data S from a known topology and run the code above
% to see if the reconstruction with our method improves over other
% approaches. Basically you need to:
%   - generate a topology: fully connected, random graph, 2D lattice, ...
%   and choose some parameter such as the parameter strenght, the size of
%   the network, the number of samples to generate. An informed choice
%   based on the literature allows one to compare more easily with other
%   works. I have uploaded some code in the dropbox folder which I have
%   written for generating different topologies, selecting different
%   coupling distribution, dilution, etc... which I used in the previous
%   publication Bulso et al. 2016.
%   - generate data from a given model. Use Metropolis-Hasting algorithm
%   - run the network reconstruction code several times
%   - compare the quality of reconstruction across different methods (e.g.
%   BIC, AIC, l1, ours)... BIC and AIC already implemented in the 
%   "decimation_logistic_Model_Selection" routine, whereas in the routine
%   "decimation_MS_log_reg_LocIsing.m" at line 173 you can find the
%   implementation of the l1 approach with a K fold CV method for selecting
%   the regularizer
% The goal of this section is to compare our method with others and see if
% it achieves better performances in network reconstruction when varying
% different parameters such as the topology, the strength of the
% connections, the size of the network, ... In particular, given the
% results in Bulso et al. 2019, I would expect this approach to outperform
% others when sample size is not extremely large and the network is highly
% connected. You can also use code in other files in the folder in order 
% to analyse data or for any other subtask in the project (of course, in
% this case you need to extract relevant part and adapt to your purposes)

%% 3rd: now you have the machine and you can use it for studying the networks
% inferred from real (or more realistic) neural data. Possible options are: 
% 1) data from the Nature paper Stensola et al. 2012 (just ask me about that 
%   when you are at this stage and I will give you the data)
% 2) data from Yoanna Sandvig for instance 
% 3) also biologically plausible synthetic data might be interesting to
% study

%% THESIS EXPECTED CONTENTS (time left is only 3 months)
%
% I think that your thesis should have the first part done, as much as you
% can of the second and a bit of the third part. 
