% TCS.m
% Topology and Coupling Strength
% Combining all the parts for Topology, Coupling Strengths, and Sparsity
% Not going to loop anthing into here, instead this is just a straightforward calculator of the J and topo array

function [Adj,J,hfield] = TCS3(topology,N,c,couplings,beta,J0,sprs,h_on)


% inputting these at a higher level
% J0 = 1;			%---"The Mean" but not exactly clear what it actually is. Keep at 1 for now
% sigJ = 1/sqrt(N); %---stddev | I used the SK model deviation: beta./sqrt(N) . probably not right, but Nicola was alright with it.

% couplings = 1;	%---Gaussian, Nicola's oiginal was positives only, I've made it add random negative values, though I'm a bit uncertain the distribution of them is correct 
% couplings = 2;   	%---Delta Function
% couplings = 3;	%---Double Delta Function. on average have the same amount of +/-1 values and will make sure all the same weight values
% couplings = 4;	%---"SK" model, my own addition based on old JH, doesn't use sigJ or J0
% couplings = 5;	%---Gaussian, but my own take on how nicola did it
% 

% topology = 1;  	%---cayley tree with coordination number c
% topology = 2;  	%---fully connected topology
% topology = 3;  	%---independent pair topology
% topology = 4;  	%---2D Ising lattice
% topology = 5;  	%---Erdos Reyni random graph
% topology = 6;  	%--- Star network
% Add +1 if you include other topologies from above
% c = 2;	  		%--- Coordination Number: Average number of conenctions each node has. c = 2,3,4,6,8 . How many children each node generates.


                    
%% Topology and Sparseness

	if topology == 1	%---cayley tree with coordination number c
		[Adj] = generateTree(N,c);
		distr = 1/sqrt(c); % for double delta
		%distr = 1./c; % Bc we are sure there are exactly 

	elseif topology == 2	%---fully connected topology | Adding sprsity
		%Adj = ones(N) - eye(N); % what nicola originally had
		Adj = double(rand(N)> sprs);
		Adj = triu(Adj,1);
		Adj = Adj + Adj';
		distr = 1/sqrt(N-1);

	elseif topology == 3	%---indipendent pair topology
		[Adj] = indipendentPair(N);    
		distr = 1/sqrt(N-1); %Really not sure what to make this distribution, going to leave it as sk model (with N-1)

	elseif topology == 4	%---2D Ising lattice
		 [Adj,c] = lattice(N,2);
		 distr = 1/sqrt(c); 
	
	elseif topology == 5	%---Erdos Reyni random graph
		p = c/(N-1); %probability of two spins being connected
		Adj = rand(N,N) < p;
		Adj = triu(Adj,1);
		Adj = Adj + Adj';
		distr = 1/sqrt(c);
		%  comment1: in the limit of Nspin tend to infty this graph gives the same
		%  result as fully connected with dilution 1-p
		%  comment2: the average vertex degree is c = p*(Nspin-1), for instance
		%  if Nspin = 50, in order to have c = 3 and 4 you need p = 0.0612 and p = 0.0816 
		%  if Nspin = 20 you can use p = 0.2 for c = 4 and p = 0.15 for c = 3; p =
		%  0.05 for c = 1
		%  if Nspin = 41, then p = 3/40 for c = 3 and p = 1/10 for c = 4
	
	elseif topology == 6
		%--- Star network
		[Adj] = star(N); 
		distr = 1/sqrt(N-1);
	elseif topology == 7
		%--- Small World - Watts Strogatz Graph
		p = c/(N-1); %Probability that an edge will switch to a random node, \beta in the Watts Strogatz model. 
						% I'm setting the probability to be the same as the probability for the average vertex degree
						% So it will scale with N and c
						% When p=0 then we have a ring graph, when p=1 it is a random graph
		% see inside why I'm dividing it by 2, it's doubling the mean node degree inside
		graphy = WattsStrogatz(N,c./2,p); %this is passedout in matlab graph format. I need to convert to array

		Adj = full(adjacency(graphy));
        distr = 1/sqrt(c);
	end


%%coupling strength
	
	if couplings == 1 %---gaussian
		%distr = sigJ;
		J = randn(N);
	 	J = beta.*Adj.*(J0 + sigJ.*(J + J')/sqrt(2)); % none of these go to negatives, causing trouble in the S generation
		% my std dev isn't any good here, it tops above 1 easily. Setting J0=0.5 makes it better but consider using the nthroot(N,3)
	 
		%% Adding the negative values
		R = ones(N);
		R(rand(N) > 0.5 ) = -1;
		J = J.*R;
	
		% Let's say I want to uniformly distribute on an interval later (e.g. interval -5,5):
		% r = -5 + (5+5)*rand(10,1)
	
		% Or maybe normally distribute with an actual mean similar to below
		% normrnd(0.5,1/2.*nthroot(N,3),[N,N]) % possibly try even higher root
		% Also fuck with this distribution: (1/nthroot(N,2)) + (1/nthroot(N,3))
	
		% My own take
		% R = ones(N);
		% R(rand(N) > 0.5 ) = -1;
		% J = beta*Adj.*normrnd(0.5,1/2.*nthroot(N,3),[N,N])
	
	elseif couplings == 2 %---delta function
	 	%distr = 1;
	 	J = beta.*J0.*Adj;

	%---double delta function
	elseif couplings == 3
		 J = randn(N);
		 %distr = 1;
		 %J =  beta.*J0.*Adj.*sign(J + J'); 
		 J =  beta.*distr.*J0.*Adj.*sign(J + J'); 
		% when use delta functions make sure that the variance of the prior
		% contains the value of beta (i.e. if sig(3) is 1 is not good to choose
		% beta = 5, instead you could pick beta = 0.5

	elseif couplings == 4 %| Mike's gaussian from the "SK" model gen: JHnorm, difference is mine is meaned around 0
		%R = double(normrnd(0,beta./sqrt(N),[N,N])); %Originally
		%distr = 1/nthroot(N,2);
		%distr = sigJ;
		J0 = 0;
		R = normrnd(J0,distr,[N,N]);
		R = triu(R,1) + triu(R,1)';
		J = beta.*Adj.*R;

	elseif couplings == 5 % My own take on Nicola's gaussian
	 	J0 = 0.5;	%mean at .5 and then random negatives added, so it's a double mean
	 	%distr = sigJ;
	 	
	 	R = normrnd(J0,distr,[N,N]);
	 	negs = ones(N)
	 	negs(rand(N) > 0.5 ) = -1;
	 	R = R.*negs;
	 	R = triu(R,1) + triu(R,1)';
	 	J = beta.*Adj.*R;
	
	end








%% h vector
if h_on == 1
	hfield = normrnd(J0,distr,[1,N]);
	hfield = hfield.*(double(rand(1,N)> sprs));
	
else
	hfield = zeros(1,N);
end



%% Sparsity and distribution from nicolas decimation_MS_log_reg_LocIsing.m
%  
%  % distribution with abs of couplings distributed between 0.5*theta and 1.5*theta
%         w = theta*(rand(1,Np)+0.5);  		Our J? no just a vector. No idea whtf he's doing
%         segno = 2*double(rand(1,Np)>.5)-1; 	randomly assigns some value
%         w = segno.*w;						assigns
%  % sparsity of couplings
%         w(randperm(Np,round(sparsity*Np)))=0;
%         
%         % normalize the couplings (ensures the field acting on y is of order 1)
%         d = Np - round(sparsity*Np);
%         if d ~= 0
%             w = w/sqrt(d);
%         end