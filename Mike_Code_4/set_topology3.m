function Adj = set_topology(topology,N,c)

if topology == 1
%---cayley tree with coordination number c
[Adj] = generateTree(N,c);
elseif topology == 2
%---fully connected topology
Adj = ones(N) - eye(N);
elseif topology == 3
%---indipendent pair topology
[Adj] = indipendentPair(N);    
elseif topology == 4
%---2D Ising lattice
 [Adj] = lattice(N,2);
elseif topology == 5
%---Erdos Reyni random graph
p = c/(N-1); %probability of two spins being connected
Adj = rand(N,N) < p;
Adj = triu(Adj,1);
Adj = Adj + Adj';
%  comment1: in the limit of Nspin tend to infty this graph gives the same
%  result as fully connected with dilution 1-p
%  comment2: the average vertex degree is c = p*(Nspin-1), for instance
%  if Nspin = 50, in order to have c = 3 and 4 you need p = 0.0612 and p = 0.0816 
%  if Nspin = 20 you can use p = 0.2 for c = 4 and p = 0.15 for c = 3; p =
%  0.05 for c = 1
%  if Nspin = 41, then p = 3/40 for c = 3 and p = 1/10 for c = 4
elseif topology == 6
%---Erdos Reyni random graph from Pachenko 2014
p = c/(N-1); %probability of two spins being connected
Adj = rand(N,N) < p;
Adj = triu(Adj,1);
Adj = Adj + Adj';
%  comment1: in the limit of Nspin tend to infty this graph gives the same
%  result as fully connected with dilution 1-p
%  comment2: the average vertex degree is c = p*(Nspin-1), for instance
%  if Nspin = 50, in order to have c = 3 and 4 you need p = 0.0612 and p = 0.0816 
%  if Nspin = 20 you can use p = 0.2 for c = 4 and p = 0.15 for c = 3; p =
%  0.05 for c = 1
%  if Nspin = 41, then p = 3/40 for c = 3 and p = 1/10 for c = 4
elseif topology == 7
%--- Star network
[Adj] = star(N); 
elseif topology == 8
%--- Small World - Watts Strogatz Graph
p = c/(N-1); %Probability that an edge will switch to a random node, \beta in the Watts Strogatz model. 
			% I'm setting the probability to be the same as the probability for the average vertex degree
			% So it will scale with N and c
			% When p=0 then we have a ring graph, when p=1 it is a random graph
graphy = WattsStrogatz(N,c,p); %this is passedout in matlab graph format. I need to convert to array
Adj = full(adjacency(graphy)); % Turns it into an adjaceny matrix (symmetric) and converts to array form
end

end