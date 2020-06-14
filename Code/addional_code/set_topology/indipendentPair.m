%%%% FUNCTION GENERATING A INDIPENDENT PAIR MADE UP NETWORK

% bug discover: if N is odd than one node will be connected to two 
% different nodes (not a big deal), whereas if N is even all pairs 
% are made up of different nodes

function [A] = indipendentPair(N)

%N = 10; % number of spin in the network

A = zeros(N); % adjaciency matrix

b = 1:N;
imax = length(b);
idx = zeros(1,2);
v = zeros(1,N);

while imax ~= 0

    for j=1:2
    idx(j) = randi([1 imax]);
    v(j) = b(idx(j));
    b(idx(j)) = [];
    imax = length(b);
      if imax == 0
        break
      end
        
    end

    A(v(1),v(2)) = 1;

end    

A = A + A';

end