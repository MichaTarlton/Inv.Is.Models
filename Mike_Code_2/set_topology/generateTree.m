%%%% FUNCTION GENERATING A TREE WITH COORDINATION NUMBER z

% the function generates a slightly modified version of the cayley tree
% since I built it in such a way that each node has "z" children. As a
% results the first node (the generator) has z connections whereas all the
% others have "z+1" except the last layer which can have a variable number
% of neibours according to the number of nodes. Therefore the ONLY
% difference between this graph and a cayley tree is the connectivity of
% the first node which should be also "z+1" but instead it is "z" (the
% other difference is that, apart from the first node, a Cayley tree with "z" neibours corresponds to
% our graph with "z-1" neibours)

function [A] = generateTree(N,z)

%N = 100; % number of spin in the network
%z = 10; % coordination number

A = zeros(N); % adjaciency matrix

b = 1:N;
imax = length(b);
l = 1;
v = zeros(N);
lmax = Inf;

while imax ~= 0;
    
    for j = 1:z^(l-1);
        idx(l,j) = randi([1 imax]);
        v(l,j) = b(idx(l,j));
        b(idx(l,j)) = [];
        imax = length(b);
        %display([l j imax]);
        if imax == 0
            lmax = l;
            break
        end
    end
    
    if l > 1
        for k = 1:z^(l-2);
            for n = z*(k-1)+1:z*k;
                if l == lmax && n > j
                    break
                end
                A(v(l-1,k),v(l,n))=1;
            end
        end
    end
    
    l = l + 1; 
    
end    

A = A + A';

end