%%%% FUNCTION GENERATING A STAR NETWORK

function [A] = star(N)

A = zeros(N);
i = randi([1 N]);
A(i,:) = 1;
A(i,i) = 0;
A = A + A';

end