%% Forward ising

function Jstruct = ising(N,jn)
Jstruct = struct;

%% Rand Dist
%for i = 1:jn
%
%	J = randi([0 1], N,N);
%	J = 2*J-1;
%	J = J - diag(diag(J));
%	%% J2 = triu(J,1) %Upper triangle minus diagonal
%	Jstruct(i).Jrand = J
%end    	

%% Gaussian Dist (depreceated)
%this is flawed you will obviously just have half aboe .5 and half under
% You need something more 
%for i = 1:jn
%	R = double(normrnd(0.5,1/N,[N,N]) > 0.5);
%	R = 2*R-1;
%	R = R - diag(diag(R)); 
%	R2 = triu(R,1) %Upper triangle minus diagonal
%   R3 = R2 + triu(R2)';
%	Jstruct(i).Jgaus = R3
%end


%% Gaussian Dist
for i = 1:jn
	R = double(normrnd(0,1/N,[N,N]));
	R = R - diag(diag(R)); 
	R2 = triu(R,1) %Upper triangle minus diagonal
   	R3 = R2 + triu(R2)';
	Jstruct(i).Jgaus = R3
end