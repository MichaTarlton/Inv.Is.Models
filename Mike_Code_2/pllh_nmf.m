
First we need to creat the submatricies for Cij, probably need to port these to main code as well

so for Cinj 

	Cinj = Cij
	for n = 1:i

%% Removing diagonal, no tdoing this after all

n =length(X)
Y = X.';
reshape(Y(eye(n)'~=1),[],n).'

%% Instead just setting diag to zero
% hmm didn't work bc indexing over

Cinj = Cij - diag(diag(Cij));


%% Making C_\i
 

n = length(Cij) 

for i =1:n 
Cni = Cij; 

%Cni(i,:) = 0 %| originally = [] 
%Cni(:,i) = 0 %| 
Cni(:,i) = [];
PLCstruct(i).Cnj = Cni;
Cni(i,:) = [];
PLCstruct(i).Cnij = Cni;
PLCstruct(i).Cninv = Cni^-1;
end

for j = 1:n

Cplline = PLCstruct(j).Cnj(j,:)*PLCstruct(j).Cninv;
PLCstruct(j).Cplline = Cplline;
Jpre(j, :) = Cplline;

end

pij = [1-mi.^2]'
PLCstruct(1).Jpre = pij.*Jpre; %| Storing Jpre and adding magnetizationg thingy


Jpl = [tril(Jpre,-1) zeros(N, 1)] + [zeros(N,1) triu(Jpre)]; %| This adds the diagonal of 0 

%for k = 1 : n
%  Jpre(k, :) = PLCstruct(k).Cplline;
%end

% Pij = diag(1-mi.^2);









% Vestigial
%%% fine fdoing it myself
%
%Cout = zeros(size(Cinj))
%
%for ii = 1:size(Cinj,1)
%
%Cni=PLCstruct(ii).Cni
%Cjk = (Cni^-1)'
%
%for jj = 1:size(Cinj,2)
%
%Cout(ii,jj) = sum(Cinj(ii,jj).*Cjk(jj,:),'all')

end
end
