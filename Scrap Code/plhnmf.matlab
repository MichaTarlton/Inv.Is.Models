%%%% PL-MH.m
% From Berg, Nguyen 2017
% eq.113

Jik = (1-mi.^2) Sum()


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
	Cni = Cij
	Cni(i,:) = []
	Cni(:,i) = []
	PLCstruct(i).Cni = Cni
 end


%% fine fdoing it myself

Cout = zeros(size(Cinj))

for ii = 1:size(Cinj,1)

Cni=PLCstruct(ii).Cni
Cjk = (Cni^-1)'

for jj = 1:size(Cinj,2)

Cout(ii,jj) = sum(Cinj(ii,jj).*Cjk(jj,:),'all')

end
end
