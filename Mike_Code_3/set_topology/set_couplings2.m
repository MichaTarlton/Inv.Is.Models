% Inputs 
% Couplings: he said to do the delta func so set to 2
% J0: the mean, set to one, not 
% sigJ: std dev | I used the SK model deviation: beta./sqrt(N)
%		This can't be right tho.
% Adj: topology as given by set_topolgy.m
%


function J = set_couplings(couplings,beta,J0,sigJ,Adj)
N = length(Adj);

%---gaussian
if couplings == 1
 J = randn(N);
 J = beta*Adj.*(J0 + sigJ*((J + J')/sqrt(2))); % none of these go to negatives, causing trouble in the S generation
	% my std dev isn't any good here, it tops above 1 easily. Setting J0=0.5 makes it better but consider using the nthroot(N,3)
 
 %% Adding the neative values
 R = ones(N);
 R(rand(N) > 0.5 ) = -1;
 R = triu(R,1) + triu(R,1)'; % to add some negative values

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

%---delta function
elseif couplings == 2 
 J = beta*J0*Adj;
%---double delta function
elseif couplings == 3
 J = randn(N);
 J = beta*J0*Adj.*sign(J + J'); 
% when use delta functions make sure that the variance of the prior
% contains the value of beta (i.e. if sig(3) is 1 is not good to choose
% beta = 5, instead you could pick beta = 0.5

%| Mike's gaussian from the "SK" model gen: JHnorm
% difference is mine is meaned around 0
elseif couplings == 4
	%R = double(normrnd(0,beta./sqrt(N),[N,N])); %Originally
	R = double(normrnd(0,beta./nthroot(N,3),[N,N]));
	R = R - diag(diag(R));
    R2 = triu(R,1) + triu(R,1)';
    J = Adj.*R2;

elseif couplings == 5
 % My own take on Nicola's gaussian
 R = ones(N);
 R(rand(N) > 0.5 ) = -1;
 J = R.*beta*Adj.*normrnd(0.5,1/nthroot(N,3),[N,N]);
 J = triu(J,1) + triu(J,1)';

end

end