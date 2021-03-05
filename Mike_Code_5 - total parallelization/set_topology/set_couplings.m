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
 J = beta*Adj.*(J0 + sigJ*((J + J')/sqrt(2)));
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
end

end