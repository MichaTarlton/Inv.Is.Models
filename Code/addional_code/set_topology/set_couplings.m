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