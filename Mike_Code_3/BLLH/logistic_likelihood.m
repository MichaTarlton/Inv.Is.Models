% MINUS NORMALIZED LOG-LIKELIHOOD FUNCTION FOR LOGISTIC REGRESSION
% # INPUT:
% 1) weight vector: w(1) = h; w(2:end) = J --> size(w) = [1,N] or [1,N+1]
% 2) X feature input: size(X) = [M,N]
% 3) Y target variable: size(Y) = [M,1]
% # OUTPUT:
% 1) likelihood function at w
% 2) gradient of the likelihood at w
% 3) Hessian matrix at w

function [l,Dl,DDl] = logistic_likelihood(w,X,Y)

M = size(X,1);
% if there are more parameters than spins than add one more spin at the
% beginning which is always equal to one to take into account for external
% field parameter
if length(w) > size(X,2) 
X = [ones(M,1),X];
end

%if M >= 1e5
	%X2 = sparse(X); % fucking around with sparse matricies to get around memory limitations, chuck out as soon as things begin to fail
	%Y2 = sparse(Y); % fucking around with sparse matricies to get around memory limitations, chuck out as soon as things begin to fail
%	% appears to work leaving for now, no not especially. Actually seems very slow.
%end

m = Y'*X/M;

% minus normalised log-likelihood function
l = - m*w' + sum(log(2*cosh(w*X')))/M;

% gradient of the likelihood function
Dl = - m + tanh(w*X')*X/M;

% hessian of the likelihood function
% causes massive array sizes to be stored, do not use unless for some reason necessary
% DDl = X'*diag(cosh(w*X').^(-2))*X/M;

end