
%% Psudo_LLH_Int.m

function [l,Dl,DDl,J,h] = Pseudo_LLH(T,S) %|Old Log_LLH

% Doing it in the most literal caveman way possible
% Berg Eq.109 & 110

% Berg Eq. 109, summarize over the number of spike samples T 
% Normally I'd plug this into another sum over i (eq. 110)

% J and h to be found by maximization function
% si a single spike series, indexed by mu for 1 to T
% si = S_hat 

%Assuming S:(NxT)

syms J h % Not certain if I should make this a symbolic function here or not
		% Not. I need t spit this out into different parts and exterior and interior.

% Attempting to do something more like how Nicola did it where the sums are just baked into the array multiplication
% This doesn't work though due to how it effectively moves the sum(mu)
%or is it elementwise mltiplication on the (S^mu)_i
J = J - diag(diag(J));

LDi = (1/T)*sum(log(0.5*(1 + S(mu,:).*tanh(h + J*S(mu,:)))),2); %| why do I ahve it set to sum along the second dimension

LD = sum(LDi); %Could just as easily incorporate it into above

% Gradient of llh
%DLD = [(S(mu,:)^2 .* sech^2(S(mu,:) J + h))/(1 + S(mu,:) tanh(S(mu,:) J + h)), (S(mu,:) sech^2(S(mu,:) J + h))/(1 + S(mu,:) tanh(S(mu,:) J + h))] %|wolfram generated
DLD= [ -(S*(tanh(h + J*S)^2 - 1))/(2*T*((S*tanh(h + J*S))/2 + 1/2)), -(S^2*(tanh(h + J*S)^2 - 1))/(2*T*((S*tanh(h + J*S))/2 + 1/2))]

% Hessian | or should this just be the laplace
DDLD = [ (S * tanh(h + J*S)*(tanh(h + J*S)^2 - 1))/(T*((S*tanh(h + J*S))/2 + 1/2)) - (S^2*(tanh(h + J*S)^2 - 1)^2)/(4*T*((S*tanh(h + J*S))/2 + 1/2)^2)...
		,(S^2*tanh(h + J*S)*(tanh(h + J*S)^2 - 1))/(T*((S*tanh(h + J*S))/2 + 1/2)) - (S^3*(tanh(h + J*S)^2 - 1)^2)/(4*T*((S*tanh(h + J*S))/2 + 1/2)^2)...
		,(S^2*tanh(h + J*S)*(tanh(h + J*S)^2 - 1))/(T*((S*tanh(h + J*S))/2 + 1/2)) - (S^3*(tanh(h + J*S)^2 - 1)^2)/(4*T*((S*tanh(h + J*S))/2 + 1/2)^2)...
		,(S^3*tanh(h + J*S)*(tanh(h + J*S)^2 - 1))/(T*((S*tanh(h + J*S))/2 + 1/2)) - (S^4*(tanh(h + J*S)^2 - 1)^2)/(4*T*((S*tanh(h + J*S))/2 + 1/2)^2)]



% M = size(X,1);
% % if there are more parameters than spins than add one more spin at the
% % beginning which is always equal to one to take into account for external
% % field parameter
% if length(w) > size(X,2) 
% X = [ones(M,1),X];
% end
% 
% m = Y'*X/M;
% 
% % minus normalised log-likelihood function
% l = - m*w' + sum(log(2*cosh(w*X')))/M;
% 
% % gradient of the likelihood function
% Dl = - m + tanh(w*X')*X/M;
% 
% % hessian of the likelihood function
% DDl = X'*diag(cosh(w*X').^(-2))*X/M;

end