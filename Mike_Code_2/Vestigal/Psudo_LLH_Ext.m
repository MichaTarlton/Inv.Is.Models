%%Psudo_LLH_Ext.m
% Part of the maximization funtion from Nguye al 17
% made to interface directly with 

function [l,Dl,DDl,Jnew,hnew] = Pseudo_LLH_Ext(T,S,J,h) %

% make inputs J and h optional then add feature here that creates random gen 
% 
%N = size(S,2)
%
%if J == dne
% J = ones(N,N)
%end
%x
%if h == dne
% h = ones(1,N)
%end

% syms J h %|so these don't actually go anywhere

% Maybe feed in PLLH_MF J and h estimates as starter assumptions

options = optimoptions(@fminunc,'Algorithm','trust-region','GradObj','on','Hessian','on','Display','off'); %copy of options taken from bulso


 		% f = @(x) logistic_likelihood(x,X(:,i),Y);
 		% x0 = 0.01*randn(1,sum(i));
 		% [w_,l_] = fminunc(f,x0,options);

f = @(J,h) Pseudo_LLH_Int(T,S,J,h)
[l,Dl,DDl,J,h] = fminunc(f,J,h,options);