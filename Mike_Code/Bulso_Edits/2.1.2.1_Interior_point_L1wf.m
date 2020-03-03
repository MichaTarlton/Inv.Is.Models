% INTERIOR POINT ALGORITHM FOR LOGISTIC REGRESSION WITH L1 without fields
% (see Koh, Kim and Boyd 2007)

% Input: training examples: x are N independent binary variables and b is
% the binary output. We have M training samples meaning b = vector of M
% components while x is a matrix of M*N components. In addition you have
% f_l which determines the choice of the regularizer. If f_l1 negative I use
% the regulariser according to Ravikumar l = abs(f_l1)*log(Nspins)/M; if l 
% positive < 2after calculating lambda max, you use f_l1 to calculate the 
% regulariser as follows l = f_l1*l_max; if f_l1 >=10 then you use the value 
% l = f_l1/1000 (the last is useful for selecting whatever value of lambda
% you want greater than l = 10/1000 = 0.01)_last update 23/01/2017;
% Output: parameters: teta(1) is the field, teta(2:N+1) are the N couplings
% and the lasting parameters are auxilliary variables used by the algorithm
% to enforce constraints (which are not returned as output);
% in addition you get adj which is the weights thresholding according to
% the Interior point method (see paper Koh et al..) and l_max according to
% the paper

% The code follows quite strictly Koh, Kim and Boyd (2007) with the
% following comments/exceptions
% 0) Remember that the parameters J corresponds to 2beta*J, I have not used
%    in this code the fields parameter h (see Interior_point_L1.m for
%    fields implementation. Hence I divide by 2 the parameters such that 
%    the output is beta*J
% 1) notice that if you change the initial condition, it does not converge;
% 2) I have changed the value of mu w.r.t. the paper because it won't
%    converge with mu = 2 so I am using mu = 1.8
% 3) I have symmetrized the hessian to ensure simmetry
% 4) I have implemented two possible set of regularisers (see the code) one
%    according with Ravikumar and the one with Kim, Koh et al.

function [teta,adj,l_max] = Interior_point_L1wf(x,b,f_l1,view)

%fileID = fopen('report_L1.txt','w');
fileID = 1;
% number of parameters Np
[M,N] = size(x);
%initialisation in case the program does not produce adj
adj = zeros(1,N);

% initial guess (according to Kim, Koh and Boyd)
teta(1:N) = 0;
teta(N+1:2*N) = 1;
if view == 1
  g=sprintf('%12.4e', teta);
  fprintf(fileID,'Initial condition teta_0 = %s\n',g);
end
% other possibles initial guess_just to see the difference
% --- notice that it DOES NOT CONVERGE (you can try and see) for other choices of the
% --- initial guess !!!! (see what happen in the case employing dual func)
%teta = (2*rand(1,2*N)-1);
%teta = 0.2*ones(1,2*N+1);

% lambda max (i.e. value of the regularizer for which you obtain all zeros)
J = teta(1:N);
z = b'.*(J*x');
l_max = max(abs((1/M)*(b'.*((1+exp(z)).^-1))*x));

% choice of the regularizer( log(Nspins)/M suggested by Revikumar, the
% other is the one adopted by Koh, Kim et al. in their paper
if f_l1 < 0 
    l = abs(f_l1)*sqrt(log(N+1)/M);
elseif f_l1 >= 0 && f_l1 < 2
    l = f_l1*l_max;
elseif f_l1 >= 10
    l = f_l1/1000;
end

% settings
t = 1/l;
eta = 1;
eps = 10^-8;
alfa = 0.01;
beta = 1/2;
it = 1;
itmax = 80;

while it < itmax && eta > eps
 % calculate function, gradient and hessian
 func = f(teta,x,b,t,l);
 g = gradient(teta,x,b,t,l);
 [Hr,D1,D2] = hessian_red(teta,x,b,t);
 % cholesky decomposition and step direction
 [L,p] = chol(Hr);
 if p ~= 0 
     %COMMENTED fprintf(fileID,'Error with L1: not positive definite Hessian\n');
     %COMMENTED fprintf(fileID,'det(Hessian) = %6.4f\n',det(Hr));
     if sum(sum(isinf(Hr))) == 0 && sum(sum(isnan(Hr))) == 0
     %COMMENTED fprintf(fileID,'minEig = %6.4f\n',min(eig(Hr)));
     end
     %COMMENTED fprintf(fileID,'b_avg = %6.4f\n',sum(b)/length(b));
     
  % thresholding condition according to Kim/Koh/Boyd
  J = teta(1:N);
  %COMMENTED fprintf(fileID,'mean_J +/- std_J = %6.4f +/- %6.4f\n',mean(J/2),std(J/2));
  z = b'.*(J*x');
  adj = abs((1/M)*(b'.*((1+exp(z)).^-1))*x);
  adj = adj > 0.9999*l;
  % divide by 2 since the parameters are 2betaJ and keep just
  % betaJ for comparison
  teta = teta/2;
  % removing variables used for regression before having checked their
  % expected convergence behaviour (see output of fprintf)
  teta(N+1:2*N) = [];
  return
 end
 gr(1:N) = g(1:N) - g(N+1:2*N)*D2*(D1^-1);
 step = zeros(2*N,1);
 step(1:N) = L\(L'\-gr');
 step(N+1:2*N) = -(D1^-1)*(g(N+1:2*N)' + D2*step(1:N));

 % backtracking using Armadijo condition and calculation of step size
 k = 0;
 m = g*step;
 while f(teta + (beta^k)*step',x,b,t,l) > func + alfa*(beta^k)*m
     k = k + 1;
 end
 % step calculation and updating parameters
 Dteta = (beta^k)*step';
 teta0 = teta;
 teta = teta + Dteta;
 % calculation of delta and the abs value of the step in parameters'space
 delta = abs(f(teta,x,b,t,l)-func);
 DTETA = sqrt(sum(Dteta.^2));

 % construction of the dual feasible point
 w = teta(1:N);
 C = min([M*l/max(abs((b'.*((1+exp(b'.*(w*x'))).^-1))*x)),1]);
 tetabar = (C/M)*((1+exp(b'.*(w*x'))).^-1);

 % evaluation of the dual gap 
 l_avg = sum(log(1+exp(-b'.*(w*x'))))/M;
 fstar = zeros(1,length(tetabar));
 for i = 1:length(tetabar)
 if tetabar(i) == 0 || tetabar(i) == 1/M
     fstar(i) = 0;
 elseif tetabar(i) < 0 || tetabar(i) > 1/M
     fstar(i) = Inf;
 elseif tetabar(i) > 0 && tetabar(i) < 1/M
     fstar(i) = (M*tetabar(i))*log(M*tetabar(i)) + (1-M*tetabar(i))*log(1-M*tetabar(i));
 end
 end
 fstar = sum(fstar)/M;
 eta = l_avg + fstar + l*sum(abs(w));
 
 % updating the value of the parameter t
 if k <= 1
     mu = 1.8;
     t_hat = 2*N/eta;
     t = max([mu*min([t_hat,t]),t]);
 end
 
 % displaying iterations and updating values
 if view == 1
   fprintf(fileID,'------------------------------------------------\n');
   fprintf(fileID,'Iteration: %d\n',it);
   fprintf(fileID,'|teta| = %10.4e\n', sqrt(sum(teta0.^2)));
   fprintf(fileID,'|Dteta| %10.4e\n', DTETA);
   fprintf(fileID,'eta %10.4e\n', eta);
   fprintf(fileID,'|abs(w)-u|_2 %10.4e\n', sqrt(sum((abs(teta(1:N))-teta(N+1:2*N)).^2)));
 end
 
 % updating and checking the number of iteration
 it = it + 1;
  if it == itmax
     fprintf(fileID,'WARNING: maximum number of iteration reached !!! \n');
  end
end

if view == 1
fprintf(fileID,'-------------------- FINISH --------------------\n');
if fileID ~= 1
fclose(fileID);
end
end

 % removing variables used for regression before having checked their
 % expected convergence behaviour (see output of fprintf)
 teta(N+1:2*N) = [];
  % thresholding condition according to Kim/Koh/Boyd
 J = teta(1:N);
 z = b'.*(J*x');
 adj = abs((1/M)*(b'.*((1+exp(z)).^-1))*x);
 adj = adj > 0.9999*l;
 % divide by 2 since the parameters are 2betaJ and keep just
 % betaJ for comparison
 teta = teta/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO MINIMIZE PLUS GRADIENT AND HESSIAN DEFINED HERE 
% simpler to generalize, you just have to change these functions and the
% arguments that you need (passed as arguments in the main routine and in
% the subroutine and maybe create or delete additional subroutine but the
% bulk of the code does not need to be changed (in summary only little
% modifications needed for generalizing)

function y = f(teta,x,b,t,l)
%get the number of samples and regression variables
[M,N] = size(x);
%divide the parameter array in fields and couplings
J = teta(1:N);
u = teta(N+1:2*N);
%calculation of the augmented weighted objective function (see Kim, Koh and
%Boyd 2007 pag.1529)
z = b'.*(J*x');
l_avg = sum(log(1+exp(-z)))/M;
Phi = -sum(log(u.^2-J.^2));
y = t*l_avg + t*l*sum(u) + Phi;
end

function dy = gradient(teta,x,b,t,l)
%get the number of samples and regression variables
[M,N] = size(x);
%divide the parameter array in fields and couplings
J = teta(1:N);
u = teta(N+1:2*N);
%calculation of the gredient of the augmented weighted objective function 
%(see Kim, Koh and Boyd 2007 pag.1531)
z = b'.*(J*x');
dy(1:N) = -(t/M)*(b'.*((1+exp(z)).^-1))*x + 2*J./(u.^2-J.^2);
dy(N+1:2*N) = t*l - 2*u./(u.^2-J.^2);
end

function [ddy_red,D1,D2] = hessian_red(teta,x,b,t)
%get the number of samples and regression variables
[M,N] = size(x);
%divide the parameter array in fields and couplings
J = teta(1:N);
u = teta(N+1:2*N);
%calculation of the hessian of the augmented weighted objective function 
%(see Kim, Koh and Boyd 2007 pag.1531)
z = b'.*(J*x');
A = (b*ones(1,N)).*x;
D0 = diag(exp(-z)./((1+exp(-z)).^2))/M;
D1 = diag(2*(u.^2+J.^2)./((u.^2-J.^2).^2));
D2 = diag(-4*(u.*J)./((u.^2-J.^2).^2));
D3 = D1 - D2*(D1^-1)*D2;

ddy_red = t*A'*D0*A+D3;

% ensure symmetry over random error and help matlab calculations
ddy_red = (ddy_red + ddy_red')/2;       
end

% function ddy = hessian(teta,x,b,t)
% %get the number of samples and regression variables
% [M,N] = size(x);
% %divide the parameter array in fields and couplings
% J = teta(1:N);
% u = teta(N+1:2*N);
% %calculation of the hessian of the augmented weighted objective function 
% %(see Kim, Koh and Boyd 2007 pag.1531)
% z = b'.*(J*x');
% A = (b*ones(1,N)).*x;
% D0 = diag(exp(-z)./((1+exp(-z)).^2))/M;
% D1 = diag(2*(u.^2+J.^2)./((u.^2-J.^2).^2));
% D2 = diag(-4*(u.*J)./((u.^2-J.^2).^2));
% 
% ddy = [t*A'*D0*A, D2;
%         D2,       D1];
% 
% % ensure symmetry over random error and help matlab calculations
% ddy = (ddy + ddy')/2;
% 
% end