% DATE 28.06.2017
% function for cross validation for selecting the regularizer lambda for l1
% logistic regression
% DESCRIPTION: the function takes data and a range of lambda and perform K
% fold cross validation. Divide the dataset in K subsets. For each subset k
% = 1:K infer the model in the dataset in which the k-th subset has been removed.
% In the same dataset infer also the parameters with maximum likelihood. (Notice
% that I could not use the maximum likelihood from the Interior Point since there
% it is used an augmented cost function. So I inferred the parameters on the 
% retrieved model. Then I test it on the k-th sample by evaluating the maximum
% likelihood there. This measures the generalisation ability of the model.
% Then average over all K subsets. In this way for each value of lambda you 
% have a function and the minimum is the optimal lambda. In the debbugging phase
% you can provide the true parameters and compare the error made by the optimal 
% lambda (found with Kfold cross validation) with that made by the lambda that 
% minimize the error. I run some tests and it seems to do a good job, at
% least when M is not too small.
% INPUT
% - X = inputs values [M,Np] where M sample size and Np number of
%   parameters
% - Y = output values [M,1] where M sample size 
% - vlambda = vector of lambdas
% - K value of the K fold cross validation
% - w_true (optional value)for debugging
% OUTPUT
% - optimal lambda
% - curve showing the mean of log-likelihood function on the testing
%   samples for all choices of lambda. the minimum is optimum lambda
% - optional for debugging

function [model,lambda_opt,CV,err] = KfoldCV(X,Y,vlambda,K,true_param)

door = 1;
try w_true = true_param; catch
    door = 0;
end

% !!! for debugging
if door == 1
err = zeros(1,numel(vlambda));
end

% M = sample size and Np = number of parameters
[M,Np] = size(X);

% initialisation
CV = zeros(1,numel(vlambda));

% check that K is a divisor of M
if mod(M,K)~=0
    error('K is not a divisor for M');
end

% Optimization algorithm: 'trust-region' with Hessian
options = optimoptions(@fminunc,'Algorithm','trust-region','GradObj','on','Hessian','on','Display','off');

% loop over all lambda
for ilam = 1:numel(vlambda)
   lambda = vlambda(ilam); 
   CV(ilam) = 0;
   
   % loop over subsets
   for k = 1:K 
       
      % select the k-th subset, training and testing part of the dataset
      sel_k = (k-1)*round(M/K)+1:k*round(M/K); others = setdiff(1:M,sel_k);  
      X_train = X(others,:); Y_train = Y(others);
      X_test = X(sel_k,:); Y_test = Y(sel_k);
       
      % model selection with l1 regularisation (Interior point routine)
      view = 0;
      [w_l1,modelL1Kim,~] = Interior_point_L1wf(X_train,Y_train,lambda,view);
      warning('off','MATLAB:singularMatrix'); % especially with small samples these warnings pop up a lot
      warning('off','MATLAB:nearlySingularMatrix');
      modelL1 = double(modelL1Kim | abs(w_l1)>= 1e-3);  % see COMMENT ON THE THRESHOLDING CONDITION below
      %modelL1 = double(modelL1Kim);
      
      % fit parameter to Xtrain and Ytrain with constraint of modelL1
      i_active = find(logical(modelL1));
      if sum(i_active)>0
        Xmodel = X_train(:,i_active); np = size(Xmodel,2);
        
        f = @(x) logistic_likelihood(x,Xmodel,Y_train);
        x0 = zeros(1,np);
        x = fminunc(f,x0,options);
        w_inf = zeros(1,Np); w_inf(i_active) = x;
        
      else
        w_inf = zeros(1,Np);
      end
      
      % measures minus normalised likelihood in [Xtest,Ytest] with modelL1 and param. inferred
      l_test = logistic_likelihood(w_inf,X_test,Y_test);
      CV(ilam) = CV(ilam) + l_test/K; 
   end
    
   % !!! for debugging 
   % check that optimal lambda minimize the error with respect to the 
   % true model (this should be true at least for large samples)
   % model selection with l1 regularisation on the whole sample
   if door == 1
       view = 0;
       [w_l1,modelL1Kim,~] = Interior_point_L1wf(X,Y,lambda,view);
       warning('off','MATLAB:singularMatrix'); % especially with small samples these warnings pop up a lot
       warning('off','MATLAB:nearlySingularMatrix');
       modelL1 = double(modelL1Kim | abs(w_l1)>= 1e-3);  % see COMMENT ON THE THRESHOLDING CONDITION below
       %modelL1 = double(modelL1Kim);
       true_model = abs(w_true)>0;
       err(ilam) = sum(modelL1~=true_model)/Np;
   end
    
end
      % optimal regularizer and model inferred with it
      [~,iopt] = min(CV);
      lambda_opt = vlambda(iopt);
      [w,modelKIM,~] = Interior_point_L1wf(X,Y,lambda_opt,view);
       warning('off','MATLAB:singularMatrix'); % especially with small samples these warnings pop up a lot
       warning('off','MATLAB:nearlySingularMatrix');
       model = double(modelKIM | abs(w)>= 1e-3);  % see COMMENT ON THE THRESHOLDING CONDITION below


      % !!! for debugging
      % check that lambda optimal is the one that approximatively minimize
      % the error with respect to the truth
      if door == 1
          figure;plot(vlambda,CV);
          
          figure;plot(vlambda,err);
          [~,imin] = min(err);
          hold on;line([vlambda(imin),vlambda(imin)],[0,1],'linestyle','-.','color','black','linewidth',2);
          hold on;line([vlambda(iopt),vlambda(iopt)],[0,1],'linestyle',':','color','black','linewidth',2);
      end

end


% COMMENT ON THE THRESHOLDING CONDITION:
% the THRESHOLDING CONDITION is the one used in Kim/Koh 2007 and it
% is based on the optimality conditions (8) (see pag.1534 of the
% paper). Sometimes does something odd when thresholding
% (expecially at small M). That's why I modified putting an
% additional condition which correct when it is missing couplings
% larger than 1e-3 in absolute value because usually the gap
% induced by the regularization is below 0.001 .. usually around
% 1e-4 or 1e-5 but this might depend on the problem at hand