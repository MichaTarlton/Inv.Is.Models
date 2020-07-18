% 7/07/18
% MODEL SELECTION WITH DECIMATION FOR LOGISTIC
% REGRESSION WITH BINARY INPUTS
% 
%	#INPUT:  - matrix X whose size should be M X Nx (# of observation times
%         number of features;
%          - vector Y of the output unit whose size is M x 1;
%          - h optional parameter: h = 1 tells you that you should consider 1-body operator on
%
% #OUTPUT:
%          - maximum likelihood estimator of parameters for all models
%          - normalised likelihood for all Np+1 non degenerate models
%         evaluated at the maximum likelihood point
%          - penalisations (BIC,AIC,MDLl, MDLu,ENT,COUNTPATTERNS see paper).
%          The penalizations come in a
%          structure form as a second output and each structure element is
%          an array of length Np+1
%          - best model: within each model selection framework
%          we evaluate the model that maximize the posterior.
%          We use the max function for the purpose, but notice that there
%          can be multiple maxima. In this case the max function returns
%          the first occourrence which is the DENSER model since the
%          models are ordered from the most (full) to the least dense (0 param)
%          model.
%          - ibest: is the index of the best model such that with it you
%          can retrieve the MLE estimates of parameters, the likelihood and
%          the penalisation that applies to the best model
%
% #LIMITATIONS: - only binary features here, for real features see the
%               routine logistic_Model_Selection
%               - if you want to include higher order you just need to
%               supply the matrix X with all combinations of features you
%               want to check
%               - the models are all only non degenerate models ... the routine
%               does not take into account models with degenerate
%               parameters
%
% NOTICE:
% 1) this can be applied to logistic regression with small to large number of
% BINARY inputs, since we are not enumerating all possible models but
% instead considering a fraction of them by decimating and starting from
% the most dense model

% 2) BIC and AIC are already knwon penalisation terms; For the others see
% the paper here a brief explanation: MDLu represents in general an upper
% bound to the absolute complexity of the model and it equals pi^n/n.
% The upper bound is reached when empirical distribution is uniform in the
% inputs. See the paper for the other complexity terms.

% 3) if you want to estimate logistic regression with many output variables
% which is a strucure in which your (independent) input units are connected
% to two or more output units, then what you need to do is to run this code
% for each output units independently and then combine the results. In fact
% the total logistic regression cost function (total likelihood) is the sum
% of the two. 

% 4) Importatly, if you have n features and want to run logistic model
% selection with all possible operators (single-body, two-body, three body
% and so on) you simply need to add more features. For instance if you have
% two features x1 and x2 and want to run only for pairwise models then just
% provide the matrix with dimensions 2xM where M is the number of
% observations. If you want to consider all possible operator then you need
% to put h = 1 for the one-body operator (or add one column with all ones)
% and another input node with values equal to x1*x2 for the three body operator.
% This is possible for all number of features but of course one should take
% into account that by increasing the number of features model selection will become slower and slower.
% Anyway the point is that this can just be fed into this routine simply as
% further features.

% 5) therefore in theory with this code you can run decimation
% model selection with jeffreys prior for all NON_DEGENERATE models of
% logistic regression with binary inputs.



function [w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,field)

% field is optional (put it only if you haven't already added a column with
% all ones)
try h = field; catch % this just makes sure there is a field entry t begin with. If none was entered, it just sets h = 0
    h = 0;
end

% check if features are binary
if ~all(X(:)==-1 | X(:)==1)
    error('Input and output features must be binary -1/1')
end

% Number of observations
M = size(X,1);

% check if the 1-body operator is on
if h == 1
    X = [ones(M,1),X];
end

% Number of parameters
Np = size(X,2);

% Binary representation of the model: this is a number of 2^Np bits and each
% bits tells you if the correspondent operator (parameter) is on or off
model = ones(Np+1,Np);

% Optimization algorithm: 'trust-region' with Hessian
options = optimoptions(@fminunc,'Algorithm','trust-region','GradObj','on','Hessian','on','Display','off');
%options = optimoptions(@fminunc,'Algorithm','trust-region','GradObj','on','Hessian','off','Display','off');

% Initialization
w_ML = zeros(Np+1,Np);
l_ML = zeros(Np+1,1);
Entropy = zeros(Np+1,1);
UniquePat = zeros(Np+1,1); %patterns equals within parity operation are considered equal
np = zeros(Np+1,1);

% Cycle over all Np+1 models
%parfor i_model = 1:Np
for i_model = 1:Np   
    % model index
    %disp(model(i_model,:))
    
    % see which parameters are 'on' in the current model
    i = logical(model(i_model,:));
    
    % store number of parameters
    np(i_model) = sum(i);
    
    % check the rank of the models
    if rank([X(:,i),Y]) < np(i_model) + 1
        disp('sgamato ... see this point')
    end
        
    % estimate the parameters with normal maximum likelihood
    f = @(x) logistic_likelihood(x,X(:,i),Y);
    x0 = 0.01*randn(1,sum(i));
    [w_,l_] = fminunc(f,x0,options);
    
    % MLE, likelihood and relative complexity at MLE
    w_ML(i_model,i) = w_';
    l_ML(i_model) = -l_;
    
    % calculate entropy and unique patterns (within parity operation)
    check = true; % you can put it to false once ensured it is ok
    [m,Hs] = countingPatterns(X(:,i),check);
    Entropy(i_model) = Hs;
    UniquePat(i_model) = m;
    
    % remove one spin and calculate next model
    [~,iabs] = min(abs(w_));
    i_rm = w_ML(i_model,:) == w_(iabs);
    model(i_model+1,:) = model(i_model,:);
    model(i_model+1,i_rm) = 0;
    
    %disp(['Model ', num2str(i_model)])
end

% zero parameters model
l_ML(i_model+1) = -log(2);

% save entropy and unique patt for all models
cost.entropy = Entropy;
cost.patterns = UniquePat;

%% MODEL SELECTION WITH SEVERAL CRITERIA (SEE PAPER)

% BIC
BIC = 0.5*np*log(M/(2*pi));
cost.BIC = BIC/M;
posterior.BIC = l_ML-cost.BIC;
[~,imax_BIC] = max(posterior.BIC);
BestModel.BIC = model(imax_BIC,:);

% AIC
AIC = np;
cost.AIC = AIC/M;
posterior.AIC = l_ML-cost.AIC;
[~,imax_AIC] = max(posterior.AIC);
BestModel.AIC = model(imax_AIC,:);

% MDL lower + BIC regularized
% regularized such that the penalisation is always increasing with n; 
% unreasonable results if it decreases with n
alpha = 2*exp(1)/pi;
mdl_l = 0.5*np.*log(0.5*alpha*pi*M./np); % using MDL lower bound + BIC
mdl_l(Np+1) = 0;
cost.MDLl = mdl_l/M;
posterior.MDLl = l_ML-cost.MDLl;
[~,imax_MDLl] = max(posterior.MDLl);
BestModel.MDLl = model(imax_MDLl,:);

% MDL upper + BIC
mdl_u = 0.5*np.*log(0.5*pi*M) - log(np); % using MDL upper bound + BIC
mdl_u(Np+1) = 0;
cost.MDLu = mdl_u/M;
posterior.MDLu = l_ML- cost.MDLu;
[~,imax_MDLu] = max(posterior.MDLu);
BestModel.MDLu = model(imax_MDLu,:);

% MDL entropy / formula proposed in the paper
Entropy0 = Entropy(1);
mdl_ent = 0.5*np.*log(exp(1)*M*Entropy./(np*Entropy0)) - log(np);
mdl_ent(Np+1) = 0;
cost.MDLent = mdl_ent/M;
posterior.MDLent = l_ML-cost.MDLent;
[~,imax_MDLent] = max(posterior.MDLent);
BestModel.MDLent = model(imax_MDLent,:);

% MDL count + BIC (same as before... count unique patterns another measure
% as entropy.. should be very close)
UEntropy = log2(UniquePat) + 1;
UEntropy0 = UEntropy(1);
mdl_count = 0.5*np.*log(exp(1)*M*UEntropy./(np*UEntropy0)) - log(np);
mdl_count(Np+1) = 0;
cost.MDLcount = mdl_count/M;
posterior.MDLcount = l_ML-cost.MDLcount;
[~,imax_MDLcount] = max(posterior.MDLcount);
BestModel.MDLcount = model(imax_MDLcount,:);

% index of the best model
IMAX.BIC = imax_BIC;
IMAX.AIC = imax_AIC;
IMAX.MDLl = imax_MDLl;
IMAX.MDLu = imax_MDLu;
IMAX.MDLent = imax_MDLent;
IMAX.MDLcount = imax_MDLcount;

end