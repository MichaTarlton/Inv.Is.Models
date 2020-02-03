% DECIMATION TESTS -  7/07/18
% Using the function decimation_logistic_Model_Selection to test performances of
% recovering structure in logistic regressions models with binary input
% (two layers neural network with one ouput unit) when large input size.
% Adapted from the same simulations on small system.


%% TESTING MODEL SELECTION APPROACH TO STRUCTURE LEARNING IN LOGISTIC REGRESSION -------------------------

function testing_MS_log_reg_decimation(trial)

% ---- # PARAMETER INITIALIZATION

% type of experiment
tipo = 1; % ising, calls function_MS_log_reg_ultimate
% tipo = 2; % uniform and localiz distr, calls function_MS_log_reg_ultimate2

% number of inputs features
Nx_v = 50;

% field is active? (yes if h = 1, no h = 0)
h = 0; %Nicola said that this is meant to stay off, there currently is not an implementation for standing fields

% number of observations: notice that K (for CV with l1) has to be an exact divisor for M (usually I take K = 5)
alpha = [5,50,200];
M = alpha*Nx_v; % only if Nx_v is scalar

% tuning parameter for localisation of distribution of input features
if tipo == 1
   vbeta = [0.01,1,1.5]; % when using "decimation_MS_log_reg_LocIsing" otherwise commented
end

% coupling strenght: mean is vtheta +- vtheta/2
theta = 1;%0.5;

% choose to compare with L1
comparewithL1 = true;

% loops on independent trials for collecting statistics
sumtime = 0;

% loop on input size
for iNx = 1:numel(Nx_v)
    Nx = Nx_v(iNx);
    
    % function for testing model selection on logistic regression models
    if tipo == 1
        time = decimation_MS_log_reg_LocIsing(Nx,h,M,theta,vbeta,trial,comparewithL1);
    elseif tipo == 2
        % I have not written the code for this
        time = decimation_MS_log_reg_LocUniform(Nx,h,M,theta,trial);
    end
    
end

% elapsed time for one trial
sumtime = sumtime + time;
fprintf('Elapsed time %4.2f seconds\n',sumtime);