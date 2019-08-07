% CALL DECIMATION MODEL SELECTION FOR "ISING" INPUTS -  7/07/18
% function for testing model selection on logistic regression models with N
% large so that a  decimation algorithm is needed to walk through the
% relevant models and do model selection. The localisation of the input is
% imposed by requiring the input distribution being a ising distribution
% with inverse temperature 1/beta

% INPUT
% - Nx = number of inputs
% - h = if a field parameter should be considered
% - M = vector of sample sizes
% - vtheta = vector of strenghts of parameters. For each elements of the
%       vector vtheta(i) the parameters are drawn from an uniform distribution
%       picked around vtheta(i) and extending from vtheta(i)/2 to 3*vtheta(i)/2
% - vbeta = distribution of input featurs according to Ising with parameters
%       beta = vbeta(i)
% - number of the trial experiment

% OUTPUT
% the output is a .mat file for each value of the sparsity, strength of
% parameter (theta) and localisation of the distribution (beta)
% In the output you can find errors, false positive and false negative with
% respect to the true parameters generating the data for different
% criteria: AIC,BIC,MDL lower bound, MDL upper bound, MDL entropy and MDL
% counting (see paper and notes)



function time = decimation_MS_log_reg_LocIsing(Nx,h,M,theta,vbeta,trial,comparewithL1)

% if you wish to compare with L1
try cmpL1 = comparewithL1; catch %Try staement executes statement and then upon error rolls over to catch block
    cmpL1 = false;
end

% call rng for reproducibility
rng(trial)

tic %starts stopwatch

% number of parameters
Np = Nx + h; %number of inputs plus field parameter (0 or 1)

% fraction of active couplings
vNp_inactive = round(linspace(0,Nx,6));

for iNp_inactive = 1:numel(vNp_inactive)
    Np_inactive = vNp_inactive(iNp_inactive);
    sparsity = Np_inactive/Np;
    
    for vb = 1:length(vbeta)
        
        % tuning parameter localisation of inputs distribution inverse temperature Ising ferromagnetic
        beta = vbeta(vb);
        
        % input nodes distributed as Ferromagnetic Ising model with temperat 1/beta
        xsi = ones(1,Nx);
        W = xsi'*xsi; % hebbian learning to store the config in the connectivity
        W = W - diag(diag(W)); % no self interactions
        W = W/Nx; % so effective field is O(1) and comparable with beta
        epsil = 0.1;
        H = epsil*xsi; % weak field order 0.1 to break the simmetry
        
        % initializations
        W_ML = zeros(Np+1,Np,numel(M));
        MODELS = zeros(Np+1,Np,numel(M));
        LIKELIHOOD = zeros(Np+1,numel(M));
        ENTROPIA = zeros(Np+1,numel(M));
        NPATTERNS = zeros(Np+1,numel(M));
        NPARAMETERS = zeros(Np+1,numel(M));
        
        RANGO = zeros(numel(M),1);
        Fp.BIC = zeros(numel(M),1);Fn.BIC = zeros(numel(M),1);Err.BIC = zeros(numel(M),1);
        BestModel.BIC = zeros(numel(M),Np);CostBest.BIC = zeros(numel(M),1);
        Posterior.BIC = zeros(Np+1,numel(M));
        Fp.AIC = zeros(numel(M),1);Fn.AIC = zeros(numel(M),1);Err.AIC = zeros(numel(M),1);
        BestModel.AIC = zeros(numel(M),Np);CostBest.AIC = zeros(numel(M),1);
        Posterior.AIC = zeros(Np+1,numel(M));
        Fp.MDLl = zeros(numel(M),1);Fn.MDLl = zeros(numel(M),1);Err.MDLl = zeros(numel(M),1);
        BestModel.MDLl = zeros(numel(M),Np);CostBest.MDLl = zeros(numel(M),1);
        Posterior.MDLl = zeros(Np+1,numel(M));
        Fp.MDLu = zeros(numel(M),1);Fn.MDLu = zeros(numel(M),1);Err.MDLu = zeros(numel(M),1);
        BestModel.MDLu = zeros(numel(M),Np);CostBest.MDLu = zeros(numel(M),1);
        Posterior.MDLu = zeros(Np+1,numel(M));
        Fp.MDLent = zeros(numel(M),1);Fn.MDLent = zeros(numel(M),1);Err.MDLent = zeros(numel(M),1);
        BestModel.MDLent = zeros(numel(M),Np);CostBest.MDLent = zeros(numel(M),1);
        Fp.MDLcount = zeros(numel(M),1);Fn.MDLcount = zeros(numel(M),1);Err.MDLcount = zeros(numel(M),1);
        Posterior.MDLent = zeros(Np+1,numel(M));
        BestModel.MDLcount = zeros(numel(M),Np);CostBest.MDLcount = zeros(numel(M),1);
        Posterior.MDLcount = zeros(Np+1,numel(M));
        Fp.L1 = zeros(numel(M),1);Fn.L1 = zeros(numel(M),1);Err.L1 = zeros(numel(M),1);
        BestModel.L1 = zeros(numel(M),Np);LambdaBestL1 = zeros(numel(M),1);
        
        
        % distribution with abs of couplings distributed between 0.5*theta and 1.5*theta
        w = theta*(rand(1,Np)+0.5);
        segno = 2*double(rand(1,Np)>.5)-1;
        w = segno.*w;
        
        % sparsity of couplings
        w(randperm(Np,round(sparsity*Np)))=0;
        
        % normalize the couplings (ensures the field acting on y is of order 1)
        d = Np - round(sparsity*Np);
        if d ~= 0
            w = w/sqrt(d);
        end
        



        
        % ---- # GENERATION OF DATA
        % NOTICE: here X is not observed with frequency p, but drawn from a distribution p. 
        % APPROXIMATION IN THE SIMULATIONS!
        % generate input distribution according to Ising distribution by means of Monte Carlo
        % generate data in the input layer according to the desired localisation
        % generate more data to avoid many Monte Carlo
        rep = true;
        XX = Metropolis_Hasting(3*M(numel(M)),Nx,beta*W,beta*H,rep);
        X = ones(M(numel(M)),Nx); Y = ones(M(numel(M)),1);
        door = true;





        % the rank of the composite matrix has to be full rank
        while rank([X,Y]) ~= Np+1 || door
            door = false;
            X = XX(randperm(3*M(numel(M)),M(numel(M))),:);
            
            % if h =1, X will have the first column filled with ones
            if h == 1
                X = [ones(M(numel(M)),1),X];
            end
            
            % generate values in the output layer
            % a) evaluate P(y(i)=1|x) for each observation --> size py_x = (M,1)
            py_x = exp(X*w')./(2*cosh(X*w'));
            % b) draw M configurations of the ouput layer from Bernoulli distribution
            % according to the probability of success p(y|x) --> size Y = (M,Ny)
            Y = 2*double(rand(M(numel(M)),1) <= py_x)-1;
        end
        
        % loop on data matrix size M
        for m = 1:numel(M)
            
            % random selection of lenght M(m) such that [X,Y] is full rank
            while RANGO(m) ~= Np + 1
                
                % random selection of length M(m)
                sel = randperm(M(numel(M)),M(m));
                
                % check the rank of the selection X(sel,:)
                RANGO(m) = rank([X(sel,:),Y(sel,:)]);
                
                % write down when the matrix is not full rank
                if RANGO(m) ~= 0 && RANGO(m) ~= Np + 1
                    fprintf('---------- NOT FULL RANK SELECTION ------------\n');
                    fprintf(['#inputs:',num2str(Nx),' #field:',num2str(h),' #inactive:',num2str(Np_inactive),...
                        ' beta:',num2str(beta),' #data:',num2str(M(m)),' trial:',num2str(trial),...
                        ' RANK:',num2str(RANGO(m)),'\n']);
                    fprintf('---------- SELECTED A NEW ONE ------------\n');
                end
                
            end
            
            % ---- # MODEL SELECTION ON NESTED MODELS THROUGH DECIMATION WITH INFORMATION CRITERIA
            % notice 1) that if there are two bests models the routine pick the
            % sparser);
            % notice 2) if there is a field you don't need to put this
            % information in the routine since in that case X already contains
            % an additional column

            [w_ML,l_ML,posterior,Cost,Best,ibest] = decimation_logistic_Model_Selection(X(sel,:),Y(sel));
            fprintf('Model selection time %4.3f\n',toc)
            
            % ---- # L1 REGULARISATION WITH K FOLD CROSS VALIDATION
            % help for selecting the regulariser
            if cmpL1 && M(m) < 5500
                vlambda = 20:20:1000; % this corresponds to lambda 0.02:0.02:1, see Interior point code
                K = 5;
                [modelL1,lambda_opt] = KfoldCV(X(sel,:),Y(sel),vlambda,K);
                fprintf('l1 time %4.3f\n',toc)
            elseif cmpL1 && M(m) > 5500
                modelL1 = zeros(1,Np);
                lambda_opt = 0;
                fprintf('l1 not evaluated here: M too large and l1 algorithm would be too slow\n')
            elseif ~cmpL1
                modelL1 = zeros(1,Np);
                lambda_opt = 0;
                fprintf('l1 comparison not selected \n')
            end
            
            % counting false positive and false negative rates
            N1 = sum(double(abs(w)>0) == 1);
            N0 = sum(double(abs(w)>0) == 0);
            
            % model selection with decimation
            W_ML(:,:,m) = w_ML;
            MODELS(:,:,m) = logical(w_ML);
            LIKELIHOOD(:,m) = l_ML;
            ENTROPIA(:,m) = Cost.entropy;
            NPATTERNS(:,m) = Cost.patterns;
            NPARAMETERS(:,m) = squeeze(sum(MODELS(:,:,m),2));
            
            % best model, penalis for best model, false positive and negative
            BestModel.BIC(m,:) = Best.BIC; CostBest.BIC(m) = Cost.BIC(ibest.BIC);
            Posterior.BIC(:,m) = posterior.BIC;
            Fp.BIC(m) = sum(Best.BIC - double(abs(w)>0) == 1)/N0;
            Fn.BIC(m) = sum(Best.BIC - double(abs(w)>0) == -1)/N1;
            
            BestModel.AIC(m,:) = Best.AIC; CostBest.AIC(m) = Cost.AIC(ibest.AIC);
            Posterior.AIC(:,m) = posterior.AIC;
            Fp.AIC(m) = sum(Best.AIC - double(abs(w)>0) == 1)/N0;
            Fn.AIC(m) = sum(Best.AIC - double(abs(w)>0) == -1)/N1;
            
            BestModel.MDLl(m,:) = Best.MDLl; CostBest.MDLl(m) = Cost.MDLl(ibest.MDLl);
            Posterior.MDLl(:,m) = posterior.MDLl;
            Fp.MDLl(m) = sum(Best.MDLl - double(abs(w)>0) == 1)/N0;
            Fn.MDLl(m) = sum(Best.MDLl - double(abs(w)>0) == -1)/N1;
            
            BestModel.MDLu(m,:) = Best.MDLu; CostBest.MDLu(m) = Cost.MDLu(ibest.MDLu);
            Posterior.MDLu(:,m) = posterior.MDLu;
            Fp.MDLu(m) = sum(Best.MDLu - double(abs(w)>0) == 1)/N0;
            Fn.MDLu(m) = sum(Best.MDLu - double(abs(w)>0) == -1)/N1;
            
            BestModel.MDLcount(m,:) = Best.MDLcount; CostBest.MDLcount(m) = Cost.MDLcount(ibest.MDLcount);
            Posterior.MDLcount(:,m) = posterior.MDLcount;
            Fp.MDLcount(m) = sum(Best.MDLcount - double(abs(w)>0) == 1)/N0;
            Fn.MDLcount(m) = sum(Best.MDLcount - double(abs(w)>0) == -1)/N1;
            
            BestModel.MDLent(m,:) = Best.MDLent; CostBest.MDLent(m) = Cost.MDLent(ibest.MDLent);
            Posterior.MDLent(:,m) = posterior.MDLent;
            Fp.MDLent(m) = sum(Best.MDLent - double(abs(w)>0) == 1)/N0;
            Fn.MDLent(m) = sum(Best.MDLent - double(abs(w)>0) == -1)/N1;
            
            BestModel.L1(m,:) = modelL1; LambdaBestL1(m) = lambda_opt;
            Fp.L1(m) = sum(modelL1 - double(abs(w)>0) == 1)/N0;
            Fn.L1(m) = sum(modelL1 - double(abs(w)>0) == -1)/N1;
            
            % special cases
            if N0 == 0
                Fp.BIC(m) = 0;
                Fp.AIC(m) = 0;
                Fp.MDLl(m) = 0;
                Fp.MDLu(m) = 0;
                Fp.MDLent(m) = 0;
                Fp.MDLcount(m) = 0;
                Fp.L1(m) = 0;
            elseif N1 == 0
                Fn.BIC(m) = 0;
                Fn.AIC(m) = 0;
                Fn.MDLl(m) = 0;
                Fn.MDLu(m) = 0;
                Fn.MDLent(m) = 0;
                Fn.MDLcount(m) = 0;
                Fn.L1(m,:) = 0;
            end
            
            % total error in reconstruction
            Err.AIC(m) = (Fn.AIC(m)*N1 + Fp.AIC(m)*N0)/(N1+N0);
            Err.BIC(m) = (Fn.BIC(m)*N1 + Fp.BIC(m)*N0)/(N1+N0);
            Err.MDLl(m) = (Fn.MDLl(m)*N1 + Fp.MDLl(m)*N0)/(N1+N0);
            Err.MDLu(m) = (Fn.MDLu(m)*N1 + Fp.MDLu(m)*N0)/(N1+N0);
            Err.MDLent(m) = (Fn.MDLent(m)*N1 + Fp.MDLent(m)*N0)/(N1+N0);
            Err.MDLcount(m) = (Fn.MDLcount(m)*N1 + Fp.MDLcount(m)*N0)/(N1+N0);
            Err.L1(m) = (Fn.L1(m)*N1 + Fp.L1(m)*N0)/(N1+N0);
            
        end
    
    fprintf('-----------------------------------------------------  \n\n');
    fprintf('Nx = %d, h = %d, s = %2.4f, J = %2.4f and beta = %2.4f \n',Nx,h,sparsity,theta,beta);
    fprintf('-----------------------------------------------------  \n\n');
    
    time = toc; % end stopwatch and record length
    
    save(['Experiments18_decimation/New_Exp1_Nx',num2str(Nx),'_h',num2str(h),'_s',...
        num2str(sparsity),'_J',num2str(theta),'_LocBeta',num2str(beta),'trial',num2str(trial),'.mat']);
    
    end
end
end