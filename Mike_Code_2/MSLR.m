%proto-llh
%

function LLH = MSLR(Nx,h,M,Sstruct)

XX = Sstruct.S_hat;

% number of parameters
Np = Nx + h;

theta = 1;

rep = true;

trial = 1;

RANGO = zeros(numel(M),1);

%% distribution with abs of couplings distributed between 0.5*theta and 1.5*theta | what does "abs" stand for #?
  w = theta*(rand(1,Np)+0.5);					% Creates a vector size of number of parameters, with values randomly generated between 0 and 1, 0.5 added, then whole vector multiplied by theta (1?)... why #?
  segno = 2*double(rand(1,Np)>.5)-1;			% Creates a vector of random +/-1 values. Segno = "sign"
  w = segno.*w;								% Ok so this is the final output vector of random connections from -1.5 to +1.5 values

% sparsity of couplings
% w(randperm(Np,round(sparsity*Np)))=0;			% decimates random values of vector of connection values, not exactly certain how yet. Not sure why using randperm

X = ones(M(numel(M)),Nx); %|For all the observations for all nodes
Y = ones(M(numel(M)),1); %|For one the observations for all nodes
door = true;

% the rank of the composite matrix has to be full rank
while rank([X,Y]) ~= Np+1 || door
    door = false;
    X = XX(randperm(M(numel(M)),M(numel(M))),:);
    
    % if h =1, X will have a ones column attached to X
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
end

% ---- # MODEL SELECTION ON NESTED MODELS THROUGH DECIMATION WITH INFORMATION CRITERIA
% notice 1) that if there are two bests models the routine pick the
% sparser);
% notice 2) if there is a field you don't need to put this
% information in the routine since in that case X already contains
% an additional column
[w_ML,l_ML,posterior,Cost,Best,ibest] = decimation_logistic_Model_Selection(X(sel,:),Y(sel));
%fprintf('Model selection time %4.3f\n',toc)

LLH.w_ML = w_ML;
LLH.l_ML = l_ML;
LLH.posterior = posterior;
LLH.Cost = Cost;
LLH.Best = Best;
LLH.ibest = ibest;

end