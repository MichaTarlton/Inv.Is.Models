%BLLH.m
% T: number of observations, length of S
% N: Number of nodes
% h_on,used for field variaable to include h value


function LLH = BLLH(T,N,h_on,SStruct)

LLH = struct('NodeModel',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});

for st = 1:size(SStruct,2) % Structs make sizing weird, this is correct though

S = SStruct(st).S_hat;
%S = S'; %Removed the top level transpose for now

%LLH = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});

disp(['Topology ', num2str(st)])
tic
NodeModel = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});	

 JconBIC   	  = [];
 JconAIC   	  = [];
 JconMDLl  	  = [];
 JconMDLu  	  = [];
 JconMDLent   = [];
 JconMDLcount = [];

 JreconBIC   	= [];
 JreconAIC   	= [];
 JreconMDLl  	= [];
 JreconMDLu  	= [];
 JreconMDLent   = [];
 JreconMDLcount = [];

 thetaBIC      = [];
 thetaAIC      = [];
 thetaMDLl     = [];
 thetaMDLu     = [];
 thetaMDLent   = [];
 thetaMDLcount = [];


	for t = 1:N
		
		disp(['T', num2str(st), ' Regressing Node ', num2str(t)])
		
		X = S(:,setdiff(1:N,t));
		Y = S(:,t);
		
		[w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,h_on);
		
		%LLH(st).NodeModel(t).w_ML 	 = w_ML;
		%LLH(st).NodeModel(t).l_ML		 = l_ML;
		%LLH(st).NodeModel(t).posterior = posterior;
		%LLH(st).NodeModel(t).cost		 = cost;
		%LLH(st).NodeModel(t).BestModel = BestModel;
		%LLH(st).NodeModel(t).IMAX		 = IMAX;
        disp('check 1')
		NodeModel(t).w_ML 	 = w_ML;
		NodeModel(t).l_ML		 = l_ML;
		NodeModel(t).posterior = posterior;
		NodeModel(t).cost		 = cost;
		NodeModel(t).BestModel = BestModel;
		NodeModel(t).IMAX		 = IMAX;
disp('check 2')
        JconBIC(t,:)    = BestModel.BIC;
        JconAIC(t,:)    = BestModel.AIC;
        JconMDLl(t,:)   = BestModel.MDLl;
        JconMDLu(t,:)   = BestModel.MDLu;
        JconMDLent(t,:)   = BestModel.MDLent; % making the connection graph at an earlier juncture
        JconMDLcount(t,:)  = BestModel.MDLcount;
disp('check 3')        
        thetaBIC(t,:)    = w_ML(IMAX.BIC,:);
		thetaAIC(t,:)    = w_ML(IMAX.AIC,:);
		thetaMDLl(t,:)   = w_ML(IMAX.MDLl,:);
		thetaMDLu(t,:)   = w_ML(IMAX.MDLu,:);
		thetaMDLent(t,:) = w_ML(IMAX.MDLent,:); % this also works I think
		thetaMDLcount(t,:) = w_ML(IMAX.MDLcount,:);
        %theta(t,:)       = LLH(t).w_ML(LLH(t).IMAX.MDLent,:); %IMAX is index of best model

		
		%LLH.JHrecon = [;LLH(t).theta]; 
	end
toc
LLH(st).NodeModel = NodeModel;


JconBIC = JconBIC(:,2:end);
JconAIC = JconAIC(:,2:end);
JconMDLl = JconMDLl(:,2:end);
JconMDLu = JconMDLu(:,2:end);
JconMDLent = JconMDLent(:,2:end);
JconMDLcount = JconMDLcount(:,2:end);

LLH(st).Jcon.BIC      = [zeros(N,1), triu(JconBIC)] +      [tril(JconBIC,-1),zeros(N,1)];
LLH(st).Jcon.AIC      = [zeros(N,1), triu(JconAIC)] +      [tril(JconAIC,-1),zeros(N,1)];
LLH(st).Jcon.MDLl     = [zeros(N,1), triu(JconMDLl)] +     [tril(JconMDLl,-1),zeros(N,1)];
LLH(st).Jcon.MDLu     = [zeros(N,1), triu(JconMDLu)] +     [tril(JconMDLu,-1),zeros(N,1)];
LLH(st).Jcon.MDLent   = [zeros(N,1), triu(JconMDLent)] +   [tril(JconMDLent,-1),zeros(N,1)];
LLH(st).Jcon.MDLcount = [zeros(N,1), triu(JconMDLcount)] + [tril(JconMDLcount,-1),zeros(N,1)];

LLH(st).theta.BIC      = thetaBIC;
LLH(st).theta.AIC      = thetaAIC;
LLH(st).theta.MDLl     = thetaMDLl;
LLH(st).theta.MDLu     = thetaMDLu;
LLH(st).theta.MDLent   = thetaMDLent;
LLH(st).theta.MDLcount = thetaMDLcount;


%LLH(st).theta 	 =  theta;

LLH(st).hrecon.BIC        =  thetaBIC(:,1);
LLH(st).hrecon.AIC        =  thetaAIC(:,1);
LLH(st).hrecon.MDLl       =  thetaMDLl(:,1);
LLH(st).hrecon.MDLu       =  thetaMDLu(:,1);
LLH(st).hrecon.MDLent     =  thetaMDLent(:,1);
LLH(st).hrecon.MDLcount   =  thetaMDLcount(:,1);
%LLH(st).hrecon   =  theta(:,1);

JreconBIC   =  thetaBIC(:,2:end);  
JreconAIC   =  thetaAIC(:,2:end);  
JreconMDLl   =  thetaMDLl(:,2:end);  
JreconMDLu   =  thetaMDLu(:,2:end);  
JreconMDLent   =  thetaMDLent(:,2:end);  
JreconMDLcount   =  thetaMDLcount(:,2:end);  
%Jrecon   =  theta(:,2:end);  

LLH(st).Jrecon.BIC      = [zeros(N,1), triu(JreconBIC)] + [tril(JreconBIC,-1),zeros(N,1)];
LLH(st).Jrecon.AIC      = [zeros(N,1), triu(JreconAIC)] + [tril(JreconAIC,-1),zeros(N,1)];
LLH(st).Jrecon.MDLl     = [zeros(N,1), triu(JreconMDLl)] + [tril(JreconMDLl,-1),zeros(N,1)];
LLH(st).Jrecon.MDLu     = [zeros(N,1), triu(JreconMDLu)] + [tril(JreconMDLu,-1),zeros(N,1)];
LLH(st).Jrecon.MDLent   = [zeros(N,1), triu(JreconMDLent)] + [tril(JreconMDLent,-1),zeros(N,1)];
LLH(st).Jrecon.MDLcount = [zeros(N,1), triu(JreconMDLcount)] + [tril(JreconMDLcount,-1),zeros(N,1)];
%LLH(st).Jrecon = [zeros(N,1), triu(Jrecon)] + [tril(Jrecon,-1),zeros(N,1)];



%trying to find which ones are connected inf v tru
% Jcon = double(not(LLH(1).Jrecon == 0)) % find the connected ones (couldalso jsut grab the best models vector)
%Jcon = Jcon(:,2:end); % Yeah this doesn'twork because you haven'tadded the zero diag yet
%Jcon = [zeros(N,1), triu(Jcon)] + [tril(Jcon,-1),zeros(N,1)];
%LLH(st).Jcon = Jcon;

LLH(st).Jasym.BIC      = tril(LLH(st).Jcon.BIC     )' + triu(LLH(st).Jcon.BIC     );
LLH(st).Jasym.AIC      = tril(LLH(st).Jcon.AIC     )' + triu(LLH(st).Jcon.AIC     );
LLH(st).Jasym.MDLl     = tril(LLH(st).Jcon.MDLl    )' + triu(LLH(st).Jcon.MDLl    );
LLH(st).Jasym.MDLu     = tril(LLH(st).Jcon.MDLu    )' + triu(LLH(st).Jcon.MDLu    );
LLH(st).Jasym.MDLent   = tril(LLH(st).Jcon.MDLent  )' + triu(LLH(st).Jcon.MDLent  );
LLH(st).Jasym.MDLcount = tril(LLH(st).Jcon.MDLcount)' + triu(LLH(st).Jcon.MDLcount);					
%LLH(st).Jasym = tril(Jcon)' + triu(Jcon);		% Checking to see if duplicate values exist (first try says yes)

LLH(st).symrate.BIC      = sum(double(LLH(st).Jasym.BIC      == 2),'all');
LLH(st).symrate.AIC      = sum(double(LLH(st).Jasym.AIC      == 2),'all');
LLH(st).symrate.MDLl     = sum(double(LLH(st).Jasym.MDLl     == 2),'all');
LLH(st).symrate.MDLu     = sum(double(LLH(st).Jasym.MDLu     == 2),'all');
LLH(st).symrate.MDLent   = sum(double(LLH(st).Jasym.MDLent   == 2),'all');
LLH(st).symrate.MDLcount = sum(double(LLH(st).Jasym.MDLcount == 2),'all');
%LLH(st).symrate = sum(double(LLH(1).Jasym == 2),'all');


LLH(st).asymrate.BIC      = sum(double(LLH(st).Jasym.BIC      == 1),'all');
LLH(st).asymrate.AIC      = sum(double(LLH(st).Jasym.AIC      == 1),'all');
LLH(st).asymrate.MDLl     = sum(double(LLH(st).Jasym.MDLl     == 1),'all');
LLH(st).asymrate.MDLu     = sum(double(LLH(st).Jasym.MDLu     == 1),'all');
LLH(st).asymrate.MDLent   = sum(double(LLH(st).Jasym.MDLent   == 1),'all');
LLH(st).asymrate.MDLcount = sum(double(LLH(st).Jasym.MDLcount == 1),'all');
%LLH(st).asymrate = sum(double(LLH(1).Jasym == 1),'all');

%Tru positives




end





%ok luckily though at a glance the values seem similar,so I may be able to just avg them out

% Ok now for comparing original connectivity to recon
% this will have to be done at a higher level, but writing here to keep track
% so the problm is that there are no zero values in my original values unless I set it to sparsify or disconnect values above a certain threshold
%
% Jcon = LLH(1).Jcon
% Jcon = double(not(Jcon == 0))
% 


end