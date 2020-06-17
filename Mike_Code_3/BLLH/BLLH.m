%BLLH.m
% T: number of observations, length of S
% N: Number of nodes
% h_on,used for field variaable to include h value


function LLH = BLLH(T,N,h_on,S);


%for t = 1:T
S = S';
LLH = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});

	for t = 1:N
	
		X = S(:,setdiff(1:N,t));
		Y = S(:,t);
		
		[w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,h_on);
		
		LLH(t).w_ML 	 = w_ML;
		LLH(t).l_ML		 = l_ML;
		LLH(t).posterior = posterior;
		LLH(t).cost		 = cost;
		LLH(t).BestModel = BestModel;
		LLH(t).IMAX		 = IMAX;


        Jcon(t,:) 		= LLH(t).BestModel.MDLent; % making the connection graph at an earlier juncture
        theta(t,:)       = LLH(t).w_ML(IMAX.MDLent,:); % this also works I think
        %theta(t,:)       = LLH(t).w_ML(LLH(t).IMAX.MDLent,:); %IMAX is index of best model

		%disp(t)
		%LLH.JHrecon = [;LLH(t).theta]; 
	end

LLH(1).theta 	 =  theta;
LLH(1).hrecon   =  theta(:,1);
Jrecon   =  theta(:,2:end);  
LLH(1).Jrecon = [zeros(N,1), triu(Jrecon)] + [tril(Jrecon,-1),zeros(N,1)];

%trying to find which ones are connected inf v tru
% Jcon = double(not(LLH(1).Jrecon == 0)) % find the connected ones (couldalso jsut grab the best models vector)
Jcon = Jcon(:,2:end); % Yeah this doesn'twork because you haven'tadded the zero diag yet
Jcon = [zeros(N,1), triu(Jcon)] + [tril(Jcon,-1),zeros(N,1)];
LLH(1).Jcon = Jcon;
					
LLH(1).Jasym = tril(Jcon)' + triu(Jcon);		% Checking to see if duplicate values exist (first try says yes)

LLH(1).symrate = sum(double(LLH(1).Jasym == 2),'all');
LLH(1).asymrate = sum(double(LLH(1).Jasym == 1),'all');



%ok luckily though at a glance the values seem similar,so I may be able to just avg them out

% Ok now for comparing original connectivity to recon
% this will have to be done at a higher level, but writing here to keep track
% so the problm is that there are no zero values in my original values unless I set it to sparsify or disconnect values above a certain threshold
%
% Jcon = LLH(1).Jcon
% Jcon = double(not(Jcon == 0))
% 


end