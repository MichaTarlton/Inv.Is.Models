%BLLH.m
% T: number of observations, length of S
% N: Number of nodes
% h_on,used for field variaable to include h value


function LLH = BLLH(T,N,h_on,S);


%for t = 1:T
S = S';
LLH = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{});

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
	
	end
end