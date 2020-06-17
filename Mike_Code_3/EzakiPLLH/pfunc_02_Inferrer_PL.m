% 2016/11/11 by T. Ezaki 
% This function infers h and J by pseudo-likelihood maximization.
% binarizedData: nodeNumber x tmax

% Making edits to work with my code
% Container of function overhauled by me at this point

%function [h,J] = pfunc_02_Inferrer_PL(binarizedData)
%function  PLLH = pfunc_02_Inferrer_PL(binarizedData,name,sname,time,beta)
function  PLLH = pfunc_02_Inferrer_PL(binarizedData,name,time,beta)

%PLLH = struct('h',{},'J',{},'name',{});
PLLH = struct;

%for ii = 1:length(Sstructvec)
    %sname=Snamevec(ii);  
    %binarizedData = Sstructvec(ii).S_hat';
    


    [nodeNumber,dataLength] = size(binarizedData);
    %iterationMax = 5000000;
    iterationMax = 100000;
    permissibleErr = 0.00000001;
    dt = 0.2;
    
    dataMean = mean(binarizedData,2);
    dataCorrelation = (binarizedData*binarizedData')/dataLength;
    
    h = zeros(nodeNumber,1);
    J = zeros(nodeNumber);
    
        for t=1:iterationMax
            %calculate Ldh
            Ldh = - dataMean + mean(tanh(J * binarizedData + h * ones(1, dataLength)),2);
            h =  h - dt * Ldh;
            
            
            LdJ = -dataCorrelation + 0.5 * binarizedData * (tanh(J * binarizedData + h * ones(1, dataLength)))'/dataLength...
                 +0.5 * (binarizedData * (tanh(J * binarizedData + h * ones(1, dataLength)))')'/dataLength;
                
            LdJ = LdJ - diag(diag(LdJ));%eliminate diagonal entities (J_{ii}=0)
            J = J - dt * LdJ;
             
            if (sqrt(norm(LdJ,'fro')^2 + norm(Ldh)^2)/nodeNumber/(nodeNumber+1) < permissibleErr)
                break;
            end
        end
    PLLH.h = h;
    PLLH.J = J;
    PLLH.name = name;
end
    %save([name,'\',time(1:5),'PLLH_','_beta',num2str(beta),'_',time(6:12),'.mat'],'PLLH');
%end

