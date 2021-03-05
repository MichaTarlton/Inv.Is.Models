% DILUTION it applies to both fully and non fully connected graph
% p is the fraction of couplings randomly set to zero


function [Adj] = dilution(Adj,p)

Nspin = length(Adj);
Ncouples = sum(sum(triu(Adj)));
Ndilution = round(p*Ncouples);
idx = zeros(Ncouples,2);

% vectorization of the couples
k = 1;
for i = 1:Nspin;
    j = 1;
    while j<i
        if Adj(i,j) == 1;
          idx(k,:) = [i j];
          k = k +1;
        end
          j = j +1;
    end
end

for q = 0:Ndilution-1;
    ii = randi(Ncouples-q);
    Adj(idx(ii,1),idx(ii,2)) = 0;
    Adj(idx(ii,2),idx(ii,1)) = 0;
    idx(ii,:) = [];
end

end
