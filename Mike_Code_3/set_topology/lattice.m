%FUNCTION FOR GENERATING A d DIMENSIONAL LATTICE (NEAREST NEIGHBOR
%INTERACTIONS: "2D ISING, 3D ISING ...)


%SI POTREBBE GENERALIZZARE AD UN d LATTICE .. VEDI PER IL MOMENTO d=2, forse
%d=3 funziona bene

%METTI ANCHE UN CONTROLLO PER CHIEDERE CHE LA RADICE D-ESIMA SIA ESATTA 

function [A] = lattice(N,d)
root = round(N^(1/d));
A = zeros(N); 

for i = 1:N-1
    if mod(i, root) ~= 0
        A(i,i+1) = 1;
    end
    j = i + root;
    if i <= N-root
        A(i,j) = 1;
    end
end

%Simmetry of the adjacency matrix
 A = A + A';  

end
% 
% %%
% N = 27;
% sqroot = sqrt(N);
% A = zeros(N);
% for i = 1:N-1;
%     if mod(i, sqroot) ~= 0
%         A(i,i+1) = 1;
%     end
%     j = i + sqroot;
%     if i <= N-sqroot
%         A(i,j) = 1;
%     end
% end
% 
% %Simmetry of the adjacency matrix
%  A = A + A';  
%  
% %%
% cubroot = N^(1/3);
% A = zeros(N);
% for i = 1:N-1;
%     if mod(i, cubroot) ~= 0
%         A(i,i+1) = 1;
%     end
%     
%     j = i + cubroot;
%     if i <= N-cubroot
%         A(i,j) = 1;
%     end
%     
%     k = i + cubroot^2;
%     if i <= N-cubroot^2
%         A(i,k) = 1;
%     end
% end
% 
% %Simmetry of the adjacency matrix
%  A = A + A'; 