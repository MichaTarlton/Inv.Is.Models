%%% EXACT PARTITION FUNCTION, MAGNETIZATIONS AND CORRELATIONS%%%
% suitable to use with very small networks

%INPUT : parameters (fields and couplings)
%OUTPUT : partition function, magnetisation and correlations

% Note: given the output you can calculate the likelihood if needed;

% Note: it works for both -1/1 and 0/1 representations;

function [Z,m,C] = Exact_Ising(h,J,vw,representation)
tic

% show plots if c is true
try c = vw;catch
    c = false;
end

% select -1/1 representation if rep = true; 0/1 otherwise
try rep = representation; catch
    rep = true;
end

% initialization
N = numel(h);
m = zeros(1,N);
C = zeros(N);
Z = 0;

%%% MAIN LOOP OVER ALL THE CONFIGURATION %%%
for i = 0:2^N-1
    
    if rep
        S = 2*(dec2bin(i,N) - '0')-1;     %Configurations
    else
        S = (dec2bin(i,N) - '0');
    end
    
    E = -0.5*S*J*S' - S*h';           %Energy
     
    W = exp(-E);                      %Boltzmann weight
    
    Z = Z + W;                        %Partition function
    
    m = m + S*W;                      %Magnetization
    
    C = C + S'*S*W;                   %Correlations
    
    if c == 1
        display([i toc])
    end
end

m = m/Z;
C = C/Z - m'*m;

%PLOT OF MAGNETIZATION AND CORRELATIONS IN THE LAST SAMPLE
if c ==1

width = 0.5;
figure;
bar(1:N,m,width)
title('Magnetizations')
ylabel('m(i)')
xlabel('i')

figure;
bar3(C,width)
title('Correlations')
zlabel('C(i,j)')
xlabel('i')
ylabel('j')

end
toc