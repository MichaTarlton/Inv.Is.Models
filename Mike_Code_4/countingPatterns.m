% 7/06/18

% COUNTING PATTERNS AND CALCULATING ENTROPY FOR MODEL SELECTION
% given the binary -1/1 input it counts how many patterns which are
% different from each other also by parity operation;
% it also calculate the entropy of the patterns without adjusting
% for parity


function [m,Hs,nu,m_nu] = countingPatterns(s,chekk)

try Check = chekk; catch
    Check = false;
end

% CALCULATE # OF PATTERNS m (PATTERNS DIFFERENT BY PARITY CONSIDERED EQUAL)
[M,N] = size(s);
randvec = 10*randn(N,1);
vec = s*randvec; % same result of the matrix product means same configuration
[clustvecp,~] = sort(abs(vec)); % put equal configurations togheter (regardless of parity)
[Kp, ~] = sort(diff(find([1,diff(clustvecp'),1]))); % K = frequency of states: each place where a transition occurs, how long is it?
m = numel(Kp); % # of patterns discounted by parity

% CALCULATE # OF PATTERNS m (PATTERNS DIFFERENT BY PARITY CONSIDERED DIFFERENT)
[clustvec,~] = sort(vec); % put equal configurations togheter
[K, ~] = sort(diff(find([1,diff(clustvec'),1]))); % K = frequency of states: each place where a transition occurs, how long is it?
nu = K(logical([1,diff(K)]))/M; % unique frequencies
m_nu = diff(find([1,diff(K),1])); % mK = frequency of frequencies

% CALCULATE ENTROPY
Hs = -sum((nu.*m_nu).*log2(nu));

if Check
% %%%%%%%%%%%% ANOTHER WAY TO FIND nu and m_nu for comparison (slower)
[~,ia,ic]=unique(s,'rows');
K2 = sort(histcounts(ic,numel(ia)));
[~,iaKu,icKu]=unique(K2);
m_nu2 = histcounts(icKu,numel(iaKu));

if sum(K-K2)~=0 || sum(m_nu-m_nu2)~=0
    disp([K;K2;m_nu;m_nu2]);
    error('Something wrong!!');
end

end