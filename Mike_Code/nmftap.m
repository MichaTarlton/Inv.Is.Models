%% nmftap
% Not sure what i'm doing but I though I would go ahead and make the m_i magnetization calculation


function nmftaps=nmftap(J,h,S)
nmftaps = struct('mk',{},'Jmk',{},'mi',{},'Cij',{});


mk = mean(S);
Jmk = mk*J; 

% make sure h and Jmk are in same direction and size
% this probably doesn't work since j != i and I don't think thisimplementation accoutn for that, however the 0 value in the J_ij might actually account for it



 mi = tanh(h+Jmk); % should result in a vector 

%% Cij = diff(mi,h); %dumbass way to do it

Cij = mi - h';

nmftaps(1).mk = mk;
nmftaps(1).mi = mi;
nmftaps(1).Jmk = Jmk;
nmftaps(1).Cij = Cij;

 %---------------------------------------------
 %% Makinf eq.3 from Roudi 2011