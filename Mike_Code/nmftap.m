%% nmftap
% Not sure what i'm doing but I though I would go ahead and make the m_i magnetization calculation




nmftap(J,h,S)



mk = mean(S);
Jmk = J*mk; 

% make sure h and Jmk are in same direction and size
% this probably doesn't work since j != i and I don't think thisimplementation accoutn for that, however the 0 value in the J_ij might actually account for it



 mi = tanh(h+Jmk) % should result in a vector 

Cij = diff(mi,h)




 %---------------------------------------------
 %% Makinf eq.3 from Roudi 2011