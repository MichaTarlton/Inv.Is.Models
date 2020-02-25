%% nmftap
% Not sure what i'm doing but I though I would go ahead and make the m_i magnetization calculation


function nmftaps=nmftap(J,h,S)
nmftaps = struct('mk',{},'Jmk',{},'mi',{},'Cij',{});

J = JHstruct(i).Jsparse;
h = JHstruct(i).Hsparse;
mi = Sstruct(i).mfinal;
chi = Sstruct(i).Cfinal; % Cfinal = S_hat'*S_hat/T;

mk = mean(S);
Jmk = mk*J; 

% make sure h and Jmk are in same direction and size
% this probably doesn't work since j != i and I don't think thisimplementation accoutn for that, however the 0 value in the J_ij might actually account for it



%mi = tanh(h+Jmk); % should result in a vector 

%% Cij = diff(mi,h); %dumbass way to do it

% Cunknonw = mi - h';
% 
% nmftaps(1).mk = mk;
% nmftaps(1).mi = mi;
% nmftaps(1).Jmk = Jmk;
% nmftaps(1).Cij = Cij;

%---------------------------------------------
%% TAP Reconstruction following Berg 2017 eq 57
	% For regular Cij
    tchk = tanh(h) - mi; 
    mtchk = mean(tchk);

    mimj = mi'*mi;
    Cij = chi - mimj;
    
 	Jtap = -2.*(Cij.^-1)./(1 + sqrt(1-8.*(Cij.^-1).*mimj))

 	htap = atan(mi') - Jtap*mi' + mi'.*(Jtap.^2)*(1-mi'.^2)

	%%	Forward Construction from J
 	Ctap = (-J - 2.*(J.^2).*mimj).^-1

%% Sessak-Monassson Reconstruction
	
	% J indie pair
	% abandon for now as you have to port in the spins and put them in the calculation



	% Jsm = -(Cij.^-1) + Jip - (Cij ./ (1-mi.^2)'*(1-mi.^2)(Cij.^2))