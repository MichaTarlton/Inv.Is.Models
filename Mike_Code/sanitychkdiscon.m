%%%sanitychk.m

%% Do sanity check, create clean inputs as defined by Yasser:
				% `Your code should get connectivity matrix J and external field vector h. 
				% What I would do first is just to for instance put to J_ij=J = 0 and h_i = h = 1 for all the spins. 
				% Then I would run the code, generate samples and see what I would get for means and correlations. 
				% If the code is right, you should get m_i = <s_i> = tanh(h) and C_ij = <s_i s_j>-<s_i><s_j> = 0.
				% You can then change h_i to different values while J= 0 And you should get the same result: m_i = tanh(h_i) and C_ij = 0.
				% You can then put J_ij = some constant for some pairs of spins, namely build a network in which spins are connected in pairs (it is called a gas of dimers). 
				% You can analytically solve for m and C (it is a good exercise) and you can also find the solution in many papers (including my PRE paper in 2011 with John Hertz) and then check your simulations results agains the analytical expression. 
				% All these would be sanity checks that there is nothing seriously wrong with the code. 



function sanitydiscon = sanitychkdiscon(jn,Sstruct,JHstruct,sparsity,time,T)
   
   sanitydiscon = struct('th',{},'tchk',{},'mtchk',{},'mimj',{},'Cij',{},'mCij',{},'mfC',{},'mfJ',{},'mfh',{},'tapJ',{},'taph',{},'tapC',{});
    
   for i = 1:jn
    	J = JHstruct(i).Jsparse;
        h = JHstruct(i).Hsparse;
        mi = Sstruct(i).mfinal;
        chi = Sstruct(i).Cfinal;
        S = Sstruct(i).S_hat;

          	
        % Check mags using MF method, eq. 51 berg, only if all values of J = 0 does this work
        % Removing for dimer check
        % actual equation is: mfmi = tanh(h + mi*J) fuck with this later
        % 
        tchk = tanh(h) - mi; 
        mtchk = mean(tchk);

        mimj = mi'*mi;
        Cij = chi - mimj;
        mCij = mean(Cij,'all');

        % Don't need this rn
        % if abs(mtchk) < 0.01
        %     saneh = 1;
        % else
        %     saneh = 0;
        % end
        % 
        % if abs(mCij) < 0.01
        %     saneCij = 1;
        % else
        %     saneCij = 0;
        % end

    %%For nMF methods of Cij
    % Mean Field methods for MF inversion from Roudi/Hertz 2011
       % Jmk = mi*J;
       % mk = tanh(h+Jmk);
       % mkmj = mk'*mk;
       % mfchi = sisj - mkmj;
       % mfmchi = mean(mfchi,'all');

    % Inferred h
    Jmk = mi*J;
    mfh = atanh(mi) - Jmk;
        
    % Inferred J 
    Pij = diag(1-mi.^2);
    Jmf = (Pij.^-1) - (Cij.^-1);

    % eq 4b. Mean Field inferred Cij from roudi 2009
    mfC = diag(1-mi.^2) + (1-mi.^2)'.*J*Cij; 
    
    % eq 4.
    %mfbCij = (1-mi.^2)'.*(eye(length(J)) + J*chi); %is it this, not likely says Nicola
    % or is it:
    % mfbCij = diag(1-mi.^2) + J*chi; 
    % mfbCij = diag(1-mi.^2) + (1-mi.^2)'.*J*chi; %This one is my closest guess so far


    %% TAP REconstruction
    Jtap = -2.*(Cij.^-1)./(1 + sqrt(1-8.*(Cij.^-1).*mimj));
    htap = atan(mi') - Jtap*mi' + mi'.*(Jtap.^2)*(1-mi'.^2);
    
    %%  Forward Construction from J
    Ctap = (-J - 2.*(J.^2).*mimj).^-1;

    %% differences between calculated and inferred
        dmfC = abs(Cij(:) - mfC(:));
        dtapC = abs(Cij(:) - Ctap(:));
        dmfJ = abs(J(:) - Jmf(:));
        dtapJ = abs(J(:) - Jtap(:));
        dmfh = abs(h(:) - mfh(:)); % So normally we would take the difference of original h with MF h but haven't done te mfh yet, see notes above
        dtaph = abs(h(:) - htap(:)); 


        sanitydiscon(i).th = tanh(h);
        sanitydiscon(i).tchk = tchk;
        sanitydiscon(i).mtchk = mtchk;
        sanitydiscon(i).mimj = mimj;
        sanitydiscon(i).chi = chi;
        sanitydiscon(i).Cij = Cij;
        sanitydiscon(i).mCij  = mCij;
        sanitydiscon(i).mfC = mfC;
        sanitydiscon(i).mfJ = Jmf;
        sanitydiscon(i).mfh = mfh;
        sanitydiscon(i).tapJ = Jtap;
        sanitydiscon(i).taph = htap;
        sanitydiscon(i).tapC = Ctap;
        sanitydiscon(i).dmfC = dmfC;
        sanitydiscon(i).dtapC = dtapC;
        sanitydiscon(i).dmfJ = dmfJ;
        sanitydiscon(i).dtapJ = dtapJ;
        sanitydiscon(i).dmfh = dmfh;
        sanitydiscon(i).dtaph = dtaph;

    end


save([time(1:5),'sanitydiscon_N',num2str(length(h)),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time(6:12),'.mat'],'sanitydiscon');
%%save(['sanity_N',num2str(length(h)),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanity');
%%save(['sanity_N',num2str(N),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanity');
end