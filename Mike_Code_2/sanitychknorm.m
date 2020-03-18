%%%sanitychknorm.m

%% Do sanity check, create clean inputs as defined by Yasser:
				% `Your code should get connectivity matrix J and external field vector h. 
				% What I would do first is just to for instance put to J_ij=J = 0 and h_i = h = 1 for all the spins. 
				% Then I would run the code, generate samples and see what I would get for means and correlations. 
				% If the code is right, you should get m_i = <s_i> = tanh(h) and C_ij = <s_i s_j>-<s_i><s_j> = 0.
				% You can then change h_i to different values while J= 0 And you should get the same result: m_i = tanh(h_i) and C_ij = 0.
				% You can then put J_ij = some constant for some pairs of spins, namely build a network in which spins are connected in pairs (it is called a gas of dimers). 
				% You can analytically solve for m and C (it is a good exercise) and you can also find the solution in many papers (including my PRE paper in 2011 with John Hertz) and then check your simulations results agains the analytical expression. 
				% All these would be sanity checks that there is nothing seriously wrong with the code. 



function sanitynorm = sanitychknorm(jn,Sstruct,JHstruct,sparsity,time,T,name)
   
   sanitynorm = struct('th',{},'tchk',{},'mtchk',{},'mimj',{},'Cij',{},'mCij',{},'mfC',{},'mfJ',{},'mfh',{},'tapJ',{},'taph',{},'tapC',{});
    
   for i = 1:jn
    	J = JHstruct(i).Jsparse;
        h = JHstruct(i).Hsparse;
        mi = Sstruct(i).mfinal;
        chi = Sstruct(i).Cfinal;
        S = Sstruct(i).S_hat;

        mimj = mi'*mi;
        Cij = chi - mimj;
        mCij = mean(Cij,'all');

    %% Forward construction
    
    % nMF Construction

    % Mean Field methods for MF inversion from Roudi/Hertz 2011
    % Check mags using MF method, eq. 51 berg, equation is: mfmi = tanh(h + mi*J)
    % Could also invert Reconstruction equation below

       % Jmk = mi*J;
       % mk = tanh(h+Jmk);
       % mkmj = mk'*mk;
       % mfchi = sisj - mkmj;
       % mfmchi = mean(mfchi,'all');

    % eq 4b. Mean Field inferred Cij from roudi 2009
    % Taking difference of C arrays is stupid as one is a permutation of the other
    % mfC = diag(1-mi.^2) + (1-mi.^2)'.*J*Cij;
    

    % TAP Construction
    %  Forward Construction from J
    Ctap = (-J - 2.*(J^2).*mimj)^-1;



    %% Infered Reconstruction

    % nMF
    % Inferred J 
    Pij = diag(1-mi.^2);
    Jmf = (Pij^-1) - (Cij^-1); % using this for now as this seems closer to what I think it is

    % Inferred h
    mfh = atanh(mi) - mi*Jmf;

    % TAP REconstruction
    Jtap = -2.*(Cij^-1)./(1 + sqrt(1-8.*(Cij^-1).*mimj));
    htap = atan(mi') - Jtap*mi' + mi'.*(Jtap^2)*(1-mi'.^2); % need to review this formula carefully


    %% differences between calculated and inferred
        %dmfC = abs(Cij(:) - mfC(:));
        dtapC = abs(Cij(:) - Ctap(:));
        dmfJ = abs(J(:) - Jmf(:));
        dtapJ = abs(J(:) - Jtap(:));
        dmfh = abs(h(:) - mfh(:)); % So normally we would take the difference of original h with MF h but haven't done te mfh yet, see notes above
        dtaph = abs(h(:) - htap(:)); 


        sanitynorm(i).th = tanh(h);
        %sanitynorm(i).tchk = tchk;
        %sanitynorm(i).mtchk = mean(tchk);
        sanitynorm(i).mimj = mimj;
        sanitynorm(i).chi = chi;
        sanitynorm(i).Cij = Cij;
        sanitynorm(i).mCij  = mCij;
        %sanitynorm(i).mfC = mfC;
        sanitynorm(i).mfJ = Jmf;
        sanitynorm(i).mfh = mfh;
        sanitynorm(i).tapJ = Jtap;
        sanitynorm(i).taph = htap;
        sanitynorm(i).tapC = Ctap;
        %sanitynorm(i).dmfC = dmfC;
        sanitynorm(i).dtapC = dtapC;
        sanitynorm(i).dmfJ = dmfJ;
        sanitynorm(i).dtapJ = dtapJ;
        sanitynorm(i).dmfh = dmfh;
        sanitynorm(i).dtaph = dtaph;

    end


save([name,'\',time(1:5),'sanitynorm_N',num2str(length(h)),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time(6:12),'.mat'],'sanitynorm');
%%save(['sanity_N',num2str(length(h)),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanity');
%%save(['sanity_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanity');
end


% Vestigial snips

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