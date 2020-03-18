%%%sanitychkferr.m

%% Do sanity check, create clean inputs as defined by Yasser:
				% `Your code should get connectivity matrix J and external field vector h. 
				% What I would do first is just to for instance put to J_ij=J = 0 and h_i = h = 1 for all the spins. 
				% Then I would run the code, generate samples and see what I would get for means and correlations. 
				% If the code is right, you should get m_i = <s_i> = tanh(h) and C_ij = <s_i s_j>-<s_i><s_j> = 0.
				% You can then change h_i to different values while J= 0 And you should get the same result: m_i = tanh(h_i) and C_ij = 0.
				% You can then put J_ij = some constant for some pairs of spins, namely build a network in which spins are connected in pairs (it is called a gas of dimers). 
				% You can analytically solve for m and C (it is a good exercise) and you can also find the solution in many papers (including my PRE paper in 2011 with John Hertz) and then check your simulations results agains the analytical expression. 
				% All these would be sanity checks that there is nothing seriously wrong with the code. 



function sanityferr = sanitychkferr(jn,Sstruct,JHstruct,sparsity,time,T,name)
   
   sanityferr = struct('th',{},'tchk',{},'mtchk',{},'mimj',{},'Cij',{},'mCij',{},'mfC',{},'mfJ',{},'mfh',{},'mferr',{},'tapJ',{},'taph',{},'tapC',{});
    
   for i = 1:jn
    	J = JHstruct(i).Jsparse;
        h = JHstruct(i).Hsparse;
        mi = Sstruct(i).mfinal;
        chi = Sstruct(i).Cfinal;
        S = Sstruct(i).S_hat;

        mimj = mi'*mi;
        Cij = chi - mimj;
        mCij = mean(Cij,'all');

    %% Mean field methods
    %% Forward

        % Check mags using MF method, eq. 51 berg, only if all values of J = 0 does this work
        % Removing for dimer check
        % actual equation is: mfmi = tanh(h + mi*J) fuck with this later
        % Roudi Check for disconnected (J=0) matrix
        tchk = tanh(h) - mi; 
       
    %%For nMF methods of Cij

    % Fix this, need to start with the gen Cij and mi then gen C, then Jmf, then hmf

    % Mean Field methods for MF inversion from Roudi/Hertz 2011
       % Jmk = mi*J;
       % mk = tanh(h+Jmk);
       % mkmj = mk'*mk;
       % mfchi = sisj - mkmj;
       % mfmchi = mean(mfchi,'all');

    %% Ferr magnetization forward mean field method
    % mferr = tanh((J.*mferr' + h)./(T1)) %From eq.3.36 in Amit
    

    N = length(h);
    mferr = rand(1,N)';
    for ii = 1:N
        mferrnew = tanh((J.*0.3*mferr + h')./(1)); %Using T=1 per nicola's instruction also J=0.3
        mferr = mferrnew;
    end
    disp(mferr')
    
    % eq 4.
    %mfbCij = (1-mi.^2)'.*(eye(length(J)) + J*chi); %is it this, not likely says Nicola
    % or is it:
    % mfbCij = diag(1-mi.^2) + J*chi; 
    % mfbCij = diag(1-mi.^2) + (1-mi.^2)'.*J*chi; %This one is my closest guess so far

    % eq 4b. Mean Field inferred Cij from roudi 2009
    % Taking difference of C arrays is stupid as one is a permutation of the other
    % mfC = diag(1-mi.^2) + (1-mi.^2)'.*J*Cij;
    mfC = rand(N,N);
    for iii = 1:N
        mfCnew = diag(1-mferr.^2) + (1-mferr.^2)'.*J*mfC; %This isn't working but who cares rn
        mfC = mfCnew;
    end
    disp(mfC)
        
    % Inferred J 
    Pij = diag(1-mi.^2);
    Jmf = (Pij^-1) - (Cij^-1); % using this for now as this seems closer to what I think it is
    % Jmf = (Pij^-1) - (mfC^-1); % Pretty certain this is wrongm the C from q 4 does not plug into the C %from eq 5

    % Inferred h mean field
    % Jmk = mi*J; % This is wrong this isn't inferred or forward. using the generated mag but the real J
    mfh = atanh(mi) - mi*Jmf; % and then the J here should be the inferred J

   

    %% TAP REconstruction
    Jtap = -2.*(Cij^-1)./(1 + sqrt(1-8.*(Cij^-1).*mimj));
    htap = atan(mi') - Jtap*mi' + mi'.*(Jtap.^2)*(1-mi'.^2); % need to review this formula carefully
    
    %%  Forward Construction from J
    % Ctap = (-J - 2.*(J.^2).*mimj)^-1;
    Ctap = (-J - 2.*(J.^2).*mimj)^-1;

    

    %% differences between calculated and inferred
        %dmfC = abs(Cij(:) - mfC(:));
        %dtapC = abs(Cij(:) - Ctap(:));
        %dmfJ = abs(J(:) - Jmf(:));
        %dtapJ = abs(J(:) - Jtap(:));
        %dmfh = abs(h(:) - mfh(:)); % So normally we would take the difference of original h with MF h but haven't done te mfh yet, see notes above
        %dtaph = abs(h(:) - htap(:)); 


        sanityferr(i).th = tanh(h);
        sanityferr(i).tchk = tchk;
        sanityferr(i).mtchk = mean(tchk);
        sanityferr(i).mimj = mimj;
        sanityferr(i).chi = chi;
        sanityferr(i).Cij = Cij;
        sanityferr(i).mCij  = mCij;
        sanityferr(i).mfC = mfC;
        sanityferr(i).mfJ = Jmf;
        sanityferr(i).mfh = mfh;
        sanityferr(i).mferr = mferr;
        sanityferr(i).tapJ = Jtap;
        sanityferr(i).taph = htap;
        sanityferr(i).tapC = Ctap;
        %sanityferr(i).dmfC = dmfC;
        %sanityferr(i).dtapC = dtapC;
        %sanityferr(i).dmfJ = dmfJ;
        %sanityferr(i).dtapJ = dtapJ;
        %sanityferr(i).dmfh = dmfh;
        %sanityferr(i).dtaph = dtaph;

    end


save([name,'\',time(1:5),'sanityferr_N',num2str(length(h)),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time(6:12),'.mat'],'sanityferr');
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