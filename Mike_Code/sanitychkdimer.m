%%%sanitychkdimer.m

%% Do sanity check, create clean inputs as defined by Yasser:
				% `Your code should get connectivity matrix J and external field vector h. 
				% What I would do first is just to for instance put to J_ij=J = 0 and h_i = h = 1 for all the spins. 
				% Then I would run the code, generate samples and see what I would get for means and correlations. 
				% If the code is right, you should get m_i = <s_i> = tanh(h) and C_ij = <s_i s_j>-<s_i><s_j> = 0.
				% You can then change h_i to different values while J= 0 And you should get the same result: m_i = tanh(h_i) and C_ij = 0.
				% You can then put J_ij = some constant for some pairs of spins, namely build a network in which spins are connected in pairs (it is called a gas of dimers). 
				% You can analytically solve for m and C (it is a good exercise) and you can also find the solution in many papers (including my PRE paper in 2011 with John Hertz) and then check your simulations results agains the analytical expression. 
				% All these would be sanity checks that there is nothing seriously wrong with the code. 



function sanitydimer = sanitychkdimer(jn,Sstruct,JHstruct,sparsity,time,T)
   sanitydimer = struct('th',{},'tchk',{},'mtchk',{},'saneh',{},'mimj',{},'chi',{},'mchi',{},'sanechi',{},
                        'Jmk',{},'mk',{},'mfchi',{},'mfmchi',{},
                        'mfaCij',{},'mfbCij',{}) 
                %struct('th',{},'tchk',{},'mtchk',{},'mimj',{},'chi',{},'mchi',{},'saneh',{},'sanechi',{});
   for i = 1:jn
    	
        J = JHstructDimer(i).Jsparse;
        h = JHstructDimer(i).Hsparse;
    	mi = SstructDimer(i).mfinal;
    	sisj = SstructDimer(i).Cfinal; % Cfinal = S_hat'*S_hat/T;
    	
        % For regular Cij
        tchk = tanh(h) - mi; 
        mtchk = mean(tchk);

        mimj = mi'*mi;
    	chi = sisj - mimj;
    	mchi = mean(chi,'all');

        %%For nMF methods of Cij
        % Brute force
        Jmk = mi*J;
        mk = tanh(h+Jmk);
        mkmj = mk'*mk;
        mfchi = sisj - mkmj;
        mfmchi = mean(mfchi,'all');

        % eq 4b.
        mfbCij = (1-mi.^2)'.*(eye(length(J)) + J*chi); %is it this, not likely says Nicola
        
        % or is it:
        % mfbCij = diag(1-mi.^2) + J*chi; 
        % mfbCij = diag(1-mi.^2) + (1-mi.^2)'.*J*chi; %This one is my closest guess so far




        if abs(mtchk) < 0.01
            saneh = 1;
        else
            saneh = 0;
        end

    	if abs(mchi) < 0.01
    		sanechi = 1;
    	else
    		sanechi = 0;
    	end

    	sanitydimer(i).th = tanh(h);
        sanitydimer(i).tchk = tchk;
        sanitydimer(i).mtchk = mtchk;
    	sanitydimer(i).mimj = mimj;
    	sanitydimer(i).chi = chi;
    	sanitydimer(i).mchi  = mchi;
        sanitydimer(i).saneh = saneh;
    	sanitydimer(i).sanechi = sanechi;
        sanitydimer(i).Jmk = Jmk;
        sanitydimer(i).mk = mk;
        sanitydimer(i).mfchi = mfchi;
        sanitydimer(i).mfmchi = mfmchi;
        sanitydimer(i).mfbCij = mfbCij;


    end

save(['sanitydimer_N',num2str(length(h)),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanitydimer');
%%save(['sanitydimer_N',num2str(N),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanitydimer');
end