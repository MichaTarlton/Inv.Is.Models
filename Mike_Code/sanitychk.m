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



function sanity = sanitychk(jn,Sstruct,JHstruct,sparsity,time,T)
   sanity = struct('th',{},'tchk',{},'mtchk',{},'mimj',{},'chi',{},'mchi',{},'saneh',{},'sanechi',{});
   for i = 1:jn
    	h = JHstruct(i).Hsparse;
    	mi = Sstruct(i).mfinal;
    	sisj = Sstruct(i).Cfinal;
    	
        tchk = tanh(h) - mi; 
        mtchk = mean(tchk);

        mimj = mi'*mi;
    	chi = sisj - mimj;
    	mchi = mean(chi,'all');

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

    	sanity(i).th = tanh(h);
        sanity(i).tchk = tchk;
        sanity(i).mtchk = mtchk;
    	sanity(i).mimj = mimj;
    	sanity(i).chi = chi;
    	sanity(i).mchi  = mchi;
        sanity(i).saneh = saneh;
    	sanity(i).sanechi = sanechi;

    end

save(['sanity_N',num2str(length(h)),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanity');
%%save(['sanity_N',num2str(N),'_T',num2str(T),'_trials',num2str(jn),'_',num2str(100*sparsity),'_',time,'.mat'],'sanity');
end