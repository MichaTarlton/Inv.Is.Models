%%%JHs.m Disconnect matrix

%% Do sanity check, create clean inputs as defined by Yasser:
				% `Your code should get connectivity matrix J and external field vector h. 
				% What I would do first is just to for instance put to J_ij=J = 0 and h_i = h = 1 for all the spins. 
				% Then I would run the code, generate samples and see what I would get for means and correlations. 
				% If the code is right, you should get m_i = <s_i> = tanh(h) and C_ij = <s_i s_j>-<s_i><s_j> = 0.
				% You can then change h_i to different values while J= 0 And you should get the same result: m_i = tanh(h_i) and C_ij = 0.
				% You can then put J_ij = some constant for some pairs of spins, namely build a network in which spins are connected in pairs (it is called a gas of dimers). 
				% You can analytically solve for m and C (it is a good exercise) and you can also find the solution in many papers (including my PRE paper in 2011 with John Hertz) and then check your simulations results agains the analytical expression. 
				% All these would be sanity checks that there is nothing seriously wrong with the code. 



function JHdiscon = JHs(N,jn,h_on,sparsity,time,T,name)
   %JHstruct = struct('Jgaus',{},'Hfield',{});
 	JHdiscon = struct('Jgaus',{},'Jsparse',{},'Hfield',{},'Hsparse',{});

 		for i = 1:jn
	R3 = zeros(N,N);
        
    	if h_on == 1
		h = ones(1,N);
		JHdiscon(i).Hfield = h;
		JHdiscon(i).Hsparse = h;
		else
		JHdiscon(i).Hfield = zeros(1,N);
		JHdiscon(i).Hsparse = zeros(1,N);
    	end

    	% To sparsify all other H values, supressing for now
    	if i > 1 
    		h = randn(1,N);
			hsparse = h.*(double(rand(1,N)> sparsity));
			JHdiscon(i).Hfield = h;
			JHdiscon(i).Hsparse = hsparse;
		end

	JHdiscon(i).Jgaus = R3;
	JHdiscon(i).Jsparse =  R3;
	save([name,'\',time(1:5),'JHsanity_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'JHdiscon');
	disp(['End JH disconnected run ',num2str(i)]) % Current sate output
	end

	
    %%save(['JHsanity_N',num2str(N),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time,'.mat'],'JHstruct');
    
end