%% J matrix for Ferromagnetic lattice

function JHferr = JF(N,jn,h_on,sparsity,time,T,name)

   %JHstruct = struct('Jgaus',{},'Hfield',{});
 	JHferr = struct('Jgaus',{},'Jsparse',{},'Hfield',{},'Hsparse',{});
	for i = 1:jn

		R = ones(N);
		R3 = triu(R,1)+triu(R,1)'; % Removes diagonal








	% R = double(normrnd(0,1/N,[N,N]));
	% R = R - diag(diag(R));
    % R = R/sqrt(N/2);                              % Normalization, not sure what sorta normalization this should be, in Nicola's code it was dependent on size of N
	% R2 = triu(R,1); %Upper triangle minus diagonal
   	% R3 = R2 + triu(R2)';
	
	% sparsity of couplings
	
	%sparsity = Np_inactive/Np;				
	% Nicola's sparsity method, Gives the fraction of active couplings, so number active out of total parameters
	% I will specify the sparsity amount here until further notice

	%%sparsity = 0.3;		%| arbitrarially setting here, consider removing since putting in at higher level
	%spars = double(rand(N)> sparsity);
	%R2s = triu(R.*spars,1);
    %R3s = R2s + triu(R2s)';			% decimates random values of vector of connection values, not exactly certain how yet. Not sure why using randperm
        
    	if h_on == 1

		% h = randn(1,N);
		% hsparse = h.(*double(rand(1,N)> sparsity));
		% JHferr(i).Hfield = h;
		% JHferr(i).Hsparse = hsparse;
		% else

		JHferr(i).Hfield = zeros(1,N);
		JHferr(i).Hsparse = zeros(1,N);
    	end

	JHferr(i).Jgaus = R3;
	JHferr(i).Jsparse =  R3;


save([name,'\',time(1:5),'JHferr_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'JHferr');	

disp(['End JH ferr run ',num2str(i)]) % Current sate output


end

	
    %%save(['JHstruct_N',num2str(N),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time,'.mat'],'JHstruct');
end