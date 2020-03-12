%% J matrix for gas of dimers

function JHdimer = JD(N,jn,h_on,sparsity,time,T,name)

   %JHstruct = struct('Jgaus',{},'Hfield',{});
 	JHdimer = struct('Jgaus',{},'Jsparse',{},'Hfield',{},'Hsparse',{});
	for i = 1:jn

		R = eye(N);
		%R = eye(N).*double(normrnd(0,1/N,[N,N]));
        R = eye(N).*double(normrnd(0,1/nthroot(N,3),[N,N]));
		R = R/sqrt(N/2); 
		R2 = R(:,randperm(N));
		R3 = triu(R2,1)+triu(R2,1)'; % Removes diagonal








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
		% JHdimer(i).Hfield = h;
		% JHdimer(i).Hsparse = hsparse;
		% else

		JHdimer(i).Hfield = zeros(1,N);
		JHdimer(i).Hsparse = zeros(1,N);
    	end

	JHdimer(i).Jgaus = R3;
	JHdimer(i).Jsparse =  R3;


save([name,'\',time(1:5),'JHDstruct_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'JHdimer');	

disp(['End JH dimer run ',num2str(i)]) % Current sate output


end

	
    %%save(['JHstruct_N',num2str(N),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time,'.mat'],'JHstruct');
end