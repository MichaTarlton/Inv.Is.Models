%BLLH2.m
% adds symmetry opti
% T: number of observations, length of S
% N: Number of nodes
% h_on,used for field variaable to include h value


function [LLH,statvecs,stats,jta] = PBLLH4(T,N,tp,beta,sprs,h_on,SStruct,Jtru,Adjset,st)

LLH = struct('NodeModel',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});


modelvec = [{'BIC'},{'AIC'},{'MDLl'},{'MDLu'},{'MDLent'},{'MDLcount'}];
lm = length(modelvec);
statfieldnames  = [{'symrate'},{'asymrate'},{'totconerr'},{'perconerr'},{'tnconerr'},{'tpconerr'},{'fnconerr'},{'fpconerr'},{'fallout'},{'recall'},{'precision'},{'TFS'},{'TFR'},{'symtotconerr'},{'symperconerr'},{'symtotconsererr'},{'symperconsererr'},{'Jrrerr'},{'Jsymrrerr'},{'Jsymconserrrerr'},{'Javgsymrrerr'}];
fm = length(statfieldnames);


%% creating storage struct cause it gets fucky laters
statvecs = struct;
stats = struct; 	%('mean',{},'stddev',{},'stderr',{},'variance',{},'range',{});
for ft1 = 1:length(statfieldnames)
	for mt1 = 1:lm
		statvecs.(modelvec{mt1}).(statfieldnames{ft1}) = [];
	end
end


%for st = 1:size(SStruct,2) % Structs make sizing weird, this is correct though
%	Adjset = JHnorm.Adjset;
%	Jtru = JHnorm.Jtopo;
%	S = SStruct.S_hat;
	%S = S'; %Removed the top level transpose for now
	
	%LLH = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});
	
	%disp(['Topology ', num2str])
	disp([' T ', num2str(T),' N ', num2str(N),'C / Sprs ', num2str(sprs),' bt ', num2str(beta),' Topo ',  num2str(tp),'trl', num2str(st)])
	tic
	NodeModel = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{});	
	
	 JconBIC   	  = [];
	 JconAIC   	  = [];
	 JconMDLl  	  = [];
	 JconMDLu  	  = [];
	 JconMDLent   = [];
	 JconMDLcount = [];
	
	 JreconBIC   	= [];
	 JreconAIC   	= [];
	 JreconMDLl  	= [];
	 JreconMDLu  	= [];
	 JreconMDLent   = [];
	 JreconMDLcount = [];
	
	 thetaBIC      = [];
	 thetaAIC      = [];
	 thetaMDLl     = [];
	 thetaMDLu     = [];
	 thetaMDLent   = [];
	 thetaMDLcount = [];
	
	
		parfor t = 1:N
			
			disp(['Sprs ', num2str(sprs),' Topo ', num2str,' N ', num2str(N),' T ', num2str(T), ' Node ', num2str(t)])
			
			X = S(:,setdiff(1:N,t));
			Y = S(:,t);
			
			[w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,1); % for my purposes h_on is always set to on here in order for the math to function correctly
			%[w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,h_on);
			
			%LLH.NodeModel(t).w_ML 	 = w_ML;
			%LLH.NodeModel(t).l_ML		 = l_ML;
			%LLH.NodeModel(t).posterior = posterior;
			%LLH.NodeModel(t).cost		 = cost;
			%LLH.NodeModel(t).BestModel = BestModel;
			%LLH.NodeModel(t).IMAX		 = IMAX;
	        %disp('check 1')
			NodeModel(t).w_ML 	 = w_ML;
			NodeModel(t).l_ML		 = l_ML;
			NodeModel(t).posterior = posterior;
			NodeModel(t).cost		 = cost;
			NodeModel(t).BestModel = BestModel;
			NodeModel(t).IMAX		 = IMAX;
			%disp('check 2')
	        JconBIC(t,:)    = BestModel.BIC;
	        JconAIC(t,:)    = BestModel.AIC;
	        JconMDLl(t,:)   = BestModel.MDLl;
	        JconMDLu(t,:)   = BestModel.MDLu;
	        JconMDLent(t,:)   = BestModel.MDLent; % making the connection graph at an earlier juncture
	        JconMDLcount(t,:)  = BestModel.MDLcount;
			%disp('check 3')        
	        thetaBIC(t,:)    = w_ML(IMAX.BIC,:);
			thetaAIC(t,:)    = w_ML(IMAX.AIC,:);
			thetaMDLl(t,:)   = w_ML(IMAX.MDLl,:);
			thetaMDLu(t,:)   = w_ML(IMAX.MDLu,:);
			thetaMDLent(t,:) = w_ML(IMAX.MDLent,:); % this also works I think
			thetaMDLcount(t,:) = w_ML(IMAX.MDLcount,:);
	        %theta(t,:)       = LLH(t).w_ML(LLH(t).IMAX.MDLent,:); %IMAX is index of best model
	
			
			%LLH.JHrecon = [;LLH(t).theta]; 
		end
		toc
		
		%disp(['Run ',num2str(jta),' / ',num2str(jtatot)])
		%jta = jta + 1;

	LLH.NodeModel = NodeModel;
	
	% Kinda forgot what this is, check the paper and above
	% oh is this the pre seperation of j and h
	% all of this could be put in the below at some more optimized time

	%LLH.theta 	 =  theta;
	LLH.theta.BIC      = thetaBIC;
	LLH.theta.AIC      = thetaAIC;
	LLH.theta.MDLl     = thetaMDLl;
	LLH.theta.MDLu     = thetaMDLu;
	LLH.theta.MDLent   = thetaMDLent;
	LLH.theta.MDLcount = thetaMDLcount;
	
	% h reconstruction values?
	LLH.hrecon.BIC        =  thetaBIC(:,1);
	LLH.hrecon.AIC        =  thetaAIC(:,1);
	LLH.hrecon.MDLl       =  thetaMDLl(:,1);
	LLH.hrecon.MDLu       =  thetaMDLu(:,1);
	LLH.hrecon.MDLent     =  thetaMDLent(:,1);
	LLH.hrecon.MDLcount   =  thetaMDLcount(:,1);
	%LLH.hrecon   =  theta(:,1);

	% J connectivity graph (seperating it from h value column)
	JconBIC = JconBIC(:,2:end);
	JconAIC = JconAIC(:,2:end);
	JconMDLl = JconMDLl(:,2:end);
	JconMDLu = JconMDLu(:,2:end);
	JconMDLent = JconMDLent(:,2:end);
	JconMDLcount = JconMDLcount(:,2:end);

	% adding the zero diagonal and storing in the struct	
	LLH.Jcon.BIC      = [zeros(N,1), triu(JconBIC)] +      [tril(JconBIC,-1),zeros(N,1)];
	LLH.Jcon.AIC      = [zeros(N,1), triu(JconAIC)] +      [tril(JconAIC,-1),zeros(N,1)];
	LLH.Jcon.MDLl     = [zeros(N,1), triu(JconMDLl)] +     [tril(JconMDLl,-1),zeros(N,1)];
	LLH.Jcon.MDLu     = [zeros(N,1), triu(JconMDLu)] +     [tril(JconMDLu,-1),zeros(N,1)];
	LLH.Jcon.MDLent   = [zeros(N,1), triu(JconMDLent)] +   [tril(JconMDLent,-1),zeros(N,1)];
	LLH.Jcon.MDLcount = [zeros(N,1), triu(JconMDLcount)] + [tril(JconMDLcount,-1),zeros(N,1)];
	
	% reconstructed J values
	JreconBIC   =  thetaBIC(:,2:end);  
	JreconAIC   =  thetaAIC(:,2:end);  
	JreconMDLl   =  thetaMDLl(:,2:end);  
	JreconMDLu   =  thetaMDLu(:,2:end);  
	JreconMDLent   =  thetaMDLent(:,2:end);  
	JreconMDLcount   =  thetaMDLcount(:,2:end);  
	%Jrecon   =  theta(:,2:end); 

	% adding zeros to the diag
	LLH.Jrecon.BIC      = [zeros(N,1), triu(JreconBIC)] + [tril(JreconBIC,-1),zeros(N,1)];
	LLH.Jrecon.AIC      = [zeros(N,1), triu(JreconAIC)] + [tril(JreconAIC,-1),zeros(N,1)];
	LLH.Jrecon.MDLl     = [zeros(N,1), triu(JreconMDLl)] + [tril(JreconMDLl,-1),zeros(N,1)];
	LLH.Jrecon.MDLu     = [zeros(N,1), triu(JreconMDLu)] + [tril(JreconMDLu,-1),zeros(N,1)];
	LLH.Jrecon.MDLent   = [zeros(N,1), triu(JreconMDLent)] + [tril(JreconMDLent,-1),zeros(N,1)];
	LLH.Jrecon.MDLcount = [zeros(N,1), triu(JreconMDLcount)] + [tril(JreconMDLcount,-1),zeros(N,1)];
	%LLH.Jrecon = [zeros(N,1), triu(Jrecon)] + [tril(Jrecon,-1),zeros(N,1)]; 
	
	for mt = 1:lm
		
		model = modelvec{mt};

		%trying to find which ones are connected inf v tru
		% Jcon = double(not(LLH(1).Jrecon == 0)) % find the connected ones (couldalso jsut grab the best models vector)
		%Jcon = Jcon(:,2:end); % Yeah this doesn'twork because you haven'tadded the zero diag yet
		%Jcon = [zeros(N,1), triu(Jcon)] + [tril(Jcon,-1),zeros(N,1)];
		%LLH.Jcon = Jcon;
		
		% the first matrix manip to find the sym and asym
		% probably don't need to store
		LLH.Jasym.(model)      = tril(LLH.Jcon.(model))' + triu(LLH.Jcon.(model));
		Jasym = LLH.Jasym.(model);
		Jrecon = LLH.Jrecon.(model);
		
		% symmetry rate
		LLH.symrate.(model)	= sum(double(LLH.Jasym.(model) == 2),'all');

		% asymmetry rate
		LLH.asymrate.(model)      = sum(double(LLH.Jasym.(model)      == 1),'all');

		% For finding connection error rates 
		% Perhaps should incorporate back into BLLH
		% totconerr: total error in infer connections
		% fpconerr : False positives
		% fnconerr: False Negatives
		
		% Total amount of error
		LLH.totconerr.(model)      = sum(double(not((LLH.Jcon.(model)      - Adjset) == 0)),'all');
		
		%Total percent connection error%
		LLH.perconerr.(model)      = LLH.totconerr.(model)      ./ N.^2;
		
		% True Negatives amount
		LLH.tnconerr.(model)      = sum(double((Adjset + LLH.Jcon.(model)       ) == 0),'all');
		
		% True positives amount
		LLH.tpconerr.(model)      = sum(double((Adjset + LLH.Jcon.(model)       ) == 2),'all');
		
		% False Negatives amount
		LLH.fnconerr.(model)      = sum(double((Adjset - LLH.Jcon.(model)       ) == 1),'all');
		
		% False positives amount
		LLH.fpconerr.(model)      = sum(double((Adjset - LLH.Jcon.(model)       ) == -1),'all');

		% True Negatives percent
		LLH.tnconper.(model)      = (sum(double((Adjset + LLH.Jcon.(model)       ) == 0),'all')) ./ N.^2;
		
		% True positives percent
		LLH.tpconper.(model)      = (sum(double((Adjset + LLH.Jcon.(model)       ) == 2),'all')) ./ N.^2;
		
		% False Negatives percent
		LLH.fnconper.(model)      = (sum(double((Adjset - LLH.Jcon.(model)       ) == 1),'all')) ./ N.^2;
		
		% False positives percent
		LLH.fpconper.(model)      = (sum(double((Adjset - LLH.Jcon.(model)       ) == -1),'all')) ./ N.^2;

		LLH.fallout.(model)     	=	LLH.fpconerr.(model)      ./ (LLH.fpconerr.(model)      + LLH.tnconerr.(model)     );			
		
		%LLH.recall.		=	tpconerr ./ (fnconerr + tpconerr);
		LLH.recall.(model)     	=	LLH.tpconerr.(model) ./ (LLH.fnconerr.(model) + LLH.tpconerr.(model));		
	
		%LLH.precision. 	=	tpconerr ./ (tpconerr + fpconerr);
		LLH.precision.(model)     	=	LLH.tpconerr.(model)      ./ (LLH.fpconerr.(model)      + LLH.tpconerr.(model)     );

		%% True False Sum			
		LLH.TFS.(model)     	=	LLH.tpconerr.(model) + LLH.fpconerr.(model);

		%%True False Ratio
		LLH.TFR.(model)     	=	(LLH.tpconerr.(model) - LLH.fpconerr.(model)) ./ (LLH.TFS.(model) + LLH.tnconerr.(model));	


		%smmetry optimization
		Jcon = LLH.Jcon.(model);
		[row,col] = find((((Jasym + Jasym') ~= 0) + Jcon) == 1);

		LLH.asymrow.(model) = row;
		LLH.asymcol.(model) = col;


		Jsymrecon = zeros(size(Jrecon));
		
		for coordno = 1:length(row)
			
			ci = row(coordno);

			wml = LLH.NodeModel(row(coordno)).w_ML(1:N,2:end) ;
			%wml2 = [zeros(N,1), triu(wml)] +      [tril(wml,-1),zeros(N,1)]; This is incorrect, do not use

			idx0 = 1:1:size(wml,2);
			zeron = zeros(size(wml,1),1);
			A = setdiff(idx0,ci:N);
			B = setdiff(idx0,1:ci-1);

			wml2 = [wml(:,A),zeron,wml(:,B)];
 

			v = nonzeros(wml2(:,col(coordno)));
			reconval = v(end);

			Jsymrecon(row(coordno),col(coordno)) = reconval;
		end

		%LLH.Jsymcon.(model) = (Jsymrecon + Jrecon) ~= 0 ;
		%LLH.Jsymrecon.(model) = Jsymrecon + Jrecon;
		Jsymrecon = Jsymrecon + Jrecon;
		LLH.Jsymrecon.(model) = Jsymrecon;
		LLH.Jsymcon.(model) = Jsymrecon ~= 0 ;

		Jsymconser = (Jasym + Jasym') == 2 ; % make a conservative version
		Jsymreconser = Jsymconser .* Jsymrecon;
		LLH.Jsymreconser.(model) 	= Jsymreconser;
		LLH.Jsymconser.(model)		= Jsymconser;

		avgcon = (triu(Jsymrecon) + tril(Jsymrecon)')/2;
		Javgsymrecon = avgcon + avgcon';

		symtotconerr = sum(double(not((Javgsymrecon - Jtru) == 0)),'all');
		symtotconsererr = sum(double(not((Jsymreconser - Jtru) == 0)),'all');
		symperconerr = symtotconerr ./ N.^2;
		symperconsererr = symtotconsererr ./ N.^2;
		
		LLH.Javgsymrecon.(model) = Javgsymrecon;
		LLH.symtotconerr.(model)      = symtotconerr;
		LLH.symperconerr.(model)      = symperconerr;
		LLH.symtotconsererr.(model)      = symtotconsererr;
		LLH.symperconsererr.(model)      = symperconsererr;

		LLH.Jrrerr.(model) 				=sqrt(sum((Jrecon 			- Jtru).^2 / sum(Jtru.^2)));
		LLH.Jsymrrerr.(model) 			=sqrt(sum((Jsymrecon 		- Jtru).^2 / sum(Jtru.^2)));
		LLH.Jsymconserrrerr.(model) 	=sqrt(sum((Jsymreconser 	- Jtru).^2 / sum(Jtru.^2)));
		LLH.Javgsymrrerr.(model)		=sqrt(sum((Javgsymrecon 	- Jtru).^2 / sum(Jtru.^2)));


		%% putting inside st for each llh struct loop
		
		for ft = 1:length(statfieldnames)
			statvecs.(model).(statfieldnames{ft}) = [statvecs.(model).(statfieldnames{ft}), LLH.(statfieldnames{ft}).(model)];
		end
	end
end


for mt2 = 1:lm 
	model = modelvec{mt2};
	for ft2 = 1:fm
		stats.mean.(model).(statfieldnames{ft2}) =		[mean(statvecs.(model).(statfieldnames{ft2}))];
		stats.stddev.(model).(statfieldnames{ft2}) =	[std(statvecs.(model).(statfieldnames{ft2}))];
		%stats.stderr.(model).(statfieldnames{ft2}) =	[std(statvecs.(model).(statfieldnames{ft2}))/sqrt(size(SStruct,2))]; % std(x)/sqrt(length(x));
		stats.stderr.(model).(statfieldnames{ft2}) =	[std(statvecs.(model).(statfieldnames{ft2}))/size(SStruct,2)]; % should be this since the N in the STDDEV is actually Nodes^2 for size of matrix
		stats.variance.(model).(statfieldnames{ft2}) =	[var(statvecs.(model).(statfieldnames{ft2}))];
		stats.range.(model).(statfieldnames{ft2}) =		[range(statvecs.(model).(statfieldnames{ft2}),'all')];
	end
end		

%% Making the stats for the LLH over trials
	%{'NodeModel'   }
    %,{'theta'       }
    %,{'hrecon'      }
    %,{'Jrecon'      }
    %,{'Jcon'        }
    %,{'Jasym'       }
    %,{'asymrow'     }
    %,{'asymcol'     }
    %,{'Jsymrecon'   }
    %,{'Jsymcon'     }
    %,{'Javgsymrecon'}


    %% inside the 
    % modelvec = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  },{'MDLcount'}]
    % This is above already


%ok luckily though at a glance the values seem similar,so I may be able to just avg them out

% Ok now for comparing original connectivity to recon
% this will have to be done at a higher level, but writing here to keep track
% so the problm is that there are no zero values in my original values unless I set it to sparsify or disconnect values above a certain threshold
%
% Jcon = LLH(1).Jcon
% Jcon = double(not(Jcon == 0))
% 
