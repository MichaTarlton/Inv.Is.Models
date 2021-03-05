%BLLH2.m
% adds symmetry opti
% T: number of observations, length of S
% N: Number of nodes
% h_on,used for field variaable to include h value


function [LLH,statvecs,stats,jta] = PBLLH2(T,N,beta,sprs,h_on,SStruct,JHnorm,jta,jtatot)

LLH = struct('NodeModel',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});


modelvec = [{'BIC'},{'AIC'},{'MDLl'},{'MDLu'},{'MDLent'},{'MDLcount'}];
lm = length(modelvec);
statfieldnames  = [{'symrate'},{'asymrate'},{'totconerr'},{'perconerr'},{'tnconerr'},{'tpconerr'},{'fnconerr'},{'fpconerr'},{'fallout'},{'recall'},{'precision'},{'TFS'},{'TFR'},{'symtotconerr'},{'symperconerr'},{'Jrrerr'},{'Jsymrrerr'},{'Javgsymrrerr'}];
fm = length(statfieldnames);


%% creating storage struct cause it gets fucky laters
statvecs = struct;
stats = struct; 	%('mean',{},'stddev',{},'stderr',{},'variance',{},'range',{});
for ft1 = 1:length(statfieldnames)
	for mt1 = 1:lm
		statvecs.(modelvec{mt1}).(statfieldnames{ft1}) = [];
	end
end


for st = 1:size(SStruct,2) % Structs make sizing weird, this is correct though
	Adjset = JHnorm(st).Adjset;
	Jtru = JHnorm(st).Jtopo;
	S = SStruct(st).S_hat;
	%S = S'; %Removed the top level transpose for now
	
	%LLH = struct('w_ML',{},'l_ML',{},'posterior',{},'cost',{},'BestModel',{},'IMAX',{},'theta',{},'hrecon',{},'Jrecon',{},'Jcon',{},'Jasym',{},'symrate',{},'asymrate',{});
	
	%disp(['Topology ', num2str(st)])
	disp(['Sprs ', num2str(sprs),' Topo ', num2str(st),' bt ', num2str(beta),' N ', num2str(N),' T ', num2str(T)])
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
			
			disp(['Sprs ', num2str(sprs),' Topo ', num2str(st),' N ', num2str(N),' T ', num2str(T), ' Node ', num2str(t)])
			
			X = S(:,setdiff(1:N,t));
			Y = S(:,t);
			
			[w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,1); % for my purposes h_on is always set to on here in order for the math to function correctly
			%[w_ML,l_ML,posterior,cost,BestModel,IMAX] = decimation_logistic_Model_Selection(X,Y,h_on);
			
			%LLH(st).NodeModel(t).w_ML 	 = w_ML;
			%LLH(st).NodeModel(t).l_ML		 = l_ML;
			%LLH(st).NodeModel(t).posterior = posterior;
			%LLH(st).NodeModel(t).cost		 = cost;
			%LLH(st).NodeModel(t).BestModel = BestModel;
			%LLH(st).NodeModel(t).IMAX		 = IMAX;
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
		
		disp(['Run ',num2str(jta),' / ',num2str(jtatot)])
		jta = jta + 1;

	LLH(st).NodeModel = NodeModel;
	
	% Kinda forgot what this is, check the paper and above
	% oh is this the pre seperation of j and h
	% all of this could be put in the below at some more optimized time

	%LLH(st).theta 	 =  theta;
	LLH(st).theta.BIC      = thetaBIC;
	LLH(st).theta.AIC      = thetaAIC;
	LLH(st).theta.MDLl     = thetaMDLl;
	LLH(st).theta.MDLu     = thetaMDLu;
	LLH(st).theta.MDLent   = thetaMDLent;
	LLH(st).theta.MDLcount = thetaMDLcount;
	
	% h reconstruction values?
	LLH(st).hrecon.BIC        =  thetaBIC(:,1);
	LLH(st).hrecon.AIC        =  thetaAIC(:,1);
	LLH(st).hrecon.MDLl       =  thetaMDLl(:,1);
	LLH(st).hrecon.MDLu       =  thetaMDLu(:,1);
	LLH(st).hrecon.MDLent     =  thetaMDLent(:,1);
	LLH(st).hrecon.MDLcount   =  thetaMDLcount(:,1);
	%LLH(st).hrecon   =  theta(:,1);

	% J connectivity graph (seperating it from h value column)
	JconBIC = JconBIC(:,2:end);
	JconAIC = JconAIC(:,2:end);
	JconMDLl = JconMDLl(:,2:end);
	JconMDLu = JconMDLu(:,2:end);
	JconMDLent = JconMDLent(:,2:end);
	JconMDLcount = JconMDLcount(:,2:end);

	% adding the zero diagonal and storing in the struct	
	LLH(st).Jcon.BIC      = [zeros(N,1), triu(JconBIC)] +      [tril(JconBIC,-1),zeros(N,1)];
	LLH(st).Jcon.AIC      = [zeros(N,1), triu(JconAIC)] +      [tril(JconAIC,-1),zeros(N,1)];
	LLH(st).Jcon.MDLl     = [zeros(N,1), triu(JconMDLl)] +     [tril(JconMDLl,-1),zeros(N,1)];
	LLH(st).Jcon.MDLu     = [zeros(N,1), triu(JconMDLu)] +     [tril(JconMDLu,-1),zeros(N,1)];
	LLH(st).Jcon.MDLent   = [zeros(N,1), triu(JconMDLent)] +   [tril(JconMDLent,-1),zeros(N,1)];
	LLH(st).Jcon.MDLcount = [zeros(N,1), triu(JconMDLcount)] + [tril(JconMDLcount,-1),zeros(N,1)];
	
	% reconstructed J values
	JreconBIC   =  thetaBIC(:,2:end);  
	JreconAIC   =  thetaAIC(:,2:end);  
	JreconMDLl   =  thetaMDLl(:,2:end);  
	JreconMDLu   =  thetaMDLu(:,2:end);  
	JreconMDLent   =  thetaMDLent(:,2:end);  
	JreconMDLcount   =  thetaMDLcount(:,2:end);  
	%Jrecon   =  theta(:,2:end); 

	% adding zeros to the diag
	LLH(st).Jrecon.BIC      = [zeros(N,1), triu(JreconBIC)] + [tril(JreconBIC,-1),zeros(N,1)];
	LLH(st).Jrecon.AIC      = [zeros(N,1), triu(JreconAIC)] + [tril(JreconAIC,-1),zeros(N,1)];
	LLH(st).Jrecon.MDLl     = [zeros(N,1), triu(JreconMDLl)] + [tril(JreconMDLl,-1),zeros(N,1)];
	LLH(st).Jrecon.MDLu     = [zeros(N,1), triu(JreconMDLu)] + [tril(JreconMDLu,-1),zeros(N,1)];
	LLH(st).Jrecon.MDLent   = [zeros(N,1), triu(JreconMDLent)] + [tril(JreconMDLent,-1),zeros(N,1)];
	LLH(st).Jrecon.MDLcount = [zeros(N,1), triu(JreconMDLcount)] + [tril(JreconMDLcount,-1),zeros(N,1)];
	%LLH(st).Jrecon = [zeros(N,1), triu(Jrecon)] + [tril(Jrecon,-1),zeros(N,1)]; 
	
	for mt = 1:lm
		
		model = modelvec{mt};

		%trying to find which ones are connected inf v tru
		% Jcon = double(not(LLH(1).Jrecon == 0)) % find the connected ones (couldalso jsut grab the best models vector)
		%Jcon = Jcon(:,2:end); % Yeah this doesn'twork because you haven'tadded the zero diag yet
		%Jcon = [zeros(N,1), triu(Jcon)] + [tril(Jcon,-1),zeros(N,1)];
		%LLH(st).Jcon = Jcon;
		
		% the first matrix manip to find the sym and asym
		% probably don't need to store
		LLH(st).Jasym.(model)      = tril(LLH(st).Jcon.(model))' + triu(LLH(st).Jcon.(model));
		Jasym = LLH(st).Jasym.(model);
		Jrecon = LLH(st).Jrecon.(model);
		
		% symmetry rate
		LLH(st).symrate.(model)	= sum(double(LLH(st).Jasym.(model) == 2),'all');

		% asymmetry rate
		LLH(st).asymrate.(model)      = sum(double(LLH(st).Jasym.(model)      == 1),'all');

		% For finding connection error rates 
		% Perhaps should incorporate back into BLLH
		% totconerr: total error in infer connections
		% fpconerr : False positives
		% fnconerr: False Negatives
		
		% Total amount of error
		LLH(st).totconerr.(model)      = sum(double(not((LLH(st).Jcon.(model)      - Adjset) == 0)),'all');
		
		%Total percent connection error%
		LLH(st).perconerr.(model)      = LLH(st).totconerr.(model)      ./ N.^2;
		
		% True Negatives amount
		LLH(st).tnconerr.(model)      = sum(double((Adjset + LLH(st).Jcon.(model)       ) == 0),'all');
		
		% True positives amount
		LLH(st).tpconerr.(model)      = sum(double((Adjset + LLH(st).Jcon.(model)       ) == 2),'all');
		
		% False Negatives amount
		LLH(st).fnconerr.(model)      = sum(double((Adjset - LLH(st).Jcon.(model)       ) == 1),'all');
		
		% False positives amount
		LLH(st).fpconerr.(model)      = sum(double((Adjset - LLH(st).Jcon.(model)       ) == -1),'all');

		LLH(st).fallout.(model)     	=	LLH(st).fpconerr.(model)      ./ (LLH(st).fpconerr.(model)      + LLH(st).tnconerr.(model)     );			
		
		%LLH(st).recall.		=	tpconerr ./ (fnconerr + tpconerr);
		LLH(st).recall.(model)     	=	LLH(st).tpconerr.(model) ./ (LLH(st).fnconerr.(model) + LLH(st).tpconerr.(model));		
	
		%LLH(st).precision. 	=	tpconerr ./ (tpconerr + fpconerr);
		LLH(st).precision.(model)     	=	LLH(st).tpconerr.(model)      ./ (LLH(st).fpconerr.(model)      + LLH(st).tpconerr.(model)     );

		%% True False Sum			
		LLH(st).TFS.(model)     	=	LLH(st).tpconerr.(model) + LLH(st).fpconerr.(model);

		%%True False Ratio
		LLH(st).TFR.(model)     	=	(LLH(st).tpconerr.(model) - LLH(st).fpconerr.(model)) ./ (LLH(st).TFS.(model) + LLH(st).tnconerr.(model));	


		%smmetry optimization
		Jcon = LLH(st).Jcon.(model);
		[row,col] = find((((Jasym + Jasym') ~= 0) + Jcon) == 1);

		LLH(st).asymrow.(model) = row;
		LLH(st).asymcol.(model) = col;


		Jsymrecon = zeros(size(Jrecon));
		for coordno = 1:length(row)
			
			ci = row(coordno);

			wml = LLH(st).NodeModel(row(coordno)).w_ML(1:N,2:end) ;
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

		LLH(st).Jsymrecon.(model) = Jsymrecon + Jrecon;
		LLH(st).Jsymcon.(model) = (Jsymrecon + Jrecon) ~= 0 ;

		avgcon = (triu(Jsymrecon) + tril(Jsymrecon)')/2;
		Javgsymrecon = avgcon + avgcon';

		symtotconerr = sum(double(not((Javgsymrecon - Jtru) == 0)),'all');
		symperconerr = symtotconerr ./ N.^2;
		
		LLH(st).Javgsymrecon.(model) = Javgsymrecon;
		LLH(st).symtotconerr.(model)      = symtotconerr;
		LLH(st).symperconerr.(model)      = symperconerr;


		LLH(st).Jrrerr.(model) 		=sqrt(sum((Jrecon 		- Jtru).^2 / sum(Jtru.^2)));
		LLH(st).Jsymrrerr.(model) 	=sqrt(sum((Jsymrecon 		- Jtru).^2 / sum(Jtru.^2)));
		LLH(st).Javgsymrrerr.(model) =sqrt(sum((Javgsymrecon 	- Jtru).^2 / sum(Jtru.^2)));


		%% putting inside st for each llh struct loop
		
		for ft = 1:length(statfieldnames)
			statvecs.(model).(statfieldnames{ft}) = [statvecs.(model).(statfieldnames{ft}), LLH(st).(statfieldnames{ft}).(model)];
		end
	end
end


for mt2 = 1:lm 
	model = modelvec{mt2};
	for ft2 = 1:fm
		stats.mean.(model).(statfieldnames{ft2}) =		[mean(statvecs.(model).(statfieldnames{ft2}))];
		stats.stddev.(model).(statfieldnames{ft2}) =	[std(statvecs.(model).(statfieldnames{ft2}))];
		stats.stderr.(model).(statfieldnames{ft2}) =	[std(statvecs.(model).(statfieldnames{ft2}))/sqrt(size(SStruct,2))]; % std(x)/sqrt(length(x));
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
