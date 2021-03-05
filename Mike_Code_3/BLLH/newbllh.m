% newbllh.m
% BLLH conversion applet


function OverStruct2 = newbllh(OverStruct)
			%Jstornew = newbllh(OverStruct)

OverStruct2 = OverStruct;
%OverStruct2 = OverStruct;

for  ali = 1:length(OverStruct)           
            
	% Topology
	
	topovec = OverStruct(ali).topologies;
	ltopo = length(topovec);
	
	
	% Sparsity
	sprsvec = OverStruct(ali).sparsity;
	lsprs = length(sprsvec);
	
	
	% Beta
	betavec = OverStruct(ali).betavec;
	lbeta = length(betavec);
	
	% T
	Tvec = OverStruct(ali).Tvec;
	lT = length(Tvec);
	
	% N
	Nvec = OverStruct(ali).Nvec;
	lN = length(Nvec);
	
	% models
	% colorvec = {'b-o','r-o','m-o','g-o','y-o','k-o'};
	% modelvec = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  },{'MDLcount'}];
	%modelvec = modelvec; 
	%lm = length(modelvec);
	
	% Jstornew = struct
	
	%for St = 1:lsprs 
	
		%fn = fieldnames(OverStruct(St).AllStruct);
		fn = fieldnames(OverStruct(ali).AllStruct);

		for ft = 1:length(fn)
	

			for Tt = 1:ltopo 
	
				LLH = OverStruct(ali).AllStruct.(fn{ft}).BLLH(Tt);
				Jcontru = OverStruct(ali).AllStruct.(fn{ft}).Jcontru{Tt};

				N = length(Jcontru);

				%LLH = OverStruct(St).AllStruct.(fn{ft}).BLLH(Tt);
				%Jcontru = OverStruct(St).AllStruct.(fn{ft}).Jcontru{Tt};


				%Total percent connection error%
				LLH.perconerr.BIC      = LLH.totconerr.BIC      ./ N.^2;
				LLH.perconerr.AIC      = LLH.totconerr.AIC      ./ N.^2;
				LLH.perconerr.MDLl     = LLH.totconerr.MDLl     ./ N.^2;
				LLH.perconerr.MDLu     = LLH.totconerr.MDLu     ./ N.^2;
				LLH.perconerr.MDLent   = LLH.totconerr.MDLent   ./ N.^2;
				LLH.perconerr.MDLcount = LLH.totconerr.MDLcount ./ N.^2;



				% True Negatives amount
				LLH.tnconerr.BIC      = sum(double((Jcontru + LLH.Jcon.BIC       ) == 0),'all');
				LLH.tnconerr.AIC      = sum(double((Jcontru + LLH.Jcon.AIC       ) == 0),'all');
				LLH.tnconerr.MDLl     = sum(double((Jcontru + LLH.Jcon.MDLl      ) == 0),'all');
				LLH.tnconerr.MDLu     = sum(double((Jcontru + LLH.Jcon.MDLu      ) == 0),'all');
				LLH.tnconerr.MDLent   = sum(double((Jcontru + LLH.Jcon.MDLent    ) == 0),'all');
				LLH.tnconerr.MDLcount = sum(double((Jcontru + LLH.Jcon.MDLcount  ) == 0),'all');
				
				% True positives amount
				LLH.tpconerr.BIC      = sum(double((Jcontru + LLH.Jcon.BIC       ) == 2),'all');
				LLH.tpconerr.AIC      = sum(double((Jcontru + LLH.Jcon.AIC       ) == 2),'all');
				LLH.tpconerr.MDLl     = sum(double((Jcontru + LLH.Jcon.MDLl      ) == 2),'all');
				LLH.tpconerr.MDLu     = sum(double((Jcontru + LLH.Jcon.MDLu      ) == 2),'all');
				LLH.tpconerr.MDLent   = sum(double((Jcontru + LLH.Jcon.MDLent    ) == 2),'all');
				LLH.tpconerr.MDLcount = sum(double((Jcontru + LLH.Jcon.MDLcount  ) == 2),'all');
					
				LLH.fallout.BIC     	=	LLH.fpconerr.BIC      ./ (LLH.fpconerr.BIC      + LLH.tnconerr.BIC     );
				LLH.fallout.AIC     	=	LLH.fpconerr.AIC      ./ (LLH.fpconerr.AIC      + LLH.tnconerr.AIC     );
				LLH.fallout.MDLl    	=	LLH.fpconerr.MDLl     ./ (LLH.fpconerr.MDLl     + LLH.tnconerr.MDLl    );
				LLH.fallout.MDLu    	=	LLH.fpconerr.MDLu     ./ (LLH.fpconerr.MDLu     + LLH.tnconerr.MDLu    );
				LLH.fallout.MDLent  	=	LLH.fpconerr.MDLent   ./ (LLH.fpconerr.MDLent   + LLH.tnconerr.MDLent  );
				LLH.fallout.MDLcount	=	LLH.fpconerr.MDLcount ./ (LLH.fpconerr.MDLcount + LLH.tnconerr.MDLcount);
			
				%LLH.recall.		=	tpconerr ./ (fnconerr + tpconerr);
				LLH.recall.BIC     	=	LLH.tpconerr.BIC ./ (LLH.fnconerr.BIC + LLH.tpconerr.BIC);
				LLH.recall.AIC     	=	LLH.tpconerr.AIC ./ (LLH.fnconerr.AIC + LLH.tpconerr.AIC);
				LLH.recall.MDLl    	=	LLH.tpconerr.MDLl ./ (LLH.fnconerr.MDLl + LLH.tpconerr.MDLl);
				LLH.recall.MDLu    	=	LLH.tpconerr.MDLu ./ (LLH.fnconerr.MDLu + LLH.tpconerr.MDLu);
				LLH.recall.MDLent  	=	LLH.tpconerr.MDLent ./ (LLH.fnconerr.MDLent + LLH.tpconerr.MDLent);
				LLH.recall.MDLcount	=	LLH.tpconerr.MDLcount ./ (LLH.fnconerr.MDLcount + LLH.tpconerr.MDLcount);
			
				%LLH.precision. 	=	tpconerr ./ (tpconerr + fpconerr);
				LLH.precision.BIC     	=	LLH.tpconerr.BIC      ./ (LLH.fpconerr.BIC      + LLH.tpconerr.BIC     );
				LLH.precision.AIC     	=	LLH.tpconerr.AIC      ./ (LLH.fpconerr.AIC      + LLH.tpconerr.AIC     );
				LLH.precision.MDLl    	=	LLH.tpconerr.MDLl     ./ (LLH.fpconerr.MDLl     + LLH.tpconerr.MDLl    );
				LLH.precision.MDLu    	=	LLH.tpconerr.MDLu     ./ (LLH.fpconerr.MDLu     + LLH.tpconerr.MDLu    );
				LLH.precision.MDLent  	=	LLH.tpconerr.MDLent   ./ (LLH.fpconerr.MDLent   + LLH.tpconerr.MDLent  );
				LLH.precision.MDLcount	=	LLH.tpconerr.MDLcount ./ (LLH.fpconerr.MDLcount + LLH.tpconerr.MDLcount);
	
                
                
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt).perconerr = {};
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt).tnconerr  = {};
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt).tpconerr  = {};
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt).fallout   = {};
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt).recall    = {};
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt).precision = {};

				%OverStruct(St).AllStruct.(fn{ft}).BLLH = LLH;
				OverStruct2(ali).AllStruct.(fn{ft}).BLLH(Tt) = LLH;
			end
		end
	%end



end

Jstor = modelgraphs9(OverStruct2);
OverStruct2(1).AllJstor = Jstor;



%OverStruct2 = OverStruct;