function Jstor = modelgraphs2(OverStruct,sparsity,betavec,Nvec,Tvec,topologies)


Jstor = struct;
models = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  },{'MDLcount'}];
%topologies = {'sk',1,2,3,4,5,6}


for topology = 1:length(topologies) % for each complete figure
	%h = h = figure;
	%hold on

	for bt = 1:length(betavec)

		figwidth = length(betavec);
		
		for Nt = 1:length(Nvec)
			
			N = Nvec(Nt);

			for Tt = 1:length(Tvec)
				
				T = Tvec(Tt);

				figheight = length(Tvec);
				%name = ['N',num2str(Nvec(Nt)),'T1E',num2str(log10(Tvec(Tt)))];
				
				BICespc			=	[];
				AICespc			=	[];
				MDLlespc		=	[];
				MDLuespc		=	[];
				MDLentespc  	=	[];
				MDLcountespc	=	[];
				BICesfp			=	[];
				AICesfp			=	[];
				MDLlesfp		=	[];
				MDLuesfp		=	[];
				MDLentesfp  	=	[];
				MDLcountesfp	=	[];
				BICesfn			=	[];
				AICesfn			=	[];
				MDLlesfn		=	[];
				MDLuesfn		=	[];
				MDLentesfn  	=	[];
				MDLcountesfn	=	[];
				BICestp			=	[];
				AICestp			=  []; 
				MDLlestp		=  []; 
				MDLuestp		=  []; 
				MDLentestp  	=  []; 
				MDLcountestp	=  []; 
				BICestn			=  []; 
				AICestn			=  []; 
				MDLlestn		=  []; 
				MDLuestn		=  []; 
				MDLentestn  	=  []; 
				MDLcountestn	=  []; 
				BICestc			=	[];
				AICestc			=	[];
				MDLlestc		=	[];
				MDLuestc		=	[];
				MDLentestc  	=	[];
				MDLcountestc	=	[];
				BICfallout		=	[];
				AICfallout			= [];
				MDLlfallout			= [];
				MDLufallout			= [];
				MDLentfallout		= [];
				MDLcountfallout		= [];
				BICrecall			= [];
				AICrecall			= [];
				MDLlrecall			= [];
				MDLurecall			= [];
				MDLentrecall		= [];
				MDLcountrecall		= [];
				BICprecision		= [];
				AICprecision		= [];
				MDLlprecision		= [];
				MDLuprecision		= [];
				MDLentprecision		= [];
				MDLcountprecision	= [];


				for sprs = 1:length(sparsity)
					
					name = ['St',num2str(sprs),'Bt',num2str(bt),'N',num2str(N),'T1E',num2str(log10(T))];

					Jstor(topology).Tt					= topologies(topology);
					Jstor(topology).Topo(bt).BT				= betavec(bt);
					Jstor(topology).Topo(bt).beta(Nt).Nt 	=  N;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).Tt			= T;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Sparsity = sparsity(sprs);
					

					
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.couplings		=  	OverStruct(sprs).AllStruct.(name).couplings;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.c				=  	OverStruct(sprs).AllStruct.(name).c;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.Jcontru		=  	OverStruct(sprs).AllStruct.(name).Jcontru(topology);
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.Jtru			=  	OverStruct(sprs).AllStruct.(name).Jtru(topology);
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.htru			=  	OverStruct(sprs).AllStruct.(name).htru(1,:); %OverStruct(sprs).AllStruct.(name).htru(topology,:); % This doesn't work for as long as we only have one h
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.S				=  	OverStruct(sprs).AllStruct.(name).S(topology);	
					
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.Jmf	=	AllStruct.(name).Jmf;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.hmf	=	AllStruct.(name).hmf;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.Jtap	=	AllStruct.(name).Jtap;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.htap	=	AllStruct.(name).htap;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.Jplmf	=	AllStruct.(name).Jplmf;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.hplmf	=	AllStruct.(name).hplmf;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.Jpllh	=	AllStruct.(name).Jpllh;
					%Jstor.Topo(topology).beta(bt).N(N).T(Tt).Spars(sprs).VarMeth.hpllh	=	AllStruct.(name).hpllh;


					for mt = 1:length(models)

						disp(['Topo',num2str(topologies{topology}),'_Bt',num2str(betavec(bt)),'_N',num2str(Nvec(Nt)),'_T1E',num2str(log10(Tvec(Tt))),'_Sprs',num2str(sparsity(sprs)),'_model_',num2str(models{mt})]);

						model = models{mt};

						Jcontru			 =		OverStruct(sprs).AllStruct.(name).Jcontru{topology};
						Jrecon		 	 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).Jrecon.(model);
						Jcon 			 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).Jcon.(model);
						hrecon			 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).hrecon.(model);


						% For fixing BLLH that's already been done without the true pos/neg values
						% you should remove this on subsequent runs
						
						% This claculation has since been moved into bllh 

						% True Negatives amount
						%tnconerr     = sum(double((Jcontru + Jcon ) == 0),'all');
						% True positives amount
						%tpconerr     = sum(double((Jcontru + Jcon ) == 2),'all');
						


						
						perconerr		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).perconerr.(model);
						tpconerr		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).tpconerr.(model);
						tnconerr		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).tnconerr.(model);
						fpconerr		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).fpconerr.(model);
						fnconerr		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).fnconerr.(model);
						totconerr		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).totconerr.(model);
						symrate			 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).symrate.(model);
						asymrate		 =  	OverStruct(sprs).AllStruct.(name).BLLH(topology).asymrate.(model);
						fallout			=  	OverStruct(sprs).AllStruct.(name).BLLH(topology).fallout.(model);
						recall			=  	OverStruct(sprs).AllStruct.(name).BLLH(topology).recall.(model);
						precision		=  	OverStruct(sprs).AllStruct.(name).BLLH(topology).precision.(model);

						% Moved calculation to BLLH and put the storageinline above above
						%fallout		=	fpconerr ./ (fpconerr + tnconerr);
						%recall		=	tpconerr ./ (fnconerr + tpconerr);
						%precision 	=	tpconerr ./ (tpconerr + fpconerr);	

						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).Jrecon	=  	Jrecon	;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).Jcon 		=  	Jcon 	;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).hrecon	=  	hrecon	;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).perconerr	=  	perconerr;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).tpconerr	=  	tpconerr;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).tnconerr	=  	tnconerr;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).fpconerr	=  	fpconerr;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).fnconerr	=  	fnconerr;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).fallout	=  	fallout	;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).recall	=  	recall	;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).precision =  	precision;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).totconerr	=  	totconerr;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).symrate	=  	symrate;
						Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.(model).asymrate	=  	asymrate;
					end

					BICespc			=	[BICespc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.perconerr];		 
					AICespc			=	[AICespc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.perconerr];		
					MDLlespc		=	[MDLlespc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.perconerr];		
					MDLuespc		=	[MDLuespc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.perconerr];		
					MDLentespc  	=	[MDLentespc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.perconerr];	
					MDLcountespc	=	[MDLcountespc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.perconerr];
					
					BICesfp			=	[BICesfp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.fpconerr];		 
					AICesfp			=	[AICesfp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.fpconerr];		
					MDLlesfp		=	[MDLlesfp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.fpconerr];		
					MDLuesfp		=	[MDLuesfp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.fpconerr];		
					MDLentesfp  	=	[MDLentesfp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.fpconerr];	
					MDLcountesfp	=	[MDLcountesfp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.fpconerr];
					
					BICesfn			=	[BICesfn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.fnconerr];		 
					AICesfn			=	[AICesfn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.fnconerr];		
					MDLlesfn		=	[MDLlesfn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.fnconerr];		
					MDLuesfn		=	[MDLuesfn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.fnconerr];		
					MDLentesfn  	=	[MDLentesfn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.fnconerr];	
					MDLcountesfn	=	[MDLcountesfn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.fnconerr];
					

					BICestp			=	[BICestp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.tpconerr];		 
					AICestp			=	[AICestp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.tpconerr];		
					MDLlestp		=	[MDLlestp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.tpconerr];		
					MDLuestp		=	[MDLuestp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.tpconerr];		
					MDLentestp  	=	[MDLentestp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.tpconerr];	
					MDLcountestp	=	[MDLcountestp, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.tpconerr];
					

					BICestn			=	[BICestn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.tnconerr];		 
					AICestn			=	[AICestn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.tnconerr];		
					MDLlestn		=	[MDLlestn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.tnconerr];		
					MDLuestn		=	[MDLuestn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.tnconerr];		
					MDLentestn  	=	[MDLentestn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.tnconerr];	
					MDLcountestn	=	[MDLcountestn, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.tnconerr];
					

					BICestc			=	[BICestc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.totconerr];		 
					AICestc			=	[AICestc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.totconerr];		
					MDLlestc		=	[MDLlestc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.totconerr];		
					MDLuestc		=	[MDLuestc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.totconerr];		
					MDLentestc  	=	[MDLentestc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.totconerr];	
					MDLcountestc	=	[MDLcountestc, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.totconerr];

					BICfallout			=	[BICfallout, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.fallout];		 
					AICfallout			=	[AICfallout, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.fallout];		
					MDLlfallout			=	[MDLlfallout, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.fallout];		
					MDLufallout			=	[MDLufallout, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.fallout];		
					MDLentfallout		=	[MDLentfallout, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.fallout];	
					MDLcountfallout		=	[MDLcountfallout, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.fallout];

					BICrecall			=	[BICrecall, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.recall];		 
					AICrecall			=	[AICrecall, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.recall];		
					MDLlrecall			=	[MDLlrecall, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.recall];		
					MDLurecall			=	[MDLurecall, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.recall];		
					MDLentrecall		=	[MDLentrecall, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.recall];	
					MDLcountrecall		=	[MDLcountrecall, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.recall];

					BICprecision		=	[BICprecision, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.BIC.precision];		 
					AICprecision		=	[AICprecision, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.AIC.precision];		
					MDLlprecision		=	[MDLlprecision, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLl.precision];		
					MDLuprecision		=	[MDLuprecision, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLu.precision];		
					MDLentprecision		=	[MDLentprecision, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLent.precision];	
					MDLcountprecision	=	[MDLcountprecision, Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(sprs).Spars.MDLcount.precision];


				end	
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICespc		 		= 	BICespc			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICespc		 		= 	AICespc			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlespc	 		= 	MDLlespc		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuespc 			= 	MDLuespc		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentespc   		= 	MDLentespc  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountespc 		= 	MDLcountespc	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICesfp		 		= 	BICesfp			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICesfp		 		= 	AICesfp			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlesfp	 		= 	MDLlesfp		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuesfp 			= 	MDLuesfp		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentesfp   		= 	MDLentesfp  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountesfp 		= 	MDLcountesfp	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICesfn		 		= 	BICesfn			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICesfn		 		= 	AICesfn			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlesfn	 		= 	MDLlesfn		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuesfn 			= 	MDLuesfn		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentesfn   		= 	MDLentesfn  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountesfn 		= 	MDLcountesfn	;

					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICestp		 		= 	BICestp			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICestp		 		= 	AICestp			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlestp	 		= 	MDLlestp		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuestp 			= 	MDLuestp		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentestp   		= 	MDLentestp  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountestp 		= 	MDLcountestp	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICestn		 		= 	BICestn			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICestn		 		= 	AICestn			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlestn	 		= 	MDLlestn		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuestn 			= 	MDLuestn		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentestn   		= 	MDLentestn  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountestn 		= 	MDLcountestn	;

					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICestc		 		= 	BICestc			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICestc		 		= 	AICestc			;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlestc	 		= 	MDLlestc		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuestc 			= 	MDLuestc		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentestc   		= 	MDLentestc  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountestc 		= 	MDLcountestc	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICfallout	 		= 	BICfallout		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICfallout	 		= 	AICfallout		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlfallout	 		= 	MDLlfallout		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLufallout			= 	MDLufallout		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentfallout		= 	MDLentfallout  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountfallout		= 	MDLcountfallout	;	
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICrecall	 		= 	BICrecall		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICrecall	 		= 	AICrecall		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlrecall	 		= 	MDLlrecall		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLurecall			= 	MDLurecall		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentrecall	 	= 	MDLentrecall  	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountrecall	 	= 	MDLcountrecall	;	
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICprecision	 	= 	BICprecision	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICprecision	 	= 	AICprecision	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlprecision	 	= 	MDLlprecision	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuprecision		= 	MDLuprecision	;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentprecision	 	= 	MDLentprecision ;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountprecision	= 	MDLcountprecision;	
			end
		end
	end	
end
