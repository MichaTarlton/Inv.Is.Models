%% organize multi all struct and grpah
% Multiallstruct should be a massive structure put together from all relevant run parameters
% comp will set what tpe of figures we want
% 		1. Perconerr v beta
%		2. AUROC v beta
%		3. AUPRC v beta
%		4. rrerr v beta
%		5. Perconerr v beta for a symmetric array
%		6. PPC v beta
% Ideally we can change up the x-axis
% I want to be able to do something more like Fig 1 from bento
%  Changing to beusable with latest version of foris

function [figstor,newstruct,dupes] = multiallgraph3(overstruct,comp,topdir,time)

multiallstruct = overstruct.list; % Or I could just pass it in as this

lmas = length(multiallstruct);

 %namevec      	=	 [multiallstruct.name];     
 Tvec         	=	 [multiallstruct.T];        
 Nvec         	=	 [multiallstruct.N];        
 betavec      	=	 [multiallstruct.beta];     
 %sprsvecvec   	=	 [multiallstruct.sprsvec];  
 %topologyvec  	=	 [multiallstruct.topology];

 %Jcontruvec   	=	 [multiallstruct.Jcontru];  
 %Jtruvec      	=	 [multiallstruct.Jtru];     
 %htruvec      	=	 [multiallstruct.htru];     
 couplingsvec 	=	 [multiallstruct.couplings];
 cvec         	=	 [multiallstruct.c];        
 %Svec         	=	 [multiallstruct.S];        
 %BLLHvec      	=	 [multiallstruct.BLLH];     


%uniqname	 	  	=	unique(namevec     );
uniqT	 	  	=	unique(Tvec        );
uniqN	 	  	=	unique(Nvec        );
uniqbeta	 	  	=	unique(betavec     );
%uniqsprsvec	  	=	unique(sprsvecvec  );	
%uniqtopology	  	=	unique(topologyvec );	
%uniqJcontru	  	=	unique(Jcontruvec  );	
%uniqJtru	 	  	=	unique(Jtruvec     );
%uniqhtru	 	  	=	unique(htruvec     );
uniqcouplings  	=	unique(couplingsvec);		
uniqc	 	  	=	unique(cvec        );
%uniqS	 	  	=	unique(Svec        );
%uniqBLLH	 	  	=	unique(BLLHvec     );


lut = length(uniqT	 );
lun = length(uniqN	 );
lub = length(uniqbeta);
luc = length(uniqc);

%models = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  },{'MDLcount'}];
models = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  }]; % don't need MDLcount. Might add cross fold instead later
lm = length(models);

topologies = [multiallstruct(1).topology{1,:}] ;% Need to fix above!
%topologies = [2,5]; %all for now, fixing above and outside

ltop = length(topologies);
toponamevec = [{'Cayley Tree'},{'Fully Connected'},{'Independent Pair'},{'2D Ising'},{'Erdős–Rényi'},{'Erdős–Rényi (Pachenko)'},{'Star Netowrk'}];
displaytoponamevec = toponamevec(topologies);



figwidth = lm;
figheight = luc;

figstor = struct('jax',{});

newstruct = struct;

dupes =[];

%cstruct = struct;

loop = 1;

jax = 1; %per individual complete figure

%for masn = 1:lmas % length of struct
for nn = 1:lun % number of nodes
	for tn = 1:lut % number of temp values
		

		
		for cn = 1:luc % number of c val, row i for subplots

			c = uniqc(cn);
			%cstruct(cn).c = c;

			for mn = 1: lm % Models, col j for subplots
			
				model = models{mn};
				figstor(jax).jax = jax; 	% for each subplot
			
				for topn = 1:ltop
				
					%topo = topologies(topn);
					topo = topn;
				
					%xxaxisvec = []; % create blank xxaxis values vector to be occuped below
					%yyaxisvec = []; % create blank yyaxis values vector to be occuped below

					%perconerr		= zeros(1,lub);
					%precision		= zeros(1,lub);
					%recall		= zeros(1,lub);
					%fallout		= zeros(1,lub);
					%Jrrerr		= zeros(1,lub);
					%Jsymrrerr		= zeros(1,lub);
					%Javgsymrrerr	= zeros(1,lub);			

					perconerr		= [];
					precision		= [];
					recall			= [];
					fallout			= [];
					Jrrerr			= [];
					Jsymrrerr		= [];
					Javgsymrrerr	= [];
					symperconerr	= [];
					TFS				= [];
					TFR			 	= [];	

					ERR_perconerr 		= [];
					ERR_precision 		= [];
					ERR_recall 			= [];
					ERR_fallout			= [];
					ERR_Jrrerr			= [];
					ERR_Jsymrrerr		= [];
					ERR_Javgsymrrerr	= [];
					ERR_symperconerr 	= [];
					ERR_TFS				= [];
					ERR_TFR			 	= [];




					for bn = 1:lub % number of beta values & we don't need another loop because the values of the subplots will follow the beta values
					
					
						idx = cvec == uniqc(cn) & Tvec == uniqT(tn) & Nvec == uniqN(nn) & betavec == uniqbeta(bn);
						
						%multiallstruct(idx)
					
						if length(multiallstruct(idx)) > 0 %& length(multiallstruct(idx)) < 2 %add this if you can't find the duplicates
							newstruct(loop).new = multiallstruct(idx);
							loop = loop + 1;
					
							%add data function. Maybe make a funtion here to add data that seems to be missing from BLLH, such as the relative reconstruction error rate as don in model graphs 4
						
							%% vector(bn) = multiallstruct(idx).BLLH(1).trial.(errvalue we want).model  %trying a more refined way of looping element into vector test and immediately switch to old way if failure	
							
							%perconerr(:,bn) = multiallstruct(idx).BLLH(1).trial.perconerr.(model);
							%precision(:,bn) = multiallstruct(idx).BLLH(1).trial.precision.(model);
							%recall(:,bn) 	= multiallstruct(idx).BLLH(1).trial.recall.(model);
							%fallout(:,bn)	= multiallstruct(idx).BLLH(1).trial.fallout.(model);
							%Jrrerr(:,bn)	= multiallstruct(idx).BLLH(1).trial.Jrrerr.(model);
							%Jsymrrerr(:,bn)		= multiallstruct(idx).BLLH(1).trial.Jsymrrerr.(model);
							%Javgsymrrerr(:,bn)	= multiallstruct(idx).BLLH(1).trial.Javgsymrrerr.(model);
							
							%perconerr 		= [perconerr, multiallstruct(idx).BLLH(1).trial.perconerr.(model)];
							%precision 		= [precision, multiallstruct(idx).BLLH(1).trial.precision.(model)];
							%recall 			= [recall, multiallstruct(idx).BLLH(1).trial.recall.(model)];
							%fallout			= [fallout, multiallstruct(idx).BLLH(1).trial.fallout.(model)];
							%Jrrerr			= [Jrrerr, multiallstruct(idx).BLLH(1).trial.Jrrerr.(model)];
							%Jsymrrerr		= [Jsymrrerr, multiallstruct(idx).BLLH(1).trial.Jsymrrerr.(model)];
							%Javgsymrrerr	= [Javgsymrrerr, multiallstruct(idx).BLLH(1).trial.Javgsymrrerr.(model)];
							%%symperconerr 	= [symperconerr, multiallstruct(idx).BLLH(1).trial.symperconerr.(model)];
							%TFS				= [TFS, multiallstruct(idx).BLLH(1).trial.TFS.(model)];
							%TFR			 	= [TFR, multiallstruct(idx).BLLH(1).trial.TFR.(model)];

							%% Now for average values taken over multiple trials
							perconerr 		= [perconerr, 	multiallstruct(idx).BLLH(topo).stats.mean.(model).perconerr];
							precision 		= [precision, 	multiallstruct(idx).BLLH(topo).stats.mean.(model).precision];
							recall 			= [recall, 		multiallstruct(idx).BLLH(topo).stats.mean.(model).recall];
							fallout			= [fallout, 	multiallstruct(idx).BLLH(topo).stats.mean.(model).fallout];
							Jrrerr			= [Jrrerr,		multiallstruct(idx).BLLH(topo).stats.mean.(model).Jrrerr];
							Jsymrrerr		= [Jsymrrerr, 	multiallstruct(idx).BLLH(topo).stats.mean.(model).Jsymrrerr];
							Javgsymrrerr	= [Javgsymrrerr,multiallstruct(idx).BLLH(topo).stats.mean.(model).Javgsymrrerr];
							symperconerr 	= [symperconerr,multiallstruct(idx).BLLH(topo).stats.mean.(model).symperconerr];
							TFS				= [TFS, 		multiallstruct(idx).BLLH(topo).stats.mean.(model).TFS];
							TFR			 	= [TFR, 		multiallstruct(idx).BLLH(topo).stats.mean.(model).TFR];

							%% Now for our measure of error
							ERR_perconerr 		= [ERR_perconerr, 	multiallstruct(idx).BLLH(topo).stats.stderr.(model).perconerr];
							ERR_precision 		= [ERR_precision, 	multiallstruct(idx).BLLH(topo).stats.stderr.(model).precision];
							ERR_recall 			= [ERR_recall, 		multiallstruct(idx).BLLH(topo).stats.stderr.(model).recall];
							ERR_fallout			= [ERR_fallout, 		multiallstruct(idx).BLLH(topo).stats.stderr.(model).fallout];
							ERR_Jrrerr			= [ERR_Jrrerr,		multiallstruct(idx).BLLH(topo).stats.stderr.(model).Jrrerr];
							ERR_Jsymrrerr		= [ERR_Jsymrrerr, 	multiallstruct(idx).BLLH(topo).stats.stderr.(model).Jsymrrerr];
							ERR_Javgsymrrerr	= [ERR_Javgsymrrerr,	multiallstruct(idx).BLLH(topo).stats.stderr.(model).Javgsymrrerr];
							ERR_symperconerr 	= [ERR_symperconerr,	multiallstruct(idx).BLLH(topo).stats.stderr.(model).symperconerr];
							ERR_TFS				= [ERR_TFS, 			multiallstruct(idx).BLLH(topo).stats.stderr.(model).TFS];
							ERR_TFR			 	= [ERR_TFR, 			multiallstruct(idx).BLLH(topo).stats.stderr.(model).TFR];



							%xxaxisvec(bn) = vector
							%yyaxisvec(bn) = vector
						end

						%for finding duplicates, in theory or at least for now I don't have multiple trial run yet
						if length(multiallstruct(idx)) > 1
							dupes = [dupes,find(idx,1)];
						end

						%%checking the parameters that need to be run
						%if c = 1
						%	%cstruct(cn).c = c; % put further up
						%	CN = [N,  ]
						%	TN = [N,  ]
						%	BN = [N,  ]							]
						%	elseif c = 2
					end
				
					%% change this to change the err values you wnat to compare
					
					if comp == 1		% perconerr
						figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =  perconerr;
						figstor(jax).fig(cn).c.(model).topo(topo).xxaxisvec = uniqbeta;
						figstor(jax).fig(cn).c.(model).topo(topo).errorvec = ERR_perconerr;

					elseif comp == 2	% AUROC
						figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =	fallout./recall;
						%figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =	trapz(fallout,recall);	%consider this as well
						figstor(jax).fig(cn).c.(model).topo(topo).xxaxisvec =	uniqbeta;
						figstor(jax).fig(cn).c.(model).topo(topo).errorvec = (ERR_recall + ERR_fallout)./2; % I should have made this earlier for correct err measure, but I'll just take average instead	

					elseif comp == 3 	% AUPRC
						figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =	precision./recall;
						%figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =	trapz(fallout,precision); % consider this as well since this is the correct value to use for "Area Under"	
						figstor(jax).fig(cn).c.(model).topo(topo).xxaxisvec =	uniqbeta;
						figstor(jax).fig(cn).c.(model).topo(topo).errorvec = (ERR_precision + ERR_fallout)./2; % see above

					elseif comp == 4	% Relative recon err, for symmetry reconstruction
						figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =	Jsymrrerr;	
						figstor(jax).fig(cn).c.(model).topo(topo).xxaxisvec =	uniqbeta;
						figstor(jax).fig(cn).c.(model).topo(topo).errorvec = ERR_Jsymrrerr;

					elseif comp == 5	%  Perconerr for a symmetric array
						figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =  symperconerr;
						figstor(jax).fig(cn).c.(model).topo(topo).xxaxisvec = uniqbeta;
						figstor(jax).fig(cn).c.(model).topo(topo).errorvec = ERR_symperconerr;

					elseif comp == 6	% PPC score: TFS/TFR
						figstor(jax).fig(cn).c.(model).topo(topo).yyaxisvec =  TFS./TFR;
						figstor(jax).fig(cn).c.(model).topo(topo).xxaxisvec = uniqbeta;
						figstor(jax).fig(cn).c.(model).topo(topo).errorvec = ERR_TFS./ERR_TFR; % I don't think this is like the one above since they are actually supposed to be done as a ratio

					end

					% just in case
					figstor(jax).fig(cn).c.(model).topo(topo).perconerr 	= perconerr;
					figstor(jax).fig(cn).c.(model).topo(topo).symperconerr 	= symperconerr;
					figstor(jax).fig(cn).c.(model).topo(topo).precision 	= precision;
					figstor(jax).fig(cn).c.(model).topo(topo).recall 		= recall;
					figstor(jax).fig(cn).c.(model).topo(topo).fallout 		= fallout;
					figstor(jax).fig(cn).c.(model).topo(topo).Jrrerr	 	= Jrrerr;
					figstor(jax).fig(cn).c.(model).topo(topo).Jsymrrerr	 	= Jsymrrerr;
					figstor(jax).fig(cn).c.(model).topo(topo).Javgsymrrerr 	= Javgsymrrerr;
					figstor(jax).fig(cn).c.(model).topo(topo).symperconerr	= symperconerr;
					figstor(jax).fig(cn).c.(model).topo(topo).TFS			= TFS;
					figstor(jax).fig(cn).c.(model).topo(topo).TFR			= TFR;

					figstor(jax).fig(cn).c.(model).topo(topo).AUROC = trapz(fallout,recall);
					figstor(jax).fig(cn).c.(model).topo(topo).AUPRC = trapz(fallout,precision);

					figstor(jax).fig(cn).c.(model).topo(topo).ERR_perconerr 	= ERR_perconerr;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_symperconerr 	= ERR_symperconerr;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_precision 	= ERR_precision;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_recall 		= ERR_recall;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_fallout 		= ERR_fallout;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_Jrrerr	 	= ERR_Jrrerr;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_Jsymrrerr	 	= ERR_Jsymrrerr;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_Javgsymrrerr 	= ERR_Javgsymrrerr;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_symperconerr	= ERR_symperconerr;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_TFS			= ERR_TFS;
					figstor(jax).fig(cn).c.(model).topo(topo).ERR_TFR			= ERR_TFR;

					figstor(jax).fig(cn).cval = c;
					figstor(jax).T = uniqT(tn);
					figstor(jax).N = uniqN(nn);
				
					%% store the coordinat value vector, not sure how I want to implement the figure creation just yet
					%% I could also just throw the figure generator below, so, idk
					%figstor(jax).Ax(topology).topo.Spars(St).([model,'xx']) = xxaxisvec;
					%figstor(jax).Ax(topology).topo.Spars(St).([model,'yy']) = yyaxisvec;
				end
			end
		end
		jax = jax+1;
	end
end

dupes = unique(dupes);

%%%%% On to the figure gen

% pulled from a an online palette
c1 = [249, 65, 68]./255;
c2 = [243, 114, 44]./255;
c3 = [248, 150, 30]./255;
c4 = [249, 199, 79]./255;
c5 = [144, 190, 109]./255;
c6 = [67, 170, 139]./255;
c7 = [87, 117, 144]./255;

% anice flat black/grey
%cbg = [53, 53, 53]./255;
cbg = [60, 60, 60]./255;

% a light grey
clg = [217, 217, 217]./255;


%for masn = 1:lmas % length of struct
%colorvec = {'b-o','r-o','m-o','g-o','c-o','k-o'};
%colorvec = {c1,c2,c3,c4,c5,c6,c7};
colorvec = {c1,c3,c4,c5,c7};


for jnx = 1:length(figstor) %total number of figure as decided by NxT
	
	%h = figure;
	h = figure('Position', get(0, 'Screensize'));
	j = 1;

	for cn = 1:luc % number of c val, row i for subplots
		c = uniqc(cn);
		
		for mn = 1: lm % Models, col j for subplots
		
			model = models{mn};
			
			subplot(figheight,figwidth,j)
			
			avgyyvec = []; %placeholder to store all y-axis values in order to be averaged over all values for all topos
			avgxxvec = [];
			avgaurocvec = [];
			avgauprcvec = []; %placeholder to store all area under the curve values in order to be averaged over all values for all topos	

			avgaurocvec = [figstor(jnx).fig(cn).c.(model).topo.AUROC];
			avgauprcvec = [figstor(jnx).fig(cn).c.(model).topo.AUPRC];

			avgyyvec  = [figstor(jnx).fig(cn).c.(model).topo.yyaxisvec];
			avgxxvec  = [figstor(jnx).fig(cn).c.(model).topo.xxaxisvec]; 


			ymax = max(avgyyvec) + 0.03;
			%ymin = min(yyaxisvec) - 0.03;
			ymin = 0;
			xmax = max(avgxxvec) + 0.03;
			xmin = min(avgxxvec) - 0.03;
			%xmin = 0;


			%axis([0 0.5 0 0.3]) 
			%axis([xmin xmax ymin ymax]) 
			

			for topn = 1:ltop
				topocolor = colorvec{topn};
				%topo = topologies(topn);
				topo = topn;
			
				%% change this to change the err values you wnat to compare
				yyaxisvec  = figstor(jnx).fig(cn).c.(model).topo(topo).yyaxisvec;
				xxaxisvec  = figstor(jnx).fig(cn).c.(model).topo(topo).xxaxisvec; 
				errvec  = figstor(jnx).fig(cn).c.(model).topo(topo).errorvec;

				axis([xmin xmax 0 ymax]) 
				%plot(xxaxisvec,yyaxisvec,topocolor,'LineWidth',2)
				%plot(xxaxisvec,yyaxisvec,'-o' ,'Color', topocolor)
				errorbar(xxaxisvec,yyaxisvec,errvec,'Marker','o','Color', topocolor)

				set(gca,'Color',cbg)
				%xlabel(num2str(j))
				hold on	

				% yyaxis right
				% ylabel('Right Side')

				%AUROC = figstor(jnx).fig(cn).c.(model).topo(topo).AUROC;
				%AUPRC = figstor(jnx).fig(cn).c.(model).topo(topo).AUPRC;

				%avgyyvec = [avgyyvec,yyaxisvec];
				
				%avgaurocvec =[avgaurocvec,AUROC];
				%avgauprcvec =[avgauprcvec,AUPRC];

			end


			%% was going to use for xlabels, but realized I can just use the yyaxisvec
			%Jsymrrerr = figstor(jax).fig(cn).c.(model).topo(topo).Jsymrrerr;
			%perconerr = figstor(jax).fig(cn).c.(model).topo(topo).perconerr;

			if comp == 1	%perconerr
				xlabel(['err rate = ',num2str(mean(avgyyvec))])
			elseif comp == 2
				%xlabel(['AUROC = ',num2str(mean(avgaurocvec))]) % The AUC method doesn't really work for this application atm, see notes for details
				xlabel(['avg ROC score = ',num2str(mean(avgyyvec))])
			elseif comp == 3 
				%xlabel(['AUPRC = ',num2str(mean(avgauprcvec))]) % The AUC method doesn't really work for this application atm, see notes for details
				xlabel(['avg PRC score = ',num2str(mean(avgyyvec))])
			elseif comp == 4
				xlabel(['\gamma_J = ',num2str(mean(avgyyvec))]) %% Maybe do an area under the curve func for this?
			elseif comp == 5
				xlabel(['sym err rate = ',num2str(mean(avgyyvec))])
			elseif comp == 6
				xlabel(['avg PPC score= ',num2str(mean(avgyyvec))])
			end
			
			j = j +1;
		end
	end
	j = j-1;

	% in theory j = figheight*figwidth
	%for line =




	axh = findobj(h,'Type','Axes');
	axh2 = axh;	
	% axh2(4).YAxisLocation = 'right'; %not sure why these aren't working suddenly
 	% axh2(1).YAxisLocation = 'right';
 
 	% ax2 = axes('XAxisLocation','top',...
	% 'YAxisLocation','right',...
	% 'Color','none',...
	% 'XColor','k','YColor','k');
	jj = 0;
	for labels = 1:figwidth
		
		title(axh(j - jj), [models{labels}],'FontSize',15,'FontWeight','bold')
		
		
		jj= jj+1;
	end
		
	for sidelabels = 1:figheight
		%sl = [15,10,5];
		sl = (figheight*figwidth)-((sidelabels-1)*figwidth);

		%ylabel(axh2(sl(sidelabels)),['C : ',num2str(uniqc(sidelabels))],'FontSize',15,'FontWeight','bold') % not finished here
            ylabel(axh2(sl),['C : ',num2str(uniqc(sidelabels))],'FontSize',15,'FontWeight','bold') % not finished here
	end

	lgnd = legend(displaytoponamevec,'Location','east');
	%set(lgnd,'color','none');
	%set(lgnd,'color',cbg,'position',[.03 .7 .06 .2]);
	set(lgnd,'color',cbg,'position',[0.92 .8 .06 .1]);

	%set(lgnd,'color',cbg);
	c = lgnd.TextColor;
	lgnd.TextColor = 'white';

	N = figstor(jnx).N;
	T = figstor(jnx).T;

	%0.03203125,0.704462326261888,0.063330077187857,0.193050475493782

	if comp == 1		% perconerr
		%[ax2,h2]=suplabel('C: Coordination No.','y');
		[ax,h2]=suplabel([{'Percent Misclassification (Accuracy)'},{'for Coordination number: C'}],'y');
		%[ax,h1]=suplabel('For \beta values: 0.1, 0.2, 0.3, 0.4');
		%[ax,h1]=suplabel('Average Error Rate for all Topologies');
		[ax,h1]=suplabel([{'Average Error Rate for all Topologies'},{'For \beta values'}]);
		%[ax,h3]=suplabel(['N=',num2str(N),' | T=1E',num2str(log10(T))] ,'t');
		%sgtitle([{'Misclassification Error for N='},{num2str(N),' | T=1E',num2str(log10(T))}])
		%sgtitle('Misclassification Error')
		sgtitle({
			['Misclassification Error ']
			['N=',num2str(N),' | T=1E',num2str(log10(T))]
			})
	elseif comp == 2	% AUROC
		[ax,h2]=suplabel([{'Recall v Fallout'},{'for Coordination number: C'}],'y'); %maybe put the ROC in the main title
		%[ax1,h1]=suplabel('For \beta values: 0.1, 0.2, 0.3, 0.4');
		[ax,h1]=suplabel([{'Average "ROC score" for all Topologies'},{'For \beta values'}]);
		%[ax,h3]=suplabel(['N=',num2str(N),' | T=1E',num2str(log10(T))] ,'t');
		sgtitle({
			['Reciever Operator Characteristic (ROC)']
			['N=',num2str(N),' | T=1E',num2str(log10(T))]
			})
	elseif comp == 3 	% AUPRC
		[ax,h2]=suplabel([{'Precision v Fallout'},{'for Coordination number: C'}],'y');			
		%[ax,h1]=suplabel('For \beta values: 0.1, 0.2, 0.3, 0.4');
		[ax,h1]=suplabel([{'Average "PRC score" for all Topologies'},{'For \beta values'}]);
		%[ax,h3]=suplabel(['N=',num2str(N),' | T=1E',num2str(log10(T))] ,'t');
		sgtitle({
			['Precision Recall Curve (PRC)']
			['N=',num2str(N),' | T=1E',num2str(log10(T))]
			})
	elseif comp == 4	% Relative recon err, for symmetry reconstruction
		[ax,h2]=suplabel([{'Relative Reconstruction Error: '},{'for Coordination number: C'}],'y');
		%[ax,h1]=suplabel('For \beta values: 0.1, 0.2, 0.3, 0.4');
		[ax,h1]=suplabel([{'Average RRerr for all Topologies'},{'For \beta values'}]);
		%[ax,h3]=suplabel(['N=',num2str(N),' | T=1E',num2str(log10(T))] ,'t');
		sgtitle({['Relative Reconstruction Error for Symmetry Optimization']
			['N=',num2str(N),' | T=1E',num2str(log10(T))]
			})
	elseif comp == 5		% perconerr with symmetry
		%[ax2,h2]=suplabel('C: Coordination No.','y');
		[ax,h2]=suplabel([{'Percent Misclassification (Accuracy)'},{'for Coordination number: C'}],'y');
		%[ax,h1]=suplabel('For \beta values: 0.1, 0.2, 0.3, 0.4');
		%[ax,h1]=suplabel('Average Error Rate for all Topologies');
		[ax,h1]=suplabel([{'Average Error Rate for all Topologies'},{'For \beta values'}]);
		%[ax,h3]=suplabel(['N=',num2str(N),' | T=1E',num2str(log10(T))] ,'t');
		%sgtitle([{'Misclassification Error for N='},{num2str(N),' | T=1E',num2str(log10(T))}])
		%sgtitle('Misclassification Error')
		sgtitle({
			['Misclassification Error for symmetry optimization']
			['N=',num2str(N),' | T=1E',num2str(log10(T))]
			})
	elseif comp == 6		% Positive Precision Curve
		%[ax2,h2]=suplabel('C: Coordination No.','y');
		[ax,h2]=suplabel([{'Positive Precision Curve (PPC)'},{'for Coordination number: C'}],'y');
		%[ax,h1]=suplabel('For \beta values: 0.1, 0.2, 0.3, 0.4');
		%[ax,h1]=suplabel('Average Error Rate for all Topologies');
		[ax,h1]=suplabel([{'Average Error Rate for all Topologies'},{'For \beta values'}]);
		%[ax,h3]=suplabel(['N=',num2str(N),' | T=1E',num2str(log10(T))] ,'t');
		%sgtitle([{'Misclassification Error for N='},{num2str(N),' | T=1E',num2str(log10(T))}])
		%sgtitle('Misclassification Error')
		sgtitle({
			['Misclassification Error for symmetry optimization']
			['N=',num2str(N),' | T=1E',num2str(log10(T))]
			})
	end

	set(h1,'FontSize',17);
	set(h2,'FontSize',17);
 	gcf;
 	set(gcf, 'InvertHardcopy', 'off')
saveas(gcf,[topdir,'\',time(1:5),'Multi_Graph ','N=',num2str(N),' T=1E',num2str(log10(T)),'Comparsion Method',num2str(comp),'_',time(6:12),'.png'])
end




%%%% Extra Code for sandboxin 


%% For fixing issue where some BLLH are missing topology 2 ()
% for nnn = 1: length(multiallstruct);
% 	if length(multiallstruct(nnn).BLLH) == 5
% 		multiallstruct(nnn).BLLH = [multiallstruct(nnn).BLLH(1),multiallstruct(nnn).BLLH(1),multiallstruct(nnn).BLLH([2:end])]
% 	end
% end



%% To delete rows from a struct:
% multiallstruct(13:24) = [];
%% fucking with loop headers, you cna remove this
%for masn = 1:lmas % length of struct
%
%	for tn = 1:lut % number of temp values
%	
%		for nn = 1:lun % number of nodes
%	
%			for mm = 1:lm	%number of models
%
%				for cn = 1:luc % number of c val
%
%					for bn = 1:lub % number of beta values
%
%						for topn = 1:ltop % number of topologies
%
%							err rates for these values and store



% multiallstruct(loop).BLLH.fieldnames:
% 	NodeModel
%     theta
%     hrecon
%     Jrecon
%     Jcon
%     Jasym
%     symrate
%     asymrate
%     totconerr
%     perconerr
%     tnconerr
%     tpconerr
%     fnconerr
%     fpconerr
%     fallout
%     recall
%     precision
%     asymrow
%     asymcol
%     Jsymrecon
%     Jsymcon% 