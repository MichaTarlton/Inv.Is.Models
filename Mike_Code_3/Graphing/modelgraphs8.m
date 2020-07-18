% modelgraphs8.m 
% This time the total correct v. total incorrect

%function figstor = modelgraphs6(OverStruct,sprsvec,betavec,Nvec,Tvec,topologies)
function figstor = modelgraphs8(OverStruct)



% mix and match rder of these to match how you want to index or compare values

% the only ones you neeed to worry about
% the bottom loop will basically be your xxaxis
% and then you want it in the order of whatever you are comparing it against
% so the next outer loop will effectively be the yyaxis

%outside of that the next two outer loops will in one way or another be the Xaxis & Yaxis
% I'm not sure if the order maters for now, but perhaops does for the order of the subplots
% so perhaps keep in order
% > Yaxis
%		> Xaxis
%			> yyaxis
%				> xxaxis for each model
% And then orders above that are not especiall important except for making sure to keep track of figure combos

% Parameter fixins
%
% Topology
%topovec = topologies;
topovec = OverStruct.topologies;
ltopo = length(topovec);


% Sparsity
%sprsvec = sprsvec;
sprsvec = OverStruct.sparsity;
lsprs = length(sprsvec);


% Beta
%betavec = betavec;
betavec = OverStruct.betavec;
lbeta = length(betavec);

% T
Tvec = OverStruct.Tvec;
lT = length(Tvec);

% N
Nvec = OverStruct.Nvec;
lN = length(Nvec);

% models
colorvec = {'b-o','r-o','m-o','g-o','y-o','k-o'};
%colorvec = {'b','r','m','g','y','k'};
modelvec = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  },{'MDLcount'}];
modelvec = modelvec; 
lm = length(modelvec);

% scores (tend to be the y subplot axis)
% so before it was the sparsity vector essentially



figwidth = length(betavec);
figheight = length(Tvec);
%figaxx = sprsvec(lsprs) + 0.2;


figstor = struct('jax',{});

jax = 1;

% making this example for tracking prcision v. recall
for topology = 1:length(topovec) % for each complete figure

	for St = 1:lsprs

		jax = 1;

		for Tt = 1:length(Tvec) % once we add trials we'll have to add another loop here
			T = Tvec(Tt);

			for bt = 1:length(betavec)
				beta = betavec(bt);
			
				
				
				figstor(jax).jax = jax;

				for mt = 1:lm
					
					model = modelvec{mt};

					% ONLY FOR THIS CASE DELeTE ELSEWHERE!
					figstor(jax).Ax(topology).topo.Spars(St).([model,'xx']) = Nvec./100;

					%commented out the x axis for now
					% xxaxisvec = []; % create blank xxaxis values vector to be occuped below
					yyaxisvec = []; % create blank yyaxis values vector to be occuped below

					for Nt = 1:length(Nvec)				
			
						%fallout			=  	OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).fallout	;
						%recall			=	OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).recall	;
						%precision		=	OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).precision;
							
						tottrue = 	OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).tpconerr + ...
										OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).tnconerr	;
						totfalse = 	OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).fpconerr + ...
										OverStruct(1).Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.(model).fnconerr	;

						ratiotf =	tottrue/totfalse;

						%since I'm doing a comparison of the T/F ratio to nodes, I don't even need the x-axis here
						% but I do need to set it down lower
						% I think I might div 10 the lower scale for correct graph ratios
						%xxaxisvec = [];
						yyaxisvec = [yyaxisvec,ratiotf];

						% in this first case doing the AUPRC per value of N, precision v recall
			
						% xxaxisvec = [xxaxisvec,recall	];
						% altyyaxisvec = [yyaxisvec,precision];
						% yyaxisvec = [yyaxisvec,fallout];

						%figstor(jax).Ax(topology).topo.Spars(St).([model,'xx']) = xxaxisvec;
						figstor(jax).Ax(topology).topo.Spars(St).([model,'yy']) = yyaxisvec;
					end
					
					

												 			
					%figstor(jax).jax = jax;
				end
						
				jax = jax + 1;
				
			end
		end
	end
end


% now for creating the graph (Maybe even set in above just below model loop)
for topology = 1:length(topovec) % for each complete figure

	for St = 1:lsprs
        h = figure;
		j = 1;

		for Tt = 1:length(Tvec) % once we add trials we'll have to add another loop here
			T = Tvec(Tt);
        
            
			for bt = 1:length(betavec)
				beta = betavec(bt);
			
			
				% figaxy = max([JBIC,JAIC,JMDLl,JMDLu,JMDLent,JMDLcount]) + 0.2; for creating subplot axis length sttandardized

         		   subplot(figheight,figwidth,j)

	 			 
	 			
	             %hold on

	            for mt = 1:lm
					model = modelvec{mt};
					modelcolor = colorvec{mt};
					 figmaxaxx = max([figstor(j).Ax(topology).topo.Spars(St).([model,'xx'])]) + 0.1;
                     figmaxaxy = max([figstor(j).Ax(topology).topo.Spars(St).([model,'yy'])]) + 0.1;
                    % figminaxx = min([figstor(j).Ax(topology).topo.Spars(St).([model,'xx'])]) + 0.1;
					% figminaxy = min([figstor(j).Ax(topology).topo.Spars(St).([model,'yy'])]) + 0.1;
                    % axis([figminaxx figmaxaxx figminaxy figmaxaxy]) 
                    axis([0 figmaxaxx 0 figmaxaxy]) 
                    % axis([0 1 0 1]) 
					plot(figstor(j).Ax(topology).topo(topology).Spars(St).([model,'xx']), figstor(j).Ax(topology).topo(topology).Spars(St).([model,'yy']),modelcolor)
					
					% yyaxis right
					% ylabel('Right Side')

					%scatter(figstor(j).Ax(topology).topo.Spars(St).([model,'xx']), figstor(j).Ax(topology).topo.Spars(St).([model,'yy']),30,modelcolor,'filled')
					hold on
				end
				
				j = j +1;

			end
		end

 			axh = findobj(h,'Type','Axes');
            axh2 = axh;
            axh2(4).YAxisLocation = 'right';
            axh2(1).YAxisLocation = 'right';

%  			ax2 = axes('XAxisLocation','top',...
%          'YAxisLocation','right',...
%          'Color','none',...
%          'XColor','k','YColor','k');

         	title(axh(6), ['\beta = ',num2str(betavec(1))],'FontSize',20,'FontWeight','bold')
         	title(axh(5), ['\beta = ',num2str(betavec(2))],'FontSize',20,'FontWeight','bold')
         	title(axh(4), ['\beta = ',num2str(betavec(3))],'FontSize',20,'FontWeight','bold')
         	%ylabel(axh(3),['T=1E',num2str(log10(Tvec(2)))],'FontWeight','bold')
         	ylabel(axh2(1),['T=1E',num2str(log10(Tvec(2)))],'FontSize',20,'FontWeight','bold')
         	%ylabel(axh(6),['T=1E',num2str(log10(Tvec(1)))],'FontWeight','bold')
         	ylabel(axh2(4),['T=1E',num2str(log10(Tvec(1)))],'FontSize',20,'FontWeight','bold')

         	%[ax1,h1]=suplabel(['AUPRC in order of values of N : ', num2str(Nvec)]);
    	    %[ax2,h2]=suplabel('Precision: $\frac{TP}{TP + FP}$','y'); figure out the latex here in a sec
    	    [ax1,h1]=suplabel(['Number of features in order:' , num2str(Nvec)]);
    	    [ax2,h2]=suplabel(['True to False reconstruction ratio'],'y');

    	    [ax,h3]=suplabel(['True - False Ratio for Topology: ',topovec{topology}, ', and Sparsity: ',num2str(sprsvec(St))] ,'t');
			set(h3,'FontSize',30)
    	    %set(h1,'FontSize',15,'FontWeight','bold','interpreter','latex')
    	    %set(h2,'FontSize',15,'FontWeight','bold','interpreter','latex')
    	    set(h1,'interpreter','latex','FontSize',25,'FontWeight','bold')
 			set(h2,'interpreter','latex','FontSize',25,'FontWeight','bold')
        	%sgtitle(['AUROC in order of values of N : ', num2str(Nvec)])

        	legend(modelvec)
    

%         	title(axh(8), ['N=',num2str(Nvec(2))])
%         	ylabel(axh(6),['T=1E',num2str(log10(Tvec(2)))])
%        	title(axh(7), ['N=',num2str(Nval(3))])
%        	ylabel(axh(3),['T=1E',num2str(log10(Tval(3)))])
% 
%         xlabel(axh(1),['N: Precision v. Recall'])
%         xlabel(axh(2),['N: Precision v. Recall'])
%         xlabel(axh(3),['N: Precision v. Recall'])
%         xlabel(axh(4),['N: Precision v. Recall'])
%         xlabel(axh(5),['N: Precision v. Recall'])
%         xlabel(axh(6),['N: Precision v. Recall'])
%         %xlabel(axh(7),['\gamma_J = ',num2str(Jstor(3).(method(hh)))])
%         %xlabel(axh(8),['\gamma_J = ',num2str(Jstor(2).(method(hh)))])
%         %xlabel(axh(9),['\gamma_J = ',num2str(Jstor(1).(method(hh)))])

	end
end


						
