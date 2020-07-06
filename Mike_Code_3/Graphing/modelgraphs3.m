%modelgraphs3.m


function modelgraphs3(Jstor,sparsity,betavec,Nvec,Tvec,topologies)

close all 

figwidth = length(betavec)
figheight = length(Tvec)

for Nt = 1:length(Nvec)



	for topology = 1:length(topologies) % for each complete figure
		h = figure;

		hold on

		j = 1;

	
		for bt = 1:length(betavec)
				
			for Tt = 1:length(Tvec)
														
				BICespc		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICespc		;
				AICespc		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICespc		;
				MDLlespc	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlespc	 	;
				MDLuespc	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuespc 	;
				MDLentespc   = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentespc   ;
				MDLcountespc = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountespc ;

				axis([0 sparsity(length(sparsity)) 0 0.2])
				subplot(figwidth,figheight,j)
				plot([sparsity], [BICespc		],'b-o')
				hold on
				plot([sparsity], [AICespc		],'r-o')
				plot([sparsity], [MDLlespc		],'m-o')
				plot([sparsity], [MDLuespc		],'g-o')
				plot([sparsity], [MDLentespc  	],'y-o')
				plot([sparsity], [MDLcountespc	],'k-o')

				%ylabel(['T = ', Tvec(Tt)])

				
				j = j+1;
	
			end
	
		end
		
		ax = findobj(h,'Type','Axes');
        
        title(ax(9), ['N=',num2str(Nval(1))])
        ylabel(ax(9),['T=1E',num2str(log10(Tval(1)))])
        title(ax(8), ['N=',num2str(Nval(2))])
        ylabel(ax(6),['T=1E',num2str(log10(Tval(2)))])
        title(ax(7), ['N=',num2str(Nval(3))])
        ylabel(ax(3),['T=1E',num2str(log10(Tval(3)))])

        xlabel(ax(1),['\gamma_J = ',num2str(Jstor(9).(method(jj)))])
        xlabel(ax(2),['\gamma_J = ',num2str(Jstor(8).(method(jj)))])
        xlabel(ax(3),['\gamma_J = ',num2str(Jstor(7).(method(jj)))])
        xlabel(ax(4),['\gamma_J = ',num2str(Jstor(6).(method(jj)))])
        xlabel(ax(5),['\gamma_J = ',num2str(Jstor(5).(method(jj)))])
        xlabel(ax(6),['\gamma_J = ',num2str(Jstor(4).(method(jj)))])
        xlabel(ax(7),['\gamma_J = ',num2str(Jstor(3).(method(jj)))])
        xlabel(ax(8),['\gamma_J = ',num2str(Jstor(2).(method(jj)))])
        xlabel(ax(9),['\gamma_J = ',num2str(Jstor(1).(method(jj)))])
    
        [ax1,h1]=suplabel('Generated through inference methods');
        [ax2,h2]=suplabel('Monte Carlo','y');
        sgtitle({'Comparative J values using ',method(jj),' inference method for \beta = ',num2str(beta)})
	end

end
