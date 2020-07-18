%modelgraphs3.m
% recoing Nicola's fig 3, comparing connection reconstruction


function modelgraphs3(Jstor,sparsity,betavec,Nvec,Tvec,topologies)


figwidth = length(betavec);
figheight = length(Tvec);
figaxx = sparsity(length(sparsity)) + 0.2;

% Our x axis for the subplots, here is taking the sparsity from the input which probably won't always be the case
% xvec = sparsity;

for Nt = 1:length(Nvec)



	for topology = 1:length(topologies) % for each complete figure
		h = figure;

		hold on

		j = 1;

	
		for bt = 1:length(betavec)
				
			for Tt = 1:length(Tvec)
														
				% BICespc		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICespc		;
				% AICespc		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICespc		;
				% MDLlespc	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlespc	 	;
				% MDLuespc	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuespc 	;
				% MDLentespc   = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentespc   ;
				% MDLcountespc = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountespc ;% 

				% figaxy = max([BICespc,AICespc,MDLlespc,MDLuespc,MDLentespc,MDLcountespc]) + 0.2;% 

				% subplot(figwidth,figheight,j)
				% axis([0 figaxx 0 figaxy])
                % hold on
				% plot([sparsity], [BICespc		],'b-o')
				% plot([sparsity], [AICespc		],'r-o')
				% plot([sparsity], [MDLlespc		],'m-o')
				% plot([sparsity], [MDLuespc		],'g-o')
				% plot([sparsity], [MDLentespc  	],'y-o')
				% plot([sparsity], [MDLcountespc	],'k-o')

				%ylabel(['T = ', Tvec(Tt)])
			
			% More generalized

				xBIC		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICrecall		;
				xAIC		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICrecall		;
				xMDLl		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlrecall		;
				xMDLu		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLurecall		;
				xMDLent 	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentrecall  	;
				xMDLcount	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountrecall	; 

				yBIC		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICprecision		;		
				yAIC		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICprecision		;		
				yMDLl		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlprecision	; 
				yMDLu		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuprecision	; 
				yMDLent 	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentprecision 	;
				yMDLcount	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountprecision; 






				 figaxy = max([yBIC,yAIC,yMDLl,yMDLu,yMDLent,yMDLcount]) + 0.2;% 

				 subplot(figwidth,figheight,j)
				 axis([0 figaxx 0 figaxy])
                 hold on
				 plot([xBIC	], [yBIC	],'b-o')
				 plot([xAIC	], [yAIC	],'r-o')
				 plot([xMDLl	], [yMDLl	],'m-o')
				 plot([xMDLu	], [yMDLu	],'g-o')
				 plot([xMDLent ], [yMDLent ],'y-o')
				% plot([xMDLcount],[yMDLcount],'k-o')

				%ylabel(['T = ', Tvec(Tt)])

				
				j = j+1;
	
			end
	
		end
		
		ax = findobj(h,'Type','Axes');
        
      	%for 1:length(ax) %LABELs SURE ARE A BITCH HUH

       		% title(ax(4), ['N=',num2str(Nvec(1))])
	    	% ylabel(ax(9),['T=1E',num2str(log10(Tvec(1)))])
	    	% title(ax(8), ['N=',num2str(Nvec(2))])
	    	% ylabel(ax(4),['T=1E',num2str(log10(Tvec(1)))])
	    	% title(ax(7), ['N=',num2str(Nvec(3))])
	    	% ylabel(ax(1),['T=1E',num2str(log10(Tvec(2)))])
	
       		% xlabel(ax(1),['\gamma_J = ',num2str(Jstor(9).(method(jj)))])
       		% xlabel(ax(2),['\gamma_J = ',num2str(Jstor(8).(method(jj)))])
       		% xlabel(ax(3),['\gamma_J = ',num2str(Jstor(7).(method(jj)))])
       		% xlabel(ax(4),['\gamma_J = ',num2str(Jstor(6).(method(jj)))])
       		% xlabel(ax(5),['\gamma_J = ',num2str(Jstor(5).(method(jj)))])
       		% xlabel(ax(6),['\gamma_J = ',num2str(Jstor(4).(method(jj)))])
       		% xlabel(ax(7),['\gamma_J = ',num2str(Jstor(3).(method(jj)))])
       		% xlabel(ax(8),['\gamma_J = ',num2str(Jstor(2).(method(jj)))])
       		% xlabel(ax(9),['\gamma_J = ',num2str(Jstor(1).(method(jj)))])
    	%end

       [ax1,h1]=suplabel('Reconstruction Error vs sparsity');
       % [ax2,h2]=suplabel('Monte Carlo','y');
       % sgtitle({'Comparative J values using ',method(jj),' inference method for \beta = ',num2str(beta)})
	end

end
