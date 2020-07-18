% modelgraphs4.m
% making a comparison graph for J reconstruction


function Jstor2 = modelgraphs4(Jstor,sparsity,betavec,Nvec,Tvec,topologies)

% First grab the correct Js

figwidth = length(betavec);
figheight = length(Tvec);
figaxx = sparsity(length(sparsity)) + 0.2;

TotJBICrrerr	= [];
TotJAICrrerr	= [];
TotJMDLlrrerr	= [];
TotJMDLurrerr	= [];
TotJMDLentrrerr  = [];
TotJMDLcountrrerr= [];

for St = 1:length(sparsity)

	for Nt = 1:length(Nvec)
	
	
	
		for topology = 1:length(topologies) % for each complete figure
			
			%h = figure;
	
			%hold on
		
			%j = 1;
			
			for bt = 1:length(betavec)
					
				for Tt = 1:length(Tvec) % once we add trials we'll have to add another loop here
	
					Jtru		= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.Jtru{1}				;								
					JBIC		= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.BIC.Jrecon			;
					JAIC		= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.AIC.Jrecon			;
					JMDLl		= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLl.Jrecon		 	;
					JMDLu		= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLu.Jrecon	 		;
					JMDLent   	= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLent.Jrecon	   		;
					JMDLcount 	= Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLcount.Jrecon	 	;
	
					Jtruvec			=	Jtru(:);
					JBICvec			=	JBIC(:);
					JAICvec			=	JAIC(:);
					JMDLlvec		=	JMDLl(:);
					JMDLuvec		=	JMDLu(:);
					JMDLentvec    	=	JMDLent(:);  
					JMDLcountvec 	=	JMDLcount(:); 
	
					%% Relative reconstruction error Nguyen 17, eq. 114
	       			% Implement theis elsewhere just coding out now
	       			% use to give a numerical value to how reliable inference method is
	       			% probably actually Jgrpahs3
	       			% sum from i where i<j ?? So essentially no self-interaction and no symmetrical interactions
	       			% by only using the upper triangle I should achieve the same thing
	       			%Jtruvec		 = sqrt(sum((Jtruvec		- Jtruvec).^2 / sum(Jtruvec.^2)));
	       			
	       			JBICrrerr		 = sqrt(sum((JBICvec		- Jtruvec).^2 / sum(Jtruvec.^2)));
	       			JAICrrerr		 = sqrt(sum((JAICvec		- Jtruvec).^2 / sum(Jtruvec.^2)));
	       			JMDLlrrerr		 = sqrt(sum((JMDLlvec		- Jtruvec).^2 / sum(Jtruvec.^2)));
	       			JMDLurrerr		 = sqrt(sum((JMDLuvec		- Jtruvec).^2 / sum(Jtruvec.^2)));
					JMDLentrrerr     = sqrt(sum((JMDLentvec  	- Jtruvec).^2 / sum(Jtruvec.^2)));
					JMDLcountrrerr 	 = sqrt(sum((JMDLcountvec	- Jtruvec).^2 / sum(Jtruvec.^2)));

					TotJBICrrerr	  	= 	[TotJBICrrerr		, JBICrrerr		];
					TotJAICrrerr	  	= 	[TotJAICrrerr		, JAICrrerr		];
					TotJMDLlrrerr	  	= 	[TotJMDLlrrerr		, JMDLlrrerr	];
					TotJMDLurrerr	  	= 	[TotJMDLurrerr		, JMDLurrerr	];
					TotJMDLentrrerr   	= 	[TotJMDLentrrerr  	, JMDLentrrerr  ];
					TotJMDLcountrrerr 	= 	[TotJMDLcountrrerr	, JMDLcountrrerr];
					
					
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.BICrerr			= 	JBICrrerr		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.AICrerr			= 	JAICrrerr		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLlrerr		= 	JMDLlrrerr		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLurerr		= 	JMDLurrerr		;
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLentrerr		= 	JMDLentrrerr   	;	
					Jstor(topology).Topo(bt).beta(Nt).N(Tt).T(St).Spars.MDLcountrerr	= 	JMDLcountrrerr 	;	
	
				
	
	
	 
	 				% figaxy = max([JBIC,JAIC,JMDLl,JMDLu,JMDLent,JMDLcount]) + 0.2;
	 
	 				% subplot(figwidth,figheight,j)
					% axis([0 figaxx 0 figaxy])
	              	 % hold on
					% plot([sparsity], [BIC		],'b-o')
					% plot([sparsity], [AIC		],'r-o')
					% plot([sparsity], [MDLl		],'m-o')
					% plot([sparsity], [MDLu		],'g-o')
					% plot([sparsity], [MDLent  	],'y-o')
					% plot([sparsity], [MDLcount	],'k-o')
	 
					%ylabel(['T = ', Tvec(Tt)])
	
% going to do my reconstruction plots first
% lmao this creates like 97 different plots
h = figure;
hold on 
for j = 1:4
subplot(2,2,j)

    % scatter(JBICvec,Jtruvec,1,'b','o')
    % scatter(JAICvec,Jtruvec,1,'r','x')
    % scatter(JMDLlvec,Jtruvec,1,'m','^')
	% scatter(JMDLuvec,Jtruvec,1,'b','o') 
	% scatter(JMDLentvec,Jtruvec,1,'k','x') 
	% scatter(JMDLcountvec,Jtruvec,1,'c','^') 

hold on
%axis([-1 1 -1 1])
axis([-0.2 0.2 -0.2 0.2])
refline(1,0)
refline
end







				%for jj = 1:6
	
					%	h = figure;
					%	hold on 
					%	
					%	for j = 1:4
					%		subplot(2,2,j)
					%		
	        		%		if jj == 1
	        		%		    scatter(JBICvec,Jtruvec,1,'b','o')
	        		%		elseif jj == 2
	        		%		    scatter(JAICvec,Jtruvec,1,'r','x')
	        		%		elseif jj == 3
	        		%		    scatter(JMDLlvec,Jtruvec,1,'m','^')
	        		%		elseif jj == 4
	        		%			scatter(JMDLuvec,Jtruvec,1,'b','o') 
	        		%		elseif jj == 5
	        		%			scatter(JMDLentvec,Jtruvec,1,'k','x') 
	        		%		elseif jj == 6
	        		%			scatter(JMDLcountvec,Jtruvec,1,'r','^') 
	        		%		else
	        		%		end
					%		
					%		hold on
	        		%		%axis([-1 1 -1 1])
	        		%		axis([-0.2 0.2 -0.2 0.2])
	        		%		refline(1,0)
	        		%		refline
	        		%		
	        		%	end
	%
	        		%end
	
					% j = j+1;
		
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
end
	Jstor.TotJBICrrerr	= TotJBICrrerr	;
	Jstor.TotJAICrrerr	= TotJAICrrerr	;
	Jstor.TotJMDLlrrerr	= TotJMDLlrrerr	;
	Jstor.TotJMDLurrerr	= TotJMDLurrerr	;
	Jstor.TotJMDLentrrerr  = TotJMDLentrrerr  ;
	Jstor.TotJMDLcountrrerr= TotJMDLcountrrerr;

	Jstor.avgJBICrrerr		= mean(TotJBICrrerr)	;
	Jstor.avgJAICrrerr		= mean(TotJAICrrerr)	;
	Jstor.avgJMDLlrrerr		= mean(TotJMDLlrrerr)	;
	Jstor.avgJMDLurrerr		= mean(TotJMDLurrerr)	;
	Jstor.avgJMDLentrrerr  	= mean(TotJMDLentrrerr)  ;
	Jstor.avgJMDLcountrrerr	= mean(TotJMDLcountrrerr);






	Jstor2 = Jstor;