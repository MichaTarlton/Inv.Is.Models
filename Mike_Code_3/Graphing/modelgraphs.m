%%modelgraphs.m
% For comparing model reconstructions
%

% function Jscorestor = modelgraphs(AllStruct,time,Snamevec,beta,topdir,betadir,Tval,Nval)
function Jscorestor = modelgraphs(OverStruct,sparsity,betavec,Nvec,Tvec,topologies)

fn = fieldnames(OverStruct);

for i = 1:length(fn) % for each set of N and T across the same beta
    name = fn{i};
    
for topo = 1:length(topology)
	Jstor(topo).name = name;
    Jstor(topo).topo = topo;
    % Jstor(i).beta 	= AllStruct.(name). not stored
	Jstor(topo).N    	= AllStruct.(name).N;
	Jstor(topo).T    	= AllStruct.(name).T;
	% Jstor(i).spars  = AllStruct.(name). not stored
	% Jstor(i).c    	= AllStruct.(name).c % coordiantion number 
	% Jstor(i).T    	= AllStruct.(name).couplings % coupling number 

	Jstor(topo).pcBIC      = AllStruct.(name).BLLH(topo).perconerr.BIC      	;
	Jstor(topo).pcAIC      = AllStruct.(name).BLLH(topo).perconerr.AIC      	;
	Jstor(topo).pcMDLl     = AllStruct.(name).BLLH(topo).perconerr.MDLl     	;
	Jstor(topo).pcMDLu     = AllStruct.(name).BLLH(topo).perconerr.MDLu     	;
	Jstor(topo).pcMDLent   = AllStruct.(name).BLLH(topo).perconerr.MDLent   	;
	Jstor(topo).pcMDLcount = AllStruct.(name).BLLH(topo).perconerr.MDLcount 	;
% https://www.mathworks.com/matlabcentral/answers/414989-logical-indexing-with-a-structure
	
end
end

h = figure%('Position', get(0, 'Screensize'));
hold on
plot([Jstor.T], [Jstor.pcBIC]      ,'b-o')
plot([Jstor.T], [Jstor.pcAIC]      ,'r-o')
plot([Jstor.T], [Jstor.pcMDLl]     ,'m-o')
plot([Jstor.T], [Jstor.pcMDLu]     ,'g-o')
plot([Jstor.T], [Jstor.pcMDLent]   ,'y-o')
plot([Jstor.T], [Jstor.pcMDLcount] ,'k-o')



%% psuedo coding
% On each mini subplot I want to compare error to sparsity (maybe just set it up to slot in anything) for multiple models
% we can set this inside of a nest of the other variables we want, but those two will be in theory out basic subcore

Jstor = OverStruct(1).Jstor



% plot size
figwidth = length(betavec)
figheight = length(Tvec)


for NT = 1:length(Nvec)

j = 1;

	for topology = 1:length(topology) % for each complete figure
		h = figure;
		hold on
	
		for bt = 1:length(betavec)
				
			for Tt = 1:length(Tvec)
														
				BICespc		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).BICespc		;
				AICespc		 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).AICespc		;
				MDLlespc	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLlespc	 	;
				MDLuespc	 = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLuespc 	;
				MDLentespc   = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLentespc   ;
				MDLcountespc = Jstor(topology).Topo(bt).beta(Nt).N(Tt).MDLcountespc ;

				
				subplot(beta,T,j)
				plot([sparsity], [BICespc		],'b-o')
				plot([sparsity], [AICespc		],'r-o')
				plot([sparsity], [MDLlespc		],'m-o')
				plot([sparsity], [MDLuespc		],'g-o')
				plot([sparsity], [MDLentespc  	],'y-o')
				plot([sparsity], [MDLcountespc	],'k-o')
			
			j = j+1;
	
			end
	
		end
	
	end

end






















%% method = ["nMF","TAP","PLMF","PLLH"];
%
%for jj = 1:4 % four different figures, one for each method
%    h = figure('Position', get(0, 'Screensize'));
%        
%    for j = 1:length(Jstor)
%            
%        %h = figure;
%        
%        subplot(2,3,j)
%
%        scatter(Jstor.pcBIC     , Jstor.T,1,'b','o')
%        scatter(Jstor.pcAIC     , Jstor.T,1,'r','o')
%        scatter(Jstor.pcMDLl    , Jstor.T,1,'m','o')
%        scatter(Jstor.pcMDLu    , Jstor.T,1,'g','o')
%        scatter(Jstor.pcMDLent  , Jstor.T,1,'l','o')
%        scatter(Jstor.pcMDLcount, Jstor.T,1,'k','o')
%
%
%        %if jj == 1
%        %    scatter(Jx1(randos),Jy(randos),1,'b','o')
%        %else if jj == 2
%        %    scatter(Jx2(randos),Jy(randos),1,'r','x')
%        %else if jj == 3
%        %    scatter(Jx3(randos),Jy(randos),1,'m','^')
%        %else if jj == 4
%        %	scatter(Jx4(randos),Jy(randos),1,'b','o') 
%        %else
%                end
%            end
%            end
%        end
%        %scatter(Jx1(:),Jy(:),[],'.','b')
%        %axis([-1 1 -1 1])
%        %axis([-0.2 0.2 -0.2 0.2])
%        %refline(1,0)
%        %refline
%        %hold on
%        
%        %ylabel('nMF')
%        %title({'J values'})
%        %if j == 1
%        %	xlabel('50')
%        %	ylabel('3')
%        %
%        %elseif j == 2
%        %	xlabel('100')
%        %	ylabel('4')
%        %
%        %elseif j == 3
%        %	xlabel('300')
%        %	ylabel('5')
%    end    
%        ax = findobj(h,'Type','Axes');
%        
%        title(ax(9), ['N=',num2str(Nval(1))])
%        ylabel(ax(9),['T=1E',num2str(log10(Tval(1)))])
%        title(ax(8), ['N=',num2str(Nval(2))])
%        ylabel(ax(6),['T=1E',num2str(log10(Tval(2)))])
%        title(ax(7), ['N=',num2str(Nval(3))])
%        ylabel(ax(3),['T=1E',num2str(log10(Tval(3)))])


        %title(ax(9),'N = 50')
        %ylabel(ax(9),'T=1E3')
        %title(ax(8),'N = 100')
        %ylabel(ax(6),'T=1E4')
        %title(ax(7),'N = 300')
        %ylabel(ax(3),'T=1E5')

%        xlabel(ax(1),['\gamma_J = ',num2str(Jstor(9).(method(jj)))])
%        xlabel(ax(2),['\gamma_J = ',num2str(Jstor(8).(method(jj)))])
%        xlabel(ax(3),['\gamma_J = ',num2str(Jstor(7).(method(jj)))])
%        xlabel(ax(4),['\gamma_J = ',num2str(Jstor(6).(method(jj)))])
%        xlabel(ax(5),['\gamma_J = ',num2str(Jstor(5).(method(jj)))])
%        xlabel(ax(6),['\gamma_J = ',num2str(Jstor(4).(method(jj)))])
%        xlabel(ax(7),['\gamma_J = ',num2str(Jstor(3).(method(jj)))])
%        xlabel(ax(8),['\gamma_J = ',num2str(Jstor(2).(method(jj)))])
%        xlabel(ax(9),['\gamma_J = ',num2str(Jstor(1).(method(jj)))])
%    
%    
    
        
%        % See the external script for this. Taken from the file exchange to add a super X and Y label to subplot
%        % otherwise a very difficult task for an otherwise unimportant aesthetic
%        [ax1,h1]=suplabel('Generated through inference methods');
%        [ax2,h2]=suplabel('Monte Carlo','y');
%        sgtitle({'Comparative J values using ',method(jj),' inference method for \beta = ',num2str(beta)})
%    
%savefig(h,[topdir,'\',betadir,'\',time(1:5),'JGraphs','_beta',num2str(beta),'_',char(method(jj)),'_',time(6:12),'.fig'],'compact')
%saveas(h,[topdir,'\',time(1:5),'JGraphs','_beta',num2str(beta),'_',char(method(jj)),'_',time(6:12),'.png'])
%
%end
end