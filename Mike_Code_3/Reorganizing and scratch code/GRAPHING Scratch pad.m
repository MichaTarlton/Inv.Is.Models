%%% GRAPHING Scratch pad


%%rgb converter

rgb = []

%matrgb = rgb./255;

% pulled from a an online palette
c1 = [249, 65, 68]./255;
c2 = [243, 114, 44]./255;
c3 = [248, 150, 30]./255;
c4 = [249, 199, 79]./255;
c5 = [144, 190, 109]./255;
c6 = [67, 170, 139]./255;
c7 = [87, 117, 144]./255;

colorvec = {c1,c2,c3,c4,c5,c6,c7};

plot(xxaxisvec,yyaxisvec, 'Color', topocolor, 'LineStyle', '-o')



toponamevec = [{'Cayley Tree'},{'Fully Connected'},{'Independent Pair'},{'2D Ising'},{'Erdős–Rényi'}];

set(gca,'Color','k')

set(gca,'Color',cbg)

% anice flat black/grey
cbg = [53, 53, 53]./255

% a light grey
clg = [217, 217, 217]./255

% For deciding output
if comp == 1	
	
elseif comp == 2
	
elseif comp == 3 

elseif comp == 4

end


%% from jgraphs4
 if jj == 1
            scatter(Jx1(randos),Jy(randos),1,'b','o')
        else if jj == 2
            scatter(Jx2(randos),Jy(randos),1,'r','x')
        else if jj == 3
            scatter(Jx3(randos),Jy(randos),1,'m','^')
        else if jj == 4
        	scatter(Jx4(randos),Jy(randos),1,'b','o') 
        else
                end
            end
            end
        end
        %scatter(Jx1(:),Jy(:),[],'.','b')
        %axis([-1 1 -1 1])
        axis([-0.2 0.2 -0.2 0.2])
        refline(1,0)
        refline
        hold on
        
        %ylabel('nMF')
        %title({'J values'})
        %if j == 1
        %	xlabel('50')
        %	ylabel('3')
        %
        %elseif j == 2
        %	xlabel('100')
        %	ylabel('4')
        %
        %elseif j == 3
        %	xlabel('300')
        %	ylabel('5')
    end    
        ax = findobj(h,'Type','Axes');
        
        title(ax(9), ['N=',num2str(Nval(1))])
        ylabel(ax(9),['T=1E',num2str(log10(Tval(1)))])
        title(ax(8), ['N=',num2str(Nval(2))])
        ylabel(ax(6),['T=1E',num2str(log10(Tval(2)))])
        title(ax(7), ['N=',num2str(Nval(3))])
        ylabel(ax(3),['T=1E',num2str(log10(Tval(3)))])


        %title(ax(9),'N = 50')
        %ylabel(ax(9),'T=1E3')
        %title(ax(8),'N = 100')
        %ylabel(ax(6),'T=1E4')
        %title(ax(7),'N = 300')
        %ylabel(ax(3),'T=1E5')

        xlabel(ax(1),['\gamma_J = ',num2str(Jstor(9).(method(jj)))])
        xlabel(ax(2),['\gamma_J = ',num2str(Jstor(8).(method(jj)))])
        xlabel(ax(3),['\gamma_J = ',num2str(Jstor(7).(method(jj)))])
        xlabel(ax(4),['\gamma_J = ',num2str(Jstor(6).(method(jj)))])
        xlabel(ax(5),['\gamma_J = ',num2str(Jstor(5).(method(jj)))])
        xlabel(ax(6),['\gamma_J = ',num2str(Jstor(4).(method(jj)))])
        xlabel(ax(7),['\gamma_J = ',num2str(Jstor(3).(method(jj)))])
        xlabel(ax(8),['\gamma_J = ',num2str(Jstor(2).(method(jj)))])
        xlabel(ax(9),['\gamma_J = ',num2str(Jstor(1).(method(jj)))])
    
    
    
        
        % See the external script for this. Taken from the file exchange to add a super X and Y label to subplot
        % otherwise a very difficult task for an otherwise unimportant aesthetic
        [ax1,h1]=suplabel('Generated through inference methods');
        [ax2,h2]=suplabel('Monte Carlo','y');
        sgtitle({'Comparative J values using ',method(jj),' inference method for \beta = ',num2str(beta)})
    
savefig(h,[topdir,'\',betadir,'\',time(1:5),'MultiGrpah','_beta',num2str(beta),'_',char(method(jj)),'_',time(6:12),'.fig'],'compact')
saveas(h,[topdir,'\',time(1:5),'JGraphs','_beta',num2str(beta),'_',char(method(jj)),'_',time(6:12),'.png'])
