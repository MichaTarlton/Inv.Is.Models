%%modelgraphs.m
% For comparing model reconstructions
%

function Jscorestor = modelgraphs(AllStruct,name,time,Snamevec,beta,topdir,betadir,Tval,Nval)


scorestor = struct('Jy',{},'Jx1',{},'Jx2',{},'Jx3',{},'Jx4',{},'randos',{});

fn = fieldnames(AllStruct);

for i = 1:length(fn)
	    
	Jstor(i).name = fn(i);

AllStruct.(name).BLLH.perconerr.BIC      	= 
AllStruct.(name).BLLH.perconerr.AIC      	= 
AllStruct.(name).BLLH.perconerr.MDLl     	= 
AllStruct.(name).BLLH.perconerr.MDLu     	= 
AllStruct.(name).BLLH.perconerr.MDLent   	= 
AllStruct.(name).BLLH.perconerr.MDLcount 	= 

	%% loop for h figures
for hh = 1:4
    g = figure('Position', get(0, 'Screensize'));
        
    for hi = 1:length(Jstor)
        
        randos = Jstor(hi).randos;
        hy  = Jstor(hi).hy;
        hx1 = Jstor(hi).hx1;
        hx2 = Jstor(hi).hx2;
        hx3 = Jstor(hi).hx3;
        hx4 = Jstor(hi).hx4;
        

        %% Relative reconstruction error Nguyen 17, eq. 114
        % Implement theis elsewhere just coding out now
        % use to give a numerical value to how reliable inference method is
        % probably actually Jgrpahs3
        % sum from i where i<j ?? So essentially no self-interaction and no symmetrical interactions
        % by only using the upper triangle I should achieve the same thing
        Jstor(hi).(method(1)) = sqrt(sum((hx1-hy).^2 / sum(hy.^2)));
        Jstor(hi).(method(2)) = sqrt(sum((hx2-hy).^2 / sum(hy.^2)));
        Jstor(hi).(method(3)) = sqrt(sum((hx3-hy).^2 / sum(hy.^2)));
        Jstor(hi).(method(4)) = sqrt(sum((hx4-hy).^2 / sum(hy.^2)));


        %h = figure;
        
        subplot(3,3,hi)
        if hh == 1
            scatter(hx1,hy,1,'b','o')
        else if hh == 2
            scatter(hx2,hy,1,'r','x')
        else if hh == 3
            scatter(hx3,hy,1,'m','^')
        else if hh == 4
            scatter(hx4,hy,1,'b','o') 
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
        %if hi == 1
        %   xlabel('50')
        %   ylabel('3')
        %
        %elseif hi == 2
        %   xlabel('100')
        %   ylabel('4')
        %
        %elseif hi == 3
        %   xlabel('300')
        %   ylabel('5')
    end    
        axh = findobj(g,'Type','Axes');
        title(axh(9), ['N=',num2str(Nval(1))])
        ylabel(axh(9),['T=1E',num2str(log10(Tval(1)))])
        title(axh(8), ['N=',num2str(Nval(2))])
        ylabel(axh(6),['T=1E',num2str(log10(Tval(2)))])
        title(axh(7), ['N=',num2str(Nval(3))])
        ylabel(axh(3),['T=1E',num2str(log10(Tval(3)))])

        xlabel(axh(1),['\gamma_J = ',num2str(Jstor(9).(method(hh)))])
        xlabel(axh(2),['\gamma_J = ',num2str(Jstor(8).(method(hh)))])
        xlabel(axh(3),['\gamma_J = ',num2str(Jstor(7).(method(hh)))])
        xlabel(axh(4),['\gamma_J = ',num2str(Jstor(6).(method(hh)))])
        xlabel(axh(5),['\gamma_J = ',num2str(Jstor(5).(method(hh)))])
        xlabel(axh(6),['\gamma_J = ',num2str(Jstor(4).(method(hh)))])
        xlabel(axh(7),['\gamma_J = ',num2str(Jstor(3).(method(hh)))])
        xlabel(axh(8),['\gamma_J = ',num2str(Jstor(2).(method(hh)))])
        xlabel(axh(9),['\gamma_J = ',num2str(Jstor(1).(method(hh)))])
    
    
    
        
        % See the external script for this. Taken from the file exchange to add a super X and Y label to subplot
        % otherwise a very difficult task for an otherwise unimportant aesthetic
        [ax1,h1]=suplabel('Generated through inference methods');
        [ax2,h2]=suplabel('Monte Carlo','y');
        sgtitle({'Comparative h values using ',method(hh),' inference method for \beta = ',num2str(beta)})
end