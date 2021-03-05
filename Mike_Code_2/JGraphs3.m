%% JGraphs.m
% For comparing inferred J values
% Point of this variant is to make the described graph as Yasser wants it
%


function Jstor = JGraphs3(AllStruct,time,Snamevec,beta,topdir,betadir,Tval,Nval)%N,jn,h_on,sparsity,time,T,name)

%PLLHout(1).J,PLLHout(2).J,PLLHout(3).J,PLLHout(4).J,PLLHout(5).J,PLLHout(6).J,PLLHout(7).J,PLLHout(8).J,PLLHout(9).J,

% Jstor = JGraphs2(JHN50T3,JHN50T4,JHN50T5,JHN100T3,JHN100T4,JHN100T5,JHN300T3,JHN300T4,JHN300T5,...
%     PLLHout.J(1),PLLHout.J(2),PLLHout.J(3),PLLHout.J(4),PLLHout.J(5),PLLHout.J(6),PLLHout.J(7),PLLHout.J(8),PLLHout.J(9),...
%     time,name,beta);

%JHstructvec = [JHN50T3,JHN50T4,JHN50T5,JHN100T3,JHN100T4,JHN100T5,JHN300T3,JHN300T4,JHN300T5];  %not sure I can do this as it's an array of structs. it works just not how I thought
%Sstructvec = [SN50T3,SN50T4,SN50T5,SN100T3,SN100T4,SN100T5,SN300T3,SN300T4,SN300T5];
%JHnamevec = {'JHN50T3','JHN50T4','JHN50T5','JHN100T3','JHN100T4','JHN100T5','JHN300T3','JHN300T4','JHN300T5'};
%Snamevec = {"SN50T3","SN50T4","SN50T5","SN100T3","SN100T4","SN100T5","SN300T3","SN300T4","SN300T5"};
%sanityvec = [sanityN50T3,sanityN50T4,sanityN50T5,sanityN100T3,sanityN100T4,sanityN100T5,sanityN300T3,sanityN300T4,sanityN300T5];
%pllhvec = {PLLHout(1).J,PLLHout(2).J,PLLHout(3).J,PLLHout(4).J,PLLHout(5).J,PLLHout(6).J,PLLHout(7).J,PLLHout(8).J,PLLHout(9).J}; %|Using cell structure for PLLH single struct method, which is probably what should be standardized at this point
%Jstor = struct('JHN50T3',{},'JHN50T4',{},'JHN50T5',{},'JHN100T3',{},'JHN100T4',{},'JHN100T5',{},'JHN300T3',{},'JHN300T4',{},'JHN300T5',{}); %| to store the plot points

Jstor = struct('Jy',{},'Jx1',{},'Jx2',{},'Jx3',{},'Jx4',{},'randos',{});

fn = fieldnames(AllStruct); %for iterating over the fields in the struct

% LT = length(T); %for setting dimensions of figure
% LN = length(N); %for setting dimensions of figure

for i = 1:length(fn)
	    
	Jstor(i).name = fn(i);
	
	%JHstruct = JHstructvec(i);
	%sanity = sanityvec(i); %suppressing for PLLH 
	
	% First for Js
	
	% note that you can get rid of the  upper tri stuff now, also the brackets
	Jypre =  AllStruct.(fn{i}).Jtru; %JHstruct.Jsparse;
	Jx1pre = AllStruct.(fn{i}).Jmf; %sanity.mfJ; %mean field
	Jx2pre = AllStruct.(fn{i}).Jtap; %sanity.tapJ; %tap
	Jx3pre = AllStruct.(fn{i}).Jplmf; %sanity.plJmf;
	Jx4pre = AllStruct.(fn{i}).Jpllh; %PLLH(i).J;
	
	N = size(Jypre,2);
	
	% to create vector out of upper triangle
	truth = triu(true(size(Jypre)),1);
	Jstor(i).Jy  = Jypre(truth);
	Jstor(i).Jx1 = Jx1pre(truth);
	Jstor(i).Jx2 = Jx2pre(truth);
	Jstor(i).Jx3 = Jx3pre(truth);
	Jstor(i).Jx4 = Jx4pre(truth);
	
	rando = randperm(N^2/2-N/2,floor(N^2/3)); %| if this breaks, it's bc someone use an uneven N value. This is for creating a random selection of the pint to be graphed due to issues when using all points
	Jstor(i).randos = rando;
	  
    Jstor(i).hy =  AllStruct.(fn{i}).htru; 
    Jstor(i).hx1 = AllStruct.(fn{i}).hmf; %
    Jstor(i).hx2 = AllStruct.(fn{i}).htap'; 
    Jstor(i).hx3 = AllStruct.(fn{i}).hplmf;
    Jstor(i).hx4 = AllStruct.(fn{i}).hpllh';  

disp(['End Jh store ',Jstor(i).name])
end

method = ["nMF","TAP","PLMF","PLLH"];

for jj = 1:4
    h = figure('Position', get(0, 'Screensize'));
        
    for j = 1:length(Jstor)
        
        randos = Jstor(j).randos;
        Jy  = Jstor(j).Jy;
        Jx1 = Jstor(j).Jx1;
        Jx2 = Jstor(j).Jx2;
        Jx3 = Jstor(j).Jx3;
        Jx4 = Jstor(j).Jx4;
        

        %% Relative reconstruction error Nguyen 17, eq. 114
        % Implement theis elsewhere just coding out now
        % use to give a numerical value to how reliable inference method is
        % probably actually Jgrpahs3
        % sum from i where i<j ?? So essentially no self-interaction and no symmetrical interactions
        % by only using the upper triangle I should achieve the same thing
        Jstor(j).(method(1)) = sqrt(sum((Jx1-Jy).^2 / sum(Jy.^2)));
        Jstor(j).(method(2)) = sqrt(sum((Jx2-Jy).^2 / sum(Jy.^2)));
        Jstor(j).(method(3)) = sqrt(sum((Jx3-Jy).^2 / sum(Jy.^2)));
        Jstor(j).(method(4)) = sqrt(sum((Jx4-Jy).^2 / sum(Jy.^2)));


        %h = figure;
        
        subplot(3,3,j)

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
    
savefig(h,[topdir,'\',betadir,'\',time(1:5),'JGraphs','_beta',num2str(beta),'_',char(method(jj)),'_',time(6:12),'.fig'],'compact')
saveas(h,[topdir,'\',time(1:5),'JGraphs','_beta',num2str(beta),'_',char(method(jj)),'_',time(6:12),'.png'])

end


%h = xlabel(1,'50')
%h = ylabel(1,'3')
%h = xlabel(2,'100')
%h = ylabel(2,'4')
%h = xlabel(3,'300')
%h = ylabel(3,'5')
%
%
%ax = findobj(h,'Type','Axes');
%for k=1:length(ax)
%    ylabel(ax(k),{'Nice'})
%    title(ax(k),{'Very Nice'})
%end

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
    
savefig(g,[topdir,'\',betadir,'\',time(1:5),'hGraphs','_beta',num2str(beta),'_',char(method(hh)),'_',time(6:12),'.fig'],'compact')
saveas(g,[topdir,'\',time(1:5),'hGraphs','_beta',num2str(beta),'_',char(method(hh)),'_',time(6:12),'.png'])

end



end

%if ~strcmpi(sname, 'JHdimer')
%scatter(Jx2(randos),Jy(randos),1,'.','r')
%end

%subplot(3,3,2)
%scatter(Jx2(randos),Jy(randos),1,'.','r')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('TAP')
%hold on
%
%subplot(3,3,3)
%scatter(Jx3(randos),Jy(randos),1,'.','m')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('PL-MF')
%hold on
%
%subplot(3,3,4)
%scatter(Jx1(randos),Jy(randos),1,'.','b')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('nMF')
%hold on
%title({'J values'})
%%if ~strcmpi(sname, 'JHdimer')
%%scatter(Jx2(randos),Jy(randos),1,'.','r')
%%end
%
%subplot(3,3,5)
%scatter(Jx2(randos),Jy(randos),1,'.','r')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('TAP')
%hold on
%
%subplot(3,3,6)
%scatter(Jx3(randos),Jy(randos),1,'.','m')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('PL-MF')
%hold on
%
%subplot(3,3,7)
%scatter(Jx1(randos),Jy(randos),1,'.','b')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('nMF')
%hold on
%title({'J values'})
%%if ~strcmpi(sname, 'JHdimer')
%%scatter(Jx2(randos),Jy(randos),1,'.','r')
%%end
%
%subplot(3,3,8)
%scatter(Jx2(randos),Jy(randos),1,'.','r')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('TAP')
%hold on
%
%subplot(3,3,9)
%scatter(Jx3(randos),Jy(randos),1,'.','m')
%%scatter(Jx1(:),Jy(:),[],'.','b')
%axis([-1 1 -1 1])
%refline(1,0)
%ylabel('PL-MF')
%hold on
%
%%title({'J values: infered v real',['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
%%xlabel(label)
%%ylabel('Real')
%
%
%sgtitle({distname,['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
%
%figure
%title({'Comparison of J values'})
%scatter(Jx1(randos),Jx2(randos),10,'.','b')
%hold on
%scatter(Jx1(randos),Jx3(randos),10,'o','r')
%scatter(Jx2(randos),Jx3(randos),10,'d','m')
%refline(1,0)

%title({'h values: infered v real',['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
%xlabel(label)
%ylabel('Real')

%hy = [JHstruct.Hsparse];
%hx1 = [sanity.mfh]; %mean field
%hx2 = [sanity.taph]; %tap
%hx3 = [sanity.plhmf];


%%% For H values
%
%%figure
%subplot(3,2,2)
%
%scatter(hx1(:),hy(:),[],'b')
%axis([-1 1 -1 1])
%refline(1,0)
%hold on
%title({'h values'})
%%if ~strcmpi(sname, 'JHdimer')
%%scatter(hx1(:),hy(:),[],'r')
%%end
%
%subplot(3,2,4)
%scatter(hx2(:),hy(:),[],'b')
%axis([-1 1 -1 1])
%refline(1,0)
%hold on
%%if ~strcmpi(sname, 'JHdimer')
%%scatter(hx2(:),hy(:),[],'r')
%%end
%
%subplot(3,2,6)
%scatter(hx3(:),hy(:),[],'b')
%axis([-1 1 -1 1])
%refline(1,0)
%hold on
%%if ~strcmpi(sname, 'JHdimer')
%%scatter(hx2(:),hy(:),[],'r')
%%end




%saveas(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
%save([name,'\',time(1:5),'JHDstruct_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'JHdimer');






