%% JGraphs.m
% For comparing inferred J values
% Point of this variant is to make the described graph as Yasser wants it
%


function Jstor = JGraphs2(JHN50T3,JHN50T4,JHN50T5,JHN100T3,JHN100T4,JHN100T5,JHN300T3,JHN300T4,JHN300T5,...
    sanityN50T3,sanityN50T4,sanityN50T5,sanityN100T3,sanityN100T4,sanityN100T5,sanityN300T3,sanityN300T4,sanityN300T5,...
    PLLH,...
    time,name,beta)%N,jn,h_on,sparsity,time,T,name)

%PLLHout(1).J,PLLHout(2).J,PLLHout(3).J,PLLHout(4).J,PLLHout(5).J,PLLHout(6).J,PLLHout(7).J,PLLHout(8).J,PLLHout(9).J,


% Jstor = JGraphs2(JHN50T3,JHN50T4,JHN50T5,JHN100T3,JHN100T4,JHN100T5,JHN300T3,JHN300T4,JHN300T5,...
%     PLLHout.J(1),PLLHout.J(2),PLLHout.J(3),PLLHout.J(4),PLLHout.J(5),PLLHout.J(6),PLLHout.J(7),PLLHout.J(8),PLLHout.J(9),...
%     time,name,beta);




JHstructvec = [JHN50T3,JHN50T4,JHN50T5,JHN100T3,JHN100T4,JHN100T5,JHN300T3,JHN300T4,JHN300T5];  %not sure I can do this as it's an array of structs. it works just not how I thought
%Sstructvec = [SN50T3,SN50T4,SN50T5,SN100T3,SN100T4,SN100T5,SN300T3,SN300T4,SN300T5];
JHnamevec = {'JHN50T3','JHN50T4','JHN50T5','JHN100T3','JHN100T4','JHN100T5','JHN300T3','JHN300T4','JHN300T5'};
%Snamevec = {"SN50T3","SN50T4","SN50T5","SN100T3","SN100T4","SN100T5","SN300T3","SN300T4","SN300T5"};
sanityvec = [sanityN50T3,sanityN50T4,sanityN50T5,sanityN100T3,sanityN100T4,sanityN100T5,sanityN300T3,sanityN300T4,sanityN300T5];
%pllhvec = {PLLHout(1).J,PLLHout(2).J,PLLHout(3).J,PLLHout(4).J,PLLHout(5).J,PLLHout(6).J,PLLHout(7).J,PLLHout(8).J,PLLHout(9).J}; %|Using cell structure for PLLH single struct method, which is probably what should be standardized at this point
%Jstor = struct('JHN50T3',{},'JHN50T4',{},'JHN50T5',{},'JHN100T3',{},'JHN100T4',{},'JHN100T5',{},'JHN300T3',{},'JHN300T4',{},'JHN300T5',{}); %| to store the plot points
Jstor = struct('Jy',{},'Jx1',{},'Jx2',{},'Jx3',{},'randos',{});
sname = 'combine J';

for i = 1:length(JHstructvec)
    
    Jstor(i).name =JHnamevec(i);
    
    JHstruct = JHstructvec(i);
    sanity = sanityvec(i); %suppressing for PLLH 
    
    % First for Js
    
    % note that you can get rid of the  upper tri stuff now, also the brackets
    Jypre =  JHstruct.Jsparse;
    Jx1pre = sanity.mfJ; %mean field
    Jx2pre = sanity.tapJ; %tap
    Jx3pre = sanity.plJmf;
    Jx4pre = PLLH(i).J;
    
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
    
end

method = ["nMF","TAP","PL-MF","PLLH"]

for jj = 1:4
    h = figure('Position', get(0, 'Screensize'));
    
    
    for j = 1:length(Jstor)
        
        randos = Jstor(j).randos;
        Jy  = Jstor(j).Jy;
        Jx1 = Jstor(j).Jx1;
        Jx2 = Jstor(j).Jx2;
        Jx3 = Jstor(j).Jx3;
        Jx4 = Jstor(j).Jx4;
        
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
        title(ax(9),'50')
        ylabel(ax(9),'3')
        title(ax(8),'100')
        ylabel(ax(6),'4')
        title(ax(7),'300')
        ylabel(ax(3),'5')
        
        % See the external script for this. Taken from the file exchange to add a super X and Y label to subplot
        % otherwise a very difficult task for an otherwise unimportant aesthetic
        [ax1,h1]=suplabel('Generated through inferrence methods');
        [ax2,h2]=suplabel('Monte Carlo','y');
        sgtitle({'Comparative J values using ',method(jj),' inference method for \beta = ',num2str(beta)})
    
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



savefig(h,[name,'\',time(1:5),'JGraphs',sname,'_beta',num2str(beta),'_',time(6:12),'.fig'],'compact')
saveas(h,[time(1:5),'JGraph',sname,'_beta',num2str(beta),'_',time(6:12),'.png'])

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






