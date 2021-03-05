%% JGraphs.m
% For comparing inferred J values


function JGraphs(JHstruct,sanity,N,T,name,time)%N,jn,h_on,sparsity,time,T,name)

JHstruct = JHstruct(1)
sanity = sanity(1)

% First for Js
% note that you can get rid of the  upper tri stuff now, also the brackets
Jypre =  JHstruct.Jsparse;
Jx1pre = sanity.mfJ; %mean field
Jx2pre = sanity.tapJ; %tap
Jx3pre = sanity.plJmf;

hy = [JHstruct.Hsparse];
hx1 = [sanity.mfh]; %mean field
hx2 = [sanity.taph]; %tap
hx3 = [sanity.plhmf];

% to create vector out of upper triangle
truth = triu(true(size(Jypre)),1);
Jy  = Jypre(truth);
Jx1 = Jx1pre(truth);
Jx2 = Jx2pre(truth);
Jx3 = Jx3pre(truth);

randos = randperm(N^2/2-N/2,floor(N^2/3)); %| if this breaks, it's bc someone use an uneven N value

sname = inputname(1);

label = 'Infered - Blue:nMF, Red:TAP';

if strcmpi(sname, 'JHnorm')
	distname = 'Normal Distribution'
elseif strcmpi(sname, 'JHdiscon')
	distname = 'Disconnected'
elseif strcmpi(sname,'JHdimer')
	distname = 'Independent Pairs'
	label = 'Infered for Independent Pairs'
	Jx1 = [sanity.Jpair]; %mean field
	Jx2 = []; %tap
	hx1 = [sanity.hpair]; %mean field
	hx2 = []; %tap
elseif strcmpi(sname,'JHferr')
	distname = 'Ferromagnetic'
else
	distname = ['Unknown: ',sname]
end

h = figure;


subplot(3,2,1)
scatter(Jx1(randos),Jy(randos),1,'.','b')
%scatter(Jx1(:),Jy(:),[],'.','b')
axis([-1 1 -1 1])
refline(1,0) 
ylabel('nMF')
hold on
title({'J values'})
%if ~strcmpi(sname, 'JHdimer')
%scatter(Jx2(randos),Jy(randos),1,'.','r')
%end

subplot(3,2,3)
scatter(Jx2(randos),Jy(randos),1,'.','r')
%scatter(Jx1(:),Jy(:),[],'.','b')
axis([-1 1 -1 1])
refline(1,0)
ylabel('TAP') 
hold on

subplot(3,2,5)
scatter(Jx3(randos),Jy(randos),1,'.','m')
%scatter(Jx1(:),Jy(:),[],'.','b')
axis([-1 1 -1 1])
refline(1,0) 
ylabel('PL-MF')
hold on

%title({'J values: infered v real',['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
%xlabel(label)
%ylabel('Real')






%% For H values

%figure
subplot(3,2,2)

scatter(hx1(:),hy(:),[],'b')
axis([-1 1 -1 1])
refline(1,0) 
hold on
title({'h values'})
%if ~strcmpi(sname, 'JHdimer')
%scatter(hx1(:),hy(:),[],'r')
%end

subplot(3,2,4)
scatter(hx2(:),hy(:),[],'b')
axis([-1 1 -1 1])
refline(1,0) 
hold on
%if ~strcmpi(sname, 'JHdimer')
%scatter(hx2(:),hy(:),[],'r')
%end

subplot(3,2,6)
scatter(hx3(:),hy(:),[],'b')
axis([-1 1 -1 1])
refline(1,0) 
hold on
%if ~strcmpi(sname, 'JHdimer')
%scatter(hx2(:),hy(:),[],'r')
%end
sgtitle({distname,['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})

figure 
title({'Comparison of J values'})
scatter(Jx1(randos),Jx2(randos),10,'.','b')
hold on
scatter(Jx1(randos),Jx3(randos),10,'o','r')
scatter(Jx2(randos),Jx3(randos),10,'d','m')
refline(1,0)

%title({'h values: infered v real',['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
%xlabel(label)
%ylabel('Real')






%savefig(h,[name,'\',time(1:5),'JGraphs',sname,'_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.fig'],'compact')
%saveas(h,[time(1:5),'JGraph',sname,'_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])

%saveas(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
%save([name,'\',time(1:5),'JHDstruct_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'JHdimer');	