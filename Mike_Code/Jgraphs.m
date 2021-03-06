%% JGraphs.m
% For comparing inferred J values


function JGraphs(JHstruct,sanity,N,T,name,time)%N,jn,h_on,sparsity,time,T,name)

% First for Js
Jy = [JHstruct.Jsparse];
Jx1 = [sanity.mfJ]; %mean field
Jx2 = [sanity.tapJ]; %tap

sname = inputname(1);

if strcmpi(sname, 'JHnorm')
	distname = 'Normal Distribution'
elseif strcmpi(sname, 'JHdiscon')
	distname = 'Disconnected'
elseif strcmpi(sname,'JHdimer')
	distname = 'Independent Pairs'
elseif strcmpi(sname,'JHferr')
	distname = 'Ferromagnetic'
else
	distname = ['Unknown: ',sname]
end

h = figure;

subplot(1,2,1)
scatter(Jx1(:),Jy(:),[],'b')
%axis([-1 1 -1 1])
refline(1,0) 
hold on 
scatter(Jx2(:),Jy(:),[],'r')
title({'J values: infered v real',['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
xlabel('Infered - Blue:nMF, Red:TAP')
ylabel('Real')

%% For H values

 hy = [JHstruct.Hsparse];
hx1 = [sanity.mfh]; %mean field
hx2 = [sanity.taph]; %tap

%figure
subplot(1,2,2)
scatter(hx1(:),hy(:),[],'b')
%axis([-1 1 -1 1])
refline(1,0) 
hold on 
scatter(hx2(:),hy(:),[],'r')
title({'h values: infered v real',['N = ',num2str(N)],['T = 1E',num2str(log10(T))]})
xlabel('Infered - Blue:nMF, Red:TAP')
ylabel('Real')
sgtitle(distname)

savefig(h,[name,'\',time(1:5),'JGraphs',sname,'_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.fig'],'compact')
saveas(h,[time(1:5),'JGraph',sname,'_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
%saveas(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
%save([name,'\',time(1:5),'JHDstruct_N',num2str(N),'_T1E',num2str(log10(T)),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'JHdimer');	