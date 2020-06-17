%% Graphs.m
% For Correlation and magnetization sanity checks


%function Graphs(SstructDisc,SstructFerr,sanitydimer,sanitydisc,sanityferr,N,T,name,time)%N,jn,h_on,sparsity,time,T,name)
function Graphs(SstructDisc,sanitydimer,sanitydisc,N,T,name,time) % Minus sanityferr

mx = [SstructDisc.mfinal];
my = [sanitydisc.th];
%fmx = [SstructFerr.mfinal];
%fmy = [sanityferr.th];
mtchk=mean([sanitydisc.mtchk]);
vtchk=var([sanitydisc.tchk]);

%h = figure;
h = figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
scatter(mx(:),my(:),[],'b')
%axis([-1 1 -1 1])

refline(1,0) 
hold on 


%scatter(fmx(:),fmy(:),[],'r')
title({'Disconnected J Magnetizations: mi',['N = ',num2str(N)],['T = 1E',num2str(log10(T))],['Variance for mi (disconnected) = ',num2str(vtchk)]})
xlabel('Generated from Monte Carlo')
ylabel('Calcualted from mean field methods. Blue: Disconnected') % Red: Ferromagnetic')

%Correlation

Cx = [sanitydimer.Cij];
Cx2 = Cx - diag(diag(Cx));
%Cy = [sanitydimer.Cpair];
%Cy1 = [sanitydimer.Cpair1].*1e-4; %seems to normalize a bit
%Cy2 = [sanitydimer.Cpair2];
Cy3 = [sanitydimer.Cpair3];

%figure
subplot(1,2,2)
%scatter(Cx(:),Cy1(:),[],'b')
%axis([-1 1 -1 1])


hold on 
%scatter(Cx(:),Cy2(:),[],'r')
scatter(Cx2(:),Cy3(:),[],'m')
refline
refline(1,0)
%FitC = polyfit(Cx2(:),Cy3(:),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
%plot(polyval(FitC,[-1:N:1]));


title({'Correlations for Independent Pairs',['N = ',num2str(N)],['T = 1E',num2str(log10(T))],['Avg. difference between M values (blue) = ',num2str(mean(abs(Cx2(:)-Cy3(:)),'all'))]})
xlabel('Monte Carlo Generated Correlations')
ylabel('Correlation: C = tanh(J)')

savefig(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.fig'],'compact')
saveas(h,[time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
%saveas(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
