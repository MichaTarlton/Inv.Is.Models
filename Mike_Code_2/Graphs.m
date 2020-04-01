%% Graphs.m
% For Correlation and magnetization sanity checks


%function Graphs(SstructDisc,SstructFerr,sanitydimer,sanitydisc,sanityferr,N,T,name,time)%N,jn,h_on,sparsity,time,T,name)
function Graphs(SstructDisc,sanitydimer,sanitydisc,N,T,name,time) % Minus sanityferr

mx = [SstructDisc.mfinal];
my = [sanitydisc.th];
%fmx = [SstructFerr.mfinal];
%fmy = [sanityferr.th];
mtchk=mean([sanitydisc.mtchk]);

%h = figure;
h = figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
scatter(mx(:),my(:),[],'b')
%axis([-1 1 -1 1])

refline(1,0) 
hold on 
%scatter(fmx(:),fmy(:),[],'r')
title({'Disconnected J Magnetizations: mi',['N = ',num2str(N)],['T = 1E',num2str(log10(T))],['Avg. difference between M values (disconnected) = ',num2str(mtchk)]})
xlabel('Generated from Monte Carlo')
ylabel('Calcualted from mean field methods. Blue: Disconnected') % Red: Ferromagnetic')

%Correlation

Cx = [sanitydimer.Cij];
%Cy = [sanitydimer.Cpair];
%Cy1 = [sanitydimer.Cpair1].*1e-4; %seems to normalize a bit
%Cy2 = [sanitydimer.Cpair2];
Cy3 = [sanitydimer.Cpair3];

%figure
subplot(1,2,2)
%scatter(Cx(:),Cy1(:),[],'b')
%axis([-1 1 -1 1])

refline(1,0)
hold on 
%scatter(Cx(:),Cy2(:),[],'r')
scatter(Cx(:),Cy3(:),[],'m')

title({'Correlations for Dimers',['N = ',num2str(N)],['T = 1E',num2str(log10(T))],['Avg. difference between M values (blue) = ',num2str(mean(abs(Cx(:)-Cy3(:)),'all'))]})
xlabel('Generated Correlations. Pink: C = tanh(J) ')
ylabel('Correlation for Indie Pair and given J')

savefig(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.fig'],'compact')
saveas(h,[time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
%saveas(h,[name,'\',time(1:5),'SanityGraphs_N',num2str(N),'_T1E',num2str(log10(T)),'_',time(6:12),'.png'])
