%% Graphs.m
% For Correlation and magnetization sanity checks


function Graphs(Sstruct,sanitydimer,sanitydisc,N,T)%N,jn,h_on,sparsity,time,T,name)

mx = [Sstruct.mfinal];
my = [sanitydisc.th];
mtchk=mean([sanitydisc.mtchk]);

figure
subplot(1,2,1)
scatter(mx(:),my(:))
%axis([-1 1 -1 1])

refline(1,0) 
hold on 
title({'Disconnected J Magnetizations: mi',['N = ',num2str(N)],['T = ',num2str(log10(T))],['Avg. difference between M values = ',num2str(mtchk)]})
xlabel('Generated from Monte Carlo')
ylabel('Calcualted from mi = tanh(h)')

%Correlation

Cx = [sanitydimer.Cij];
%Cy = [sanitydimer.Cpair];
Cy1 = [sanitydimer.Cpair1].*1e-4; %seems to normalizee a bit
Cy2 = [sanitydimer.Cpair2];
Cy3 = [sanitydimer.Cpair3];

%figure
subplot(1,2,2)
scatter(Cx(:),Cy1(:),[],'b')
%axis([-1 1 -1 1])

refline(1,0)
hold on 
scatter(Cx(:),Cy2(:),[],'r')
scatter(Cx(:),Cy3(:),[],'m')

title({'Correlations for Dimers',['N = ',num2str(N)],['T = ',num2str(log10(T))],['Avg. difference between M values (blue) = ',num2str(mean(abs(Cx(:)-Cy1(:)),'all'))]})
xlabel('Generated Correlations. Blue: eq8 w/o irrational, Red: eq8 w irrational, Pink: C = tanh(J) ')
ylabel('Correlation for Indie Pair and given J')

