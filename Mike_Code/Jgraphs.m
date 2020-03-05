%% JGraphs.m
% For comparing inferred J values


function JGraphs(JHstruct,sanity,N,T)%N,jn,h_on,sparsity,time,T,name)

% First for Js
Jy = [JHstruct.Jsparse];
Jx1 = [sanity.mfJ]; %mean field
Jx2 = [sanity.tapJ]; %tap

figure
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

