%%% diffchk.m

% mean checks on differences

function diffchkstruct = diffchk(jn,N,T,sparsity,time,sanityinput,sanitydimer,sanitydiscon)

%---------------------------------------------------------------------------------------------- old
% sanity = struct('sanity',{},'sanitydimer',{},'sanitydiscon',{}) %,'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{});
% sanity = struct('sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{},'sanity',{},'sanitydimer',{},'sanitydiscon',{});

% for i = 1:jn

%% Mean Field checks
% ass info for what should equal what under roudi's suggestions
% Connectivity matrix
% sanity(i).sanity.dmfC = mean([sanityinput.dmfC],'all');
% sanity(i).sanitydimer.dmfC = mean([sanitydimer.dmfC],'all');
% sanity(i).sanitydiscon.dmfC = mean([sanitydiscon.dmfC],'all');
% 
% % MF J  matrix
% sanity(i).sanity.dmfJ = mean([sanityinput.dmfJ],'all');
% sanity(i).sanitydimer.dmfJ = mean([sanitydimer.dmfJ],'all');
% sanity(i).sanitydiscon.dmfJ = mean([sanitydiscon.dmfJ],'all');
% 
% % MF h values
% sanity(i).sanity.dmfh = mean([sanityinput.dmfh],'all');
% sanity(i).sanitydimer.dmfh = mean([sanitydimer.dmfh],'all');
% sanity(i).sanitydiscon.dmfh = mean([sanitydiscon.dmfh],'all');
% 
% %% TAP checks
% %TAP Connectivity Matrix
% sanity(i).sanity.dtapC = mean([sanityinput.dtapC],'all');
% sanity(i).sanitydimer.dtapC = mean([sanitydimer.dtapC],'all');
% sanity(i).sanitydiscon.dtapC = mean([sanitydiscon.dtapC],'all');
% 
% % TAP J Matrix
% sanity(i).sanity.dtapJ = mean([sanityinput.dtapJ],'all');
% sanity(i).sanitydimer.dtapJ = mean([sanitydimer.dtapJ],'all');
% sanity(i).sanitydiscon.dtapJ = mean([sanitydiscon.dtapJ],'all');
% 
% % TAP h values
% sanity(i).sanity.dtaph = mean([sanityinput.dtaph],'all');
% sanity(i).sanitydimer.dtaph = mean([sanitydimer.dtaph],'all');
% sanity(i).sanitydiscon.dtaph = mean([sanitydiscon.dtaph],'all');
%---------------------------------------------------------------------------------------------- old

% %% Mean Field checks
% % ass info for what should equal what under roudi's suggestions
% % Connectivity matrix
% sanity.sanity.dmfC = 		mean([sanityinput.dmfC],'all');
% sanity.sanitydimer.dmfC = 	mean([sanitydimer.dmfC],'all');
% sanity.sanitydiscon.dmfC = 	mean([sanitydiscon.dmfC],'all');
% % MF Jatrix
% sanity.sanity.dmfJ = 		mean([sanityinput.dmfJ],'all');
% sanity.sanitydimer.dmfJ = 	mean([sanitydimer.dmfJ],'all');
% sanity.sanitydiscon.dmfJ = 	mean([sanitydiscon.dmfJ],'all');
% % MF hlues
% sanity.sanity.dmfh = 		mean([sanityinput.dmfh],'all');
% sanity.sanitydimer.dmfh = 	mean([sanitydimer.dmfh],'all');
% sanity.sanitydiscon.dmfh = 	mean([sanitydiscon.dmfh],'all');
% %% TAPecks
% %TAP Cectivity Matrix
% sanity.sanity.dtapC = 		mean([sanityinput.dtapC],'all');
% sanity.sanitydimer.dtapC = 	mean([sanitydimer.dtapC],'all');
% sanity.sanitydiscon.dtapC = mean([sanitydiscon.dtapC],'all');
% % TAP atrix
% sanity.sanity.dtapJ = 		mean([sanityinput.dtapJ],'all');
% sanity.sanitydimer.dtapJ = 	mean([sanitydimer.dtapJ],'all');
% sanity.sanitydiscon.dtapJ = mean([sanitydiscon.dtapJ],'all');
% % TAP alues
% sanity.sanity.dtaph = 		mean([sanityinput.dtaph],'all');
% sanity.sanitydimer.dtaph = 	mean([sanitydimer.dtaph],'all');
% sanity.sanitydiscon.dtaph = mean([sanitydiscon.dtaph],'all');


%----------------------------------------------------------------------------------------------
a = [sanityinput.dmfC];
b = [sanitydimer.dmfC];
c = [sanitydiscon.dmfC];
d = [sanityinput.dmfJ];
e = [sanitydimer.dmfJ];
f = [sanitydiscon.dmfJ];
g = [sanityinput.dmfh];
h = [sanitydimer.dmfh];
k = [sanitydiscon.dmfh];
l = [sanityinput.dtapC];
m = [sanitydimer.dtapC];
n = [sanitydiscon.dtapC];
o = [sanityinput.dtapJ];
p = [sanitydimer.dtapJ];
q = [sanitydiscon.dtapJ];
r = [sanityinput.dtaph];
s = [sanitydimer.dtaph];
t = [sanitydiscon.dtaph];

a(~isfinite(a)) = 0;
b(~isfinite(b)) = 0;
c(~isfinite(c)) = 0;
d(~isfinite(d)) = 0;
e(~isfinite(e)) = 0;
f(~isfinite(f)) = 0;
g(~isfinite(g)) = 0;
h(~isfinite(h)) = 0;
k(~isfinite(k)) = 0;
l(~isfinite(l)) = 0;
m(~isfinite(m)) = 0;
n(~isfinite(n)) = 0;
o(~isfinite(o)) = 0;
p(~isfinite(p)) = 0;
q(~isfinite(q)) = 0;
r(~isfinite(r)) = 0;
s(~isfinite(s)) = 0;
t(~isfinite(t)) = 0;




%% Mean Field checks
% ass info for what should equal what under roudi's suggestions
% Connectivity matrix
diffchk.sanity.dmfC = 			mean(a,'all');
diffchk.sanitydimer.dmfC = 		mean(b,'all');
diffchk.sanitydiscon.dmfC = 	mean(c,'all');
diffchk.sanity.dmfJ = 			mean(d,'all');
diffchk.sanitydimer.dmfJ = 		mean(e,'all');
diffchk.sanitydiscon.dmfJ = 	mean(f,'all');
diffchk.sanity.dmfh = 			mean(g,'all');
diffchk.sanitydimer.dmfh = 		mean(h,'all');
diffchk.sanitydiscon.dmfh = 	mean(k,'all');
diffchk.sanity.dtapC = 			mean(l,'all');
diffchk.sanitydimer.dtapC = 	mean(m,'all');
diffchk.sanitydiscon.dtapC =	mean(n,'all');
diffchk.sanity.dtapJ = 			mean(o,'all');
diffchk.sanitydimer.dtapJ = 	mean(p,'all');
diffchk.sanitydiscon.dtapJ =	mean(q,'all');
diffchk.sanity.dtaph = 			mean(r,'all');
diffchk.sanitydimer.dtaph =		mean(s,'all');
diffchk.sanitydiscon.dtaph =	mean(t,'all');
	
	
	
	
	

save([time(1:5),'DiffChkr_N',num2str(N),'_T',num2str(T),'_trials',num2str(jn),'_sprs',num2str(100*sparsity),'_',time(6:12),'.mat'],'diffchk');

diffchkstruct = diffchk;
end
