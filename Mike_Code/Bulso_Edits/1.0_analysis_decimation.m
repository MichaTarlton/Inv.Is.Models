%% ANALYSIS

%Nx = 50;
Nx = 100;
h = 0;
J = 1;
Ntrials = 100;
new_exp = true;

if ~new_exp
    allname = dir(['Experiments18_decimation/1st_round/Exp1_Nx',num2str(Nx),'_h',num2str(h),'*']);
    load(['Experiments18_decimation/1st_round/',allname(1).name],'vNp_inactive','Np','vbeta','M')
elseif new_exp && Nx == 50
    allname = dir(['Experiments18_decimation/last_round/New_Exp1_Nx',num2str(Nx),'_h',num2str(h),'*']);
    load(['Experiments18_decimation/last_round/',allname(1).name],'vNp_inactive','Np','vbeta','M')
elseif new_exp && Nx == 100
    allname = dir(['Experiments18_decimation/last_round_N=100/New_Exp1_Nx',num2str(Nx),'_h',num2str(h),'*']);
    load(['Experiments18_decimation/last_round_N=100/',allname(1).name],'vNp_inactive','Np','vbeta','M')
end

%vbeta = [0.01,1,1.5,2];
vbeta = [0.01,1,1.5];

% initialization
Error_trials = zeros(7,numel(M),numel(vbeta),numel(vNp_inactive),Ntrials);
Fp_trials = zeros(7,numel(M),numel(vbeta),numel(vNp_inactive),Ntrials);
Fn_trials = zeros(7,numel(M),numel(vbeta),numel(vNp_inactive),Ntrials);
Entropy_trials = zeros(Nx+1,numel(M),numel(vbeta),numel(vNp_inactive),Ntrials);
Unique_trials = zeros(Nx+1,numel(M),numel(vbeta),numel(vNp_inactive),Ntrials);

% main loop to unpack results
splevels = vNp_inactive/Np;
for isp = 1:numel(splevels)
    sp = splevels(isp);
    
    for ibeta = 1:numel(vbeta)
        beta = vbeta(ibeta);
        
        if ~new_exp
            name = dir(['Experiments18_decimation/1st_round/Exp1_Nx',num2str(Nx),'_h',num2str(h),...
                '_s',num2str(sp),'_J',num2str(J),'_LocBeta',num2str(beta),'trial*']);
        elseif new_exp && Nx == 50
            name = dir(['Experiments18_decimation/last_round/New_Exp1_Nx',num2str(Nx),'_h',num2str(h),...
                '_s',num2str(sp),'_J',num2str(J),'_LocBeta',num2str(beta),'trial*']);
        elseif new_exp && Nx == 100
            name = dir(['Experiments18_decimation/last_round_N=100/New_Exp1_Nx',num2str(Nx),'_h',num2str(h),...
                '_s',num2str(sp),'_J',num2str(J),'_LocBeta',num2str(beta),'trial*']);
        end

        for itrial = 1:Ntrials

            if ~new_exp
                load(['Experiments18_decimation/1st_round/',name(itrial).name],'Err','Fn','Fp',...
                    'ENTROPIA','NPATTERNS')
            elseif new_exp && Nx == 50
                load(['Experiments18_decimation/last_round/',name(itrial).name],'Err','Fn','Fp',...
                    'ENTROPIA','NPATTERNS')
            elseif new_exp && Nx == 100
                load(['Experiments18_decimation/last_round_N=100/',name(itrial).name],'Err','Fn','Fp',...
                    'ENTROPIA','NPATTERNS')
            end
            
            for im = 1:numel(M)
                
                % Errors
                Error_trials(1,im,ibeta,isp,itrial) = Err.BIC(im);
                Error_trials(2,im,ibeta,isp,itrial) = Err.AIC(im);
                Error_trials(3,im,ibeta,isp,itrial) = Err.MDLl(im);
                Error_trials(4,im,ibeta,isp,itrial) = Err.MDLu(im);
                Error_trials(5,im,ibeta,isp,itrial) = Err.MDLent(im);
                Error_trials(6,im,ibeta,isp,itrial) = Err.MDLcount(im);
                Error_trials(7,im,ibeta,isp,itrial) = Err.L1(im);
                
                % False positive
                Fp_trials(1,im,ibeta,isp,itrial) = Fp.BIC(im);
                Fp_trials(2,im,ibeta,isp,itrial) = Fp.AIC(im);
                Fp_trials(3,im,ibeta,isp,itrial) = Fp.MDLl(im);
                Fp_trials(4,im,ibeta,isp,itrial) = Fp.MDLu(im);
                Fp_trials(5,im,ibeta,isp,itrial) = Fp.MDLent(im);
                Fp_trials(6,im,ibeta,isp,itrial) = Fp.MDLcount(im);
                Fp_trials(7,im,ibeta,isp,itrial) = Fp.L1(im);
                
                % False negative
                Fn_trials(1,im,ibeta,isp,itrial) = Fn.BIC(im);
                Fn_trials(2,im,ibeta,isp,itrial) = Fn.AIC(im);
                Fn_trials(3,im,ibeta,isp,itrial) = Fn.MDLl(im);
                Fn_trials(4,im,ibeta,isp,itrial) = Fn.MDLu(im);
                Fn_trials(5,im,ibeta,isp,itrial) = Fn.MDLent(im);
                Fn_trials(6,im,ibeta,isp,itrial) = Fn.MDLcount(im);
                Fn_trials(7,im,ibeta,isp,itrial) = Fn.L1(im);
                
                % Entropy & unique patters
                Entropy_trials(:,im,ibeta,isp,itrial) = ENTROPIA(:,im);
                Unique_trials(:,im,ibeta,isp,itrial) = NPATTERNS(:,im);
            end
        end
    end
end


%% plotting results

Errors = mean(Error_trials,5);
stdErrors = std(Error_trials,[],5);
Fp = mean(Fp_trials,5);
Fn = mean(Fn_trials,5);

methods = {"BIC","AIC","MDLl","MDLu","MDLent","MDLcount","L1"};
i_methods = logical([1,1,0,0,1,0,1]);
i_methodsLargeM = logical([1,1,0,0,1,0,1]);
%i_methodsLargeM = logical([0,0,0,0,1,0,0]);
v_spars = logical([ones(1,5),0]);
%v_spars = logical(ones(1,6));
vbeta2 = vbeta(1:3);

figure;
%hold on;
ibilde = 0;
for im = 1:numel(M)
    for ibeta = 1:numel(vbeta2)
        ibilde = ibilde + 1;
        subplot(numel(M),numel(vbeta2),ibilde);
        
        % remove l1 for large M
        i_methods = logical([1,1,0,0,1,0,1]);
        %i_methods = logical([0,0,0,0,1,0,0]);
        if M(im) > 3000
            i_methods = i_methodsLargeM;
        end
        
        % loop over methods
        for imet = 1:numel(i_methods)
            if i_methods(imet)
            hold on
            errorbar(splevels(v_spars),squeeze(Errors(imet,im,ibeta,v_spars)),squeeze(stdErrors(imet,im,ibeta,v_spars)),'linewidth',1);
            end
        end
        %ylabel('Errors')
        %xlabel('sparsity')
        axis([min(splevels(v_spars))-.1,max(splevels(v_spars))+.1,-0.01,1.05])
        %title(['N=',num2str(Nx),' M=',num2str(M(im)),' Beta=',num2str(vbeta(ibeta))]);
        %legend(methods{i_methods})
        box on;grid on;
        ylim([-0.05,1.05])
        xlim([min(splevels(v_spars))-.1,max(splevels(v_spars))+.1])
        set(gca,'fontsize',16);
        
    end
end
%title('decimation')

%% also false positive and false negative

figure;
ibilde = 0;
for im = 1:numel(M)
    for ibeta = 1:numel(vbeta)
        ibilde = ibilde + 1;
        subplot(numel(M),numel(vbeta),ibilde);
        
        plot(squeeze(Fp(i_methods,im,ibeta,:))')
        ylabel('False Positive')
        xlabel('sparsity')
        
        title(['N=',num2str(Nx),' M=',num2str(M(im)),' Beta=',num2str(vbeta(ibeta))]);
        legend(methods{i_methods})
    end
end


figure;
ibilde = 0;
for im = 1:numel(M)
    for ibeta = 1:numel(vbeta)
        ibilde = ibilde + 1;
        subplot(numel(M),numel(vbeta),ibilde);
        
        plot(squeeze(Fn(i_methods,im,ibeta,:))')
        ylabel('False Negative')
        xlabel('sparsity')
        
        title(['N=',num2str(Nx),' M=',num2str(M(im)),' Beta=',num2str(vbeta(ibeta))]);
        legend(methods{i_methods})
    end
end

%% entropy

figure;
n = 0:Nx;
ibilde = 0;
for im = 1:numel(M)
    for ibeta = 1:numel(vbeta2)
        ibilde = ibilde + 1;
        subplot(numel(M),numel(vbeta2),ibilde);
        for u = 1:numel(vNp_inactive)
            for it = 1:Ntrials
                hold on;
                plot(n,squeeze(sort(Entropy_trials(:,im,ibeta,u,it)))/log2(M(im)));
            end
        end
        hold on;
        plot(n,n/log2(M(im)),'k:');
        plot(n,ones(1,numel(n)),'k--');
        
        %title(['N=',num2str(Nx),' M=',num2str(M(im)),' Beta=',num2str(vbeta2(ibeta))]);
        %ylabel('Entropy')
        %xlabel('active parameters')
        box on; grid on;
        axis([min(n),max(n),0,1.1])
        set(gca,'fontsize',16);
    end
end

%% unique patterns

figure;
n = 0:Nx;
ibilde = 0;
for im = 1:numel(M)
    for ibeta = 1:numel(vbeta2)
        ibilde = ibilde + 1;
        subplot(numel(M),numel(vbeta2),ibilde);
        for u = 1:numel(vNp_inactive)
            for it = 1:Ntrials
                hold on;
                plot(n,squeeze(sort(Unique_trials(:,im,ibeta,u,it)))/M(im));
            end
        end
        hold on;
        plot(n,n/M(im),'k:');
        plot(n,2.^(n-1)/M(im),'k:');
        plot(n,ones(1,numel(n)),'k--');
        
        %title(['N=',num2str(Nx),' M=',num2str(M(im)),' Beta=',num2str(vbeta(ibeta))]);
        %ylabel('Unique patterns')
        %xlabel('active parameters')
        box on; grid on;
        axis([min(n),max(n),0,1.1])
        set(gca,'fontsize',16);
    end
end





