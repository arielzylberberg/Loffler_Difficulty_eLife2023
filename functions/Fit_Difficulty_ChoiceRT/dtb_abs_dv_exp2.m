function [choice, decision_time, coh_model, idx_coh_model, uni_coh] = ...
        dtb_abs_dv_exp2(kappa, B0, a, d, coh, B0mini, model_flag,trcond)

dt = 0.005;
t = 0:dt:5;

ntr_per_coh = 1000;

% this could be computed outside...
[uni_coh] = unique(coh,'rows');
coh_model = repmat(uni_coh, ntr_per_coh,1);
idx_coh_model = repmat([1:size(uni_coh,1)]', ntr_per_coh,1);

ntr = size(coh_model, 1);
ntimes = length(t);
noise1 = randn(ntr,ntimes);
noise2 = randn(ntr,ntimes);

coh1 = coh_model(:,1);
coh2 = coh_model(:,2);

if trcond == 2 % known color (just convert everything to positive drifts for simplicity; i.e., negative values = negative evidence)
    coh1 = abs(coh1);
    coh2 = abs(coh2);
end

mu1 = kappa * coh1 *dt;
mu2 = kappa * coh2 *dt;
sigma1 = sqrt(dt);
sigma2 = sqrt(dt);

e1 = bsxfun(@times, noise1, sigma1) + mu1;
e2 = bsxfun(@times, noise2, sigma2) + mu2;


switch model_flag
    case 2
        
        dv1 = cumsum(e1,2);
        dv2 = cumsum(e2,2);

        if trcond == 1 % unknown color = abs(DV's)
            dv = abs(dv1) - abs(dv2);
        elseif trcond == 2 % known color = raw DV's
            dv = dv1 - dv2;
        end
        
    case 5 % Mike's mini decision model - only 1 DV needs to cross bound
        
        % get the signs
        dv1 = cumsum(e1,2);
        dv2 = cumsum(e2,2);
        
        % limit the mini decision to 2 second max - then sharp collapse
        Bup_mini = expand_bounds(t,B0mini,1000,2,'Logistic')';
        
        % four combinations
        
        tm1 = single_bound_cross_dynamic((dv1>0 & dv2>Bup_mini) | (dv1>Bup_mini & dv2>0), 1);
        tm2 = single_bound_cross_dynamic((dv1>0 & dv2<-Bup_mini) | (dv1>Bup_mini & dv2<0), 1);
        tm3 = single_bound_cross_dynamic((dv1<0 & dv2>Bup_mini) | (dv1<-Bup_mini & dv2>0), 1);
        tm4 = single_bound_cross_dynamic((dv1<0 & dv2<-Bup_mini) | (dv1<-Bup_mini & dv2<0), 1);
        
        signos = [1,1;1,-1;-1,1;-1,-1]; % signs corresponding to the four combinations
        
        [winner_mini,~,decision_time_mini] = winning_race([tm1, tm2, tm3, tm4]*dt);
        
        w = signos(winner_mini,:);

        
        % now for each trial, set the dv to zero up to the end of the mini
        % decision, so that it cannot cross before that
        dt_step = findclose(t, decision_time_mini);
        
        for i=1:ntr
            dv(i,1:dt_step(i)-1) = 0;
        end
        

    case 6 % difference-of-confidence model
        
        dv1 = cumsum(e1,2);
        dv2 = cumsum(e2,2);

%         prior_over_coh = ones(size(nanunique(coh)); % assumes uniform prior
        ucoh = nanunique(coh(:));
        conf_s1 = calc_conf(dv1, dt, ucoh, kappa);
        conf_s2 = calc_conf(dv2, dt, ucoh, kappa);

        if trcond == 1 % unknown color = abs(DV's)
            dv = abs(log(conf_s1./(1-conf_s1))) - abs(log(conf_s2./(1-conf_s2)));
%             dv = abs(dv1) - abs(dv2);
        elseif trcond == 2 % known color = raw DV's
%             dv = dv1 - dv2;
            dv = log(conf_s1./(1-conf_s1)) - log(conf_s2./(1-conf_s2));
        end

    case 7 % another version of the conf model: threshold on difference of probabilities, not log-odds
        
        dv1 = cumsum(e1,2);
        dv2 = cumsum(e2,2);

%         prior_over_coh = ones(size(nanunique(coh)); % assumes uniform prior
        ucoh = nanunique(coh(:));
        conf_s1 = calc_conf(dv1, dt, ucoh, kappa);
        conf_s2 = calc_conf(dv2, dt, ucoh, kappa);

        if trcond == 1 % unknown color = abs(DV's)
            dv = max(conf_s1,1-conf_s1) - max(conf_s2, 1-conf_s2);
%             dv = abs(dv1) - abs(dv2);
        elseif trcond == 2 % known color = raw DV's
%             dv = dv1 - dv2;
            dv = conf_s1 - conf_s2;
        end

end


USfunc = 'Logistic';
Bup = expand_bounds(t,B0,a,d,USfunc);
t_step1 = single_bound_cross_dynamic(dv, Bup);
t_step2 = single_bound_cross_dynamic(-dv, Bup);
[~,isWinner,decision_time] = winning_race([t_step1,t_step2]*dt);


% assume seriality
serial_flag = 1;
if serial_flag
    decision_time = decision_time*2;
end


% [winner,isWinner,dec_time,Bup,cev] = dtb_multi(ev,t,a,d,B0,Breflect,USfunc);
choice = isWinner(:,1);

end


%%%%%%%%%%%%%%%%%%%%%%%%

function pright_is_correct = calc_conf(dv, dt, uni_coh, kappa)


[ntrials,nt] = size(dv);

aux = repmat(dt * [1:nt],ntrials,1);
ts = aux(:);

w = double(uni_coh>0);
w(uni_coh==0) = 0.5; % distribute across both choices
w = w(:)';

ddvv = dv(:);

ncoh = length(uni_coh);
pdist_over_coh = nan(nt*ntrials, ncoh);

sigma = sqrt(ts);

for i=1:ncoh
    mu = kappa * ts * uni_coh(i);
    pdist_over_coh(:,i) = normpdf(ddvv, mu, sigma);
end

pright_is_correct = pdist_over_coh * w(:) ./sum(pdist_over_coh,2);
pright_is_correct = reshape(pright_is_correct, [ntrials, nt]);

% for i=1:nt
%     acc_ev = dv(:,i);
%     sigma = sqrt(dt*i);
% 
%     for j=1:length(uni_coh)
%         mu = kappa * dt * i * uni_coh(j);
%         pdist_over_coh(j,i,:) = normpdf(acc_ev, mu, sigma);
%     end
% end

% pright_is_correct = reshape(w*pdist_over_coh(:,:), [ntrials,nt]) ./ squeeze(nansum(squeeze(pdist_over_coh),1));



end


