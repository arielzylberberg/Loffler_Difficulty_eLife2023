function [choice, decision_time, coh_model, idx_coh_model, uni_coh] = ...
    dtb_abs_dv(kappa, B0, a, d, coh, B0mini, model_flag)

% coh1 = 0.5;
% coh2 = 0.5;
% kappa = 10;
%
% a = 1;
% d = 0;
% B0 = 1.5;


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

mu1 = kappa * coh1 *dt;
mu2 = kappa * coh2 *dt;
sigma1 = sqrt(dt);
sigma2 = sqrt(dt);

e1 = bsxfun(@times, noise1, sigma1) + mu1;
e2 = bsxfun(@times, noise2, sigma2) + mu2;

USfunc = 'Logistic';

switch model_flag
    case 1 % Race model
        
        dv1 = cumsum(e1,2);
        dv2 = cumsum(e2,2);
        
    case 2 % Difference model [Difference in abs(DV's)]
        
        dv1 = cumsum(e1,2);
        dv2 = cumsum(e2,2);
        dv = abs(dv1) - abs(dv2);
        
    case 3 % Two-step model
        
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
        
        dv = dv1.*w(:,1) - dv2.*w(:,2);
        
        % now for each trial, set the dv to zero up to the end of the mini
        % decision, so that it cannot cross before that
        dt_step = findclose(t, decision_time_mini);
        
        for i=1:ntr
            dv(i,1:dt_step(i)-1) = 0;
        end
        
    case 4 % Absolute momentary evidence model
        dv = cumsum(abs(e1) - abs(e2),2);
        
end


Bup = expand_bounds(t,B0,a,d,USfunc);

if model_flag == 1 % first DV to hit bound wins
    
    % for left stimulus (DV1)
    t_step1_dv1 = single_bound_cross_dynamic(dv1, Bup);
    t_step2_dv1 = single_bound_cross_dynamic(-dv1, Bup);
    [~,isWinner_dv1,decision_time_dv1] = winning_race([t_step1_dv1,t_step2_dv1]*dt);
    
    % for right stimulus (DV2)
    t_step1_dv2 = single_bound_cross_dynamic(dv2, Bup);
    t_step2_dv2 = single_bound_cross_dynamic(-dv2, Bup);
    [~,isWinner_dv2,decision_time_dv2] = winning_race([t_step1_dv2,t_step2_dv2]*dt);
    
    [decision_time, idx] = min([decision_time_dv1 decision_time_dv2],[],2);
    
    choice = -1*idx+2; % 0 if right stimulus is winning, 1 if left stimulus is winning
    
else
    
    Bup = expand_bounds(t,B0,a,d,USfunc);
    t_step1 = single_bound_cross_dynamic(dv, Bup);
    t_step2 = single_bound_cross_dynamic(-dv, Bup);
    [~,isWinner,decision_time] = winning_race([t_step1,t_step2]*dt);
    
    % [winner,isWinner,dec_time,Bup,cev] = dtb_multi(ev,t,a,d,B0,Breflect,USfunc);
    choice = isWinner(:,1);
    
end


% assume seriality (time-sharing)
serial_flag = 1;
if serial_flag
    decision_time = decision_time*2;
end

end
