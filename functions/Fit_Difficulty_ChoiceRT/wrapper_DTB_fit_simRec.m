function [err,m,fig_handle] = wrapper_DTB_fit_simRec(theta,coh_pair,choices,RT,model_flag,nTrialsCond)

kappa   = theta(1);
ndt_mu  = theta(2);
B0      = theta(3);
a       = theta(4);
d       = theta(5);

if model_flag==3 || model_flag==5
    B0mini = theta(6);
else
    B0mini = nan;
    theta(6) = nan;
end



[m.choice, m.decision_time, m.coh_model, m.idx_coh_model, m.uni_coh] = ...
    dtb_abs_dv_simRec(kappa, B0, a, d, coh_pair, B0mini, model_flag,nTrialsCond);


%%

optim_flag = 2;
switch optim_flag 
    case 1 % RT only

        p_RT = nan(size(RT));        
        n = size(m.uni_coh,1);
        filt = RT>ndt_mu;
        for i=1:n

            I = all(coh_pair==m.uni_coh(i,:),2) & filt;
            J = m.idx_coh_model == i;
            pd = fitdist(m.decision_time(J),'kernel','Kernel','epanechnikov','support','positive');
            rt = RT(I)-ndt_mu;
            p_RT(I) = pd.pdf(rt);

        end

        pPred = p_RT;
        pPred(pPred<=0 | isnan(pPred)) = eps;
        err = -sum(log(pPred));
        
    case 2
        
        p_RT_choice = nan(size(RT,1),1);        
        n = size(m.uni_coh,1);
        filt = RT>ndt_mu;
        for i=1:n
            for k=1:2
                K = all(coh_pair==m.uni_coh(i,:),2) & filt;
                I = K & choices==(k-1);
                J = m.idx_coh_model == i & m.choice==(k-1);
                if sum(J)>1 && sum(I)>0 && sum(~isnan(m.decision_time(J)))>1
                    pd = fitdist(m.decision_time(J),'kernel','Kernel','epanechnikov','support','positive');
                    rt = RT(I)-ndt_mu;
                    p_RT_choice(I) = pd.pdf(rt) * nanmean(choices(K)==(k-1));
                end
                
            end
        end

        pPred = p_RT_choice;
        pPred(pPred<=0 | isnan(pPred)) = eps;
        err = -sum(log(pPred));
        
end

% non-dec time
% G = normpdf(m.t,ndt_mu,ndt_sigma);
% ind = find(G>0,1,'last');
% G = G(1:ind)/sum(G);


%%

% err_method = 1;
% 
% switch err_method
%     case 1 % what I've been doing, mainly
%         p_resp_times_1 = conv2(dec_times(:,:,1),G);
%         p_resp_times_2 = conv2(dec_times(:,:,2),G);
% 
%         aRT = clip(RT,0,max(m.t));
%         response_step = round(aRT/dt);
% 
%         IND = sub2ind(size(p_resp_times_1),[1:ntr]',response_step);
%         pPred = p_resp_times_1(IND).*(choices==0) + p_resp_times_2(IND).*(choices==1);
%         pPred_aux = 0.5 * (p_resp_times_1(IND)+p_resp_times_2(IND));
% 
%         pPred(pPred==0) = pPred_aux(pPred==0)/sims_per_trial; %might be too low
%         % last resource:
%         pPred(pPred==0) = eps;
%         
%         err = -sum(log(pPred));
%         
%             
% end

% %%
%plot_flag = 1;
m.response_time = m.decision_time + ndt_mu;
%fig_handle = do_plot(plot_flag, choices, coh_pair, RT, m,1);

%% print
var_names = {'kappa', 'ndt_mu', 'B0', 'a', 'd','B0mini'};
fprintf_params(var_names,err,theta);


end

