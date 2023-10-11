function [err,m] = wrapper_DTB_fit(theta,coh_pair,choices,RT,model_flag)

kappa   = theta(1);
ndt_mu  = theta(2);
B0      = theta(3);
a       = theta(4);
d       = theta(5);

if model_flag==3
    B0mini = theta(6);
else
    B0mini = nan;
    theta(6) = nan;
end


%%
[m.choice, m.decision_time, m.coh_model, m.idx_coh_model, m.uni_coh] = ...
    dtb_abs_dv(kappa, B0, a, d, coh_pair, B0mini, model_flag);


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

m.response_time = m.decision_time + ndt_mu;

%% print
var_names = {'kappa', 'ndt_mu', 'B0', 'a', 'd','B0mini'};
fprintf_params(var_names,err,theta);


end

