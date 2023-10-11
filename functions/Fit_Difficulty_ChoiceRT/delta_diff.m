function p_dim1_is_easiest = delta_diff(kappa,ev1,t1,ev2,t2,uni_coh)
% calculate the probability that dimension 1 is easier than dimension 2, 
% given the tuple (t,e) for each dimension.
% kappa: signal-to-noise from the DTB model
% evx: accumulated evidence for dimension x
% tx: time spent accumulating evidence evx
% uni_coh: unique values of signed coherence. If left empty, it creates the
% usual ones

if nargin<6 || isempty(uni_coh)
    uni_coh = [0.0000001,0.512.*0.5.^[0:4]];
    uni_coh = sort([-uni_coh,uni_coh]);
end

prior_coh = ones(size(uni_coh));
prior_coh = prior_coh/sum(prior_coh);

if t1==0
    p_coh1_e = prior_coh;
else
    p_e_coh1 = normpdf(ev1,kappa*uni_coh*t1,sqrt(t1));
    p_coh1_e = p_e_coh1.*prior_coh;
    p_coh1_e = p_coh1_e/sum(p_coh1_e);
end

if t2==0
    p_coh2_e = prior_coh;
else
    p_e_coh2 = normpdf(ev2,kappa*uni_coh*t2,sqrt(t2));
    p_coh2_e = p_e_coh2.*prior_coh;
    p_coh2_e = p_coh2_e/sum(p_coh2_e);
end

p = p_coh2_e'*p_coh1_e;
[X,Y] = meshgrid(uni_coh,uni_coh);
w = double(abs(X)>abs(Y)) + 1/2*(abs(X)==abs(Y));

p_dim1_is_easiest = sum(sum(p.*w));

end

