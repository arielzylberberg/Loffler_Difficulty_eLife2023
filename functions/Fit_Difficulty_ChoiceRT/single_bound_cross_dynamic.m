function t_step = single_bound_cross_dynamic(dv, Bup)

a = dv>=Bup';
t_step = sum(cumprod(~a,2),2)+1;
nt = size(dv,2);
I = t_step==(nt+1);
t_step(I) = nan;

