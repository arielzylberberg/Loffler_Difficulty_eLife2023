theta = [];
for i=1:20
    d = load(['./fits/optim_suj',num2str(i),'_model2.mat']);
    theta = [theta; d.theta];
end

mean(theta)