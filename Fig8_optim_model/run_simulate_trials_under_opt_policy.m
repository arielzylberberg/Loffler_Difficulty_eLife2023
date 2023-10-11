function run_simulate_trials_under_opt_policy(datafolder)
% e.g., run_simulate_trials_under_opt_policy('./unknown_color');

if nargin==0
    datafolder = pwd;
end

%%

d = importdata(fullfile(datafolder,'output.txt'));
for i = 1:length(d.textdata)
    eval([d.textdata{i},'=d.data(1:end,i);']);
end

%%

run(fullfile(datafolder,'params.m'));


%%

[RT,coh,actionSelected,dim_queried, DecTime, DV] = fn_simulate_trials_under_policy(bestAction,pars);


%%

task_choice = ismember(actionSelected,[0:3]);
motion_choice(:,1) = ismember(actionSelected,[0,1,4,6]);
motion_choice(:,2) = ismember(actionSelected,[0,2,4,5]);

easier_task = abs(coh(:,1))>abs(coh(:,2));
I = coh(:,1)==coh(:,2);
easier_task(I) = rand(sum(I),1)>0.5;

% task_choice = ismember(actionSelected,[0,1]);
% motion_choice = zeros(size(easier_task));
% motion_choice(task_choice==1 & actionSelected==0) = 1; % rightward
% motion_choice(task_choice==0 & actionSelected==2) = 1; % rightward
I = task_choice==1;
coh_chosen = nan(size(easier_task));
coh_unchosen = nan(size(easier_task));
coh_chosen(I) = coh(I,1);
coh_unchosen(I) = coh(I,2);
coh_chosen(I==0) = coh(I==0,2);
coh_unchosen(I==0) = coh(I==0,1);

c_task = task_choice==easier_task;
c_motion = (coh_chosen>0 & motion_choice==1) | (coh_chosen<0 & motion_choice==0);




% num switches
d = dim_queried;
I = ismember(dim_queried,[8,9]);
d(~I) = nan;

x = abs(diff(d,[],2));
num_switches = nansum(x,2);
t_first_switch = nan(size(dim_queried,1),1);
for i=1:length(t_first_switch)
    ff = find(x(i,:)==1, 1);
    if isempty(ff)
        t_first_switch(i) = nan;
    else
        t_first_switch(i) = ff;
    end
end

ntr = length(actionSelected);
max_switches = 15;
t_switches = nan(ntr,max_switches);
for i=1:ntr
    cont = 0;
    for j=2:size(dim_queried,2)
        if (dim_queried(i,j-1)~=dim_queried(i,j)) && ~isnan(dim_queried(i,j))
            %         if (dim_queried(i,j-1)==4 && dim_queried(i,j)==5) || ...
            %             (dim_queried(i,j-1)==5 && dim_queried(i,j)==4)
            cont = cont+1;
            t_switches(i,cont) = j-1; % switched AFTER this step
        end
    end
end



savefilename = fullfile(datafolder,'simulatedTrials');

save(savefilename,'RT','coh','actionSelected','dim_queried','DV',...
    'DecTime','c_task','c_motion','t_first_switch','num_switches','t_switches',...
    'coh_chosen','coh_unchosen','task_choice','motion_choice');



end


