function run_prep_decision_space_plot(datadir)

if nargin==0
    datadir = './unknown_color/';
end

d = importdata(fullfile(datadir,'output.txt'));
s = load(fullfile(datadir,'simulatedTrials.mat'));
run(fullfile(datadir,'params.m')); 
x = -pars.MAXDV:pars.DVDELTA:pars.MAXDV;


%%


u = 2:4:16;
ut = 2 * u * pars.deltatime;
N = round(length(u));

for i=1:N

    I = d.data(:,2)==(u(i)) & d.data(:,4)==(u(i)) & d.data(:,5)==0;
    D = [d.data(I,1), d.data(I,3), d.data(I,7)];
    n = sqrt(size(D,1));
    
    
    y = reshape(D(:,3),[n,n]);
    y = ismember(y,[8,9]); % listen
    
    

    dat(i).x = x;
    dat(i).y = y;
    dat(i).time = ut(i);
    
end

save(fullfile(datadir,'for_decision_space_plot'),'dat');



end
