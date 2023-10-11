function [RT,coh,actionSelected,dim_queried, DecTime, DV] = fn_simulate_trials_under_policy(bestAction,pars)


MAXDV       = pars.MAXDV;
DVDELTA     = pars.DVDELTA;
EV_STEP     = pars.EV_STEP;
NSTEPS      = pars.NSTEPS;
SIGMA       = pars.SIGMA;
% COH         = pars.COH;


%%

DV_VALUES = -MAXDV:DVDELTA:MAXDV;
NDV = length(DV_VALUES);

nTrials = 200000;
RT      = nan(nTrials,1);

evStep1 = repmat(EV_STEP,ceil(nTrials/length(EV_STEP)),1);
evStep1 = evStep1(:);
evStep1 = shuffle(evStep1);
evStep1 = evStep1(1:nTrials);
evStep2 = shuffle(evStep1);

coh = [evStep1,evStep2];
DecTime = [];
dim_queried = nan(nTrials,2*NSTEPS);

actionSelected = nan(nTrials,1);
for iTrial = 1:nTrials
    disp(['------'])
    disp(['trial: ',num2str(iTrial)])

    coh1 = evStep1(iTrial);
    coh2 = evStep2(iTrial);
    
    disp(['coh1: ',num2str(coh1)])
    disp(['coh2: ',num2str(coh2)])
    
    endtrial = 0;
    dv = [];
    
    iStep1 = 0; %?????
    iStep2 = 0;
    
    ba = 8+(rand>0.5); % assume first listen at random
%     ba = 8; %assume I start listening to task1
    allEv1 = 0;
    allEv2 = 0;
    
    dim_queried(iTrial,iStep1+iStep2+1) = ba;
    
    while not(endtrial)
        
        dim_queried(iTrial,iStep1+iStep2+1) = ba;
        
        if ba==8 % listen task1
            
            iStep1 = iStep1 + 1;
            currentEv1  = randn*SIGMA + coh1;
            allEv1      = [allEv1;currentEv1];
            disp(['listen 1: ', num2str(currentEv1)])
            
        elseif ba==9 % listen task2
            iStep2 = iStep2 + 1;
            currentEv2  = randn*SIGMA + coh2;
            allEv2      = [allEv2;currentEv2];
            disp(['listen 2: ', num2str(currentEv2)])
            
        else
            
            actionSelected(iTrial) = ba;
            endtrial = 1;
            RT(iTrial) = iStep1 + iStep2;
            disp(['saccade: ',num2str(ba)])
        
        end
        
        if ~endtrial
            currentDv1 = DV_VALUES(findclose(DV_VALUES,sum(allEv1)));
%             iDV1 = find(currentDv1==DV_VALUES) - 1;
            iDV1 = find(currentDv1==DV_VALUES);
            
            currentDv2 = DV_VALUES(findclose(DV_VALUES,sum(allEv2)));
%             iDV2 = find(currentDv2==DV_VALUES) - 1;
            iDV2 = find(currentDv2==DV_VALUES);
            
            iQueried = double(ba==9);
%             index = (iDV2-1)*(NSTEPS) + (iStep2+1) + NSTEPS*NDV*NSTEPS*(iDV1-1) + NDV*NSTEPS*iStep1 + 1;
            index = (iDV2-1)*(NSTEPS) + (iStep2) + NSTEPS*NDV*NSTEPS*(iDV1-1) + NDV*NSTEPS*iStep1 + 1 + ...
                iQueried*NSTEPS^2*NDV^2;
            
            % index = find(d.data(:,1)==currentDv1 & d.data(:,2)==iStep1 & d.data(:,3)==currentDv2 & d.data(:,4)==iStep2 & d.data(:,5)==iQueried);
            try
                ba = bestAction(index);
            catch
                aa
            end
            
            dv = [dv; currentDv1, currentDv2];
        
        end
        
        
%         try
%             ba = bestAction(inds);
%         catch
%             aa
%         end
%         if ba<=3
%             actionSelected(iTrial) = ba;
%             endtrial = 1;
%             RT(iTrial) = currentStep1+1;
%             disp('saccade')
%         end
    end
    DV{iTrial} = dv;
    DecTime = [DecTime; iStep1, iStep2];
end




if isfield(pars,'switchDimTime')
    d = dim_queried;
    I = ismember(dim_queried,[8,9]);
    d(~I) = nan;

    x = abs(diff(d,[],2));
    num_switches = nansum(x,2);
    RT = RT * pars.deltatime + pars.nonDecisionTime + pars.switchDimTime*num_switches;
    
else
    RT = RT * pars.deltatime + pars.nonDecisionTime;
end



DecTime = DecTime * pars.deltatime;

coh = coh/(pars.KAPPA*pars.deltatime);

%%
% c1 = sign(evStep1)==sign(1-2*actionSelected);
% 
% c1 = double(c1);
% c1(isnan(actionSelected))=nan;
% 
% [~,~,ind] = unique(evStep1);
% coh1 = COH(ind);
% 
% doPlot = false;
% if doPlot
%     p = publish_plot(1,2);
%     p.next();
%     [t,x,s] = curva_media(1-actionSelected,coh1,[],0);
%     errorbar(t,x,s);
%     xlabel('coh')
%     ylabel('p-right')
%     
%     p.next();
%     [t,x,s] = curva_media(RT,coh1,[],0);
%     errorbar(t,x,s);
%     xlabel('coh')
%     ylabel('rt (sec)')
%     
%     p.format('FontSize',14)
% end






