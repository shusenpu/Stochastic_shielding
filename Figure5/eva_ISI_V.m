function [num_ISI,mean_ISI1,std_ISI1,var_ISI1,cv_ISI1] = eva_ISI_V(S,t)
% S=S(:);
% t=t(:);
dt=t(2)-t(1); % msec

%the threshold is subject to change ---by Pu Oct/04/2018
%vthresh=40; % mV. 
    vthresh=-10; % mV.

    spikeidxLE1=find((S(1:(end-1))<=vthresh).*(S(2:end)>vthresh));

    %the right way should be: (linear interpolation of spike time)
    stimeLE=t(spikeidxLE1)+(dt*(vthresh-S(spikeidxLE1))./(S(spikeidxLE1+1)-S(spikeidxLE1)));
    %ISI of the spike
    ISI=diff(stimeLE);
    
    %%% PT revision 8/28/18: discard first 10 spikes
    %ISI=ISI(:);
    ISI=ISI(11:end);
    %%% PT revision 8/28/18: discard first 10 spikes
    
    mean_ISI1=mean(ISI);
    std_ISI1=std(ISI);
    var_ISI1=var(ISI);
    cv_ISI1=sqrt(var_ISI1)/mean_ISI1;
    num_ISI=length(ISI);

end
