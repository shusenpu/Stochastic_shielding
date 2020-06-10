function [ISI] = eva_ISI(V,t)
dt=t(2)-t(1);
vthresh=-10;
t=t(:);
V=V(:);
spikeidx=find((V(1:(end-1))<=vthresh).*(V(2:end)>vthresh));

%discard the first few spikes
spikeidx=spikeidx(11:end-1); %ensure a full cycle at the end
% min_idx=nan(length(spikeidx)-1,1);
% for i=1:length(spikeidx)-1
%     
%     [~,n]=min(V(spikeidx(i):spikeidx(i+1)));
%     
%     min_idx(i)=spikeidx(i)+n;   
% end

% figure, plot(t,V_MC)
% hold on
% plot(t(min_idx),V_MC(min_idx),'b+')
% plot(t(spikeidx_MC),V_MC(spikeidx_MC),'r+')

spiketime=t(spikeidx)+(dt*(vthresh-V(spikeidx))./(V(spikeidx+1)-V(spikeidx)));
%ISI of the spike
ISI=diff(spiketime);

end

