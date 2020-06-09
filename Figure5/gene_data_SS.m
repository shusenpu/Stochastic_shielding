

% Written by Shusen Pu and Peter Thomas
% June, 2020
%sample script to generate sample SS traces

% to make the plot in Figure 5 at least 1000 of this script is needed
% we used parall computing to manage this task

% generate 1000 different trajectories


%%

%simulation time
t=0:0.008:7000;
t=t(:);

%parameters generating spikes
Ifunc=@(t)10; 
Area=100;
Nt=length(t); 

%use the full strength of nouise
eps=1;

%number of simulations
Nsim=16;

%save data
Na_SS=nan(20,Nsim);
K_SS=nan(8,Nsim);
seed=nan(Nsim,1);

parfor i=1:Nsim
    seed(i)=sum(clock)*10^6+i;
    rng(seed(i))
    
    KNoise1=randn(8,Nt+1);
    NaNoise1=randn(20,Nt+1);
    K_SS_temp=nan(8,1);
    Na_SS_temp=nan(20,1);

    %+++++++++ SS for K+++++++++++++++++
    Na_mask=zeros(20,1);
    for j=1:8
        %include the jth edge of K
        K_mask=zeros(8,1);
        K_mask(j)=eps;
        K1_temp=EM_V(t, Ifunc, Area, NaNoise1, KNoise1,Na_mask,K_mask);
        [num_ISI,mean_ISI1,std_ISI1,var_ISI1,cv_ISI1] = eva_ISI_V(K1_temp,t);
        K_SS_temp(j)=var_ISI1;      
    end
    
        %+++++++++ SS for Na+++++++++++++++++
    K_mask=zeros(8,1);
    for k=1:20
        %include the jth edge of Na
        Na_mask=zeros(20,1);      
        Na_mask(k)=eps;
        K1_temp=EM_V(t, Ifunc, Area, NaNoise1, KNoise1,Na_mask,K_mask);
        [num_ISI,mean_ISI1,std_ISI1,var_ISI1,cv_ISI1] = eva_ISI_V(K1_temp,t);
        Na_SS_temp(k)=var_ISI1;      
    end
    
    Na_SS(:,i)=Na_SS_temp;
    K_SS(:,i)=K_SS_temp;
   
end

%%
save Na_SS Na_SS
save K_SS K_SS
