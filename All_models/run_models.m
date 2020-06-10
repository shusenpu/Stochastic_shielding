
%% illustration case
dt=0.01;
% t=0:dt:1000;
 Nsim=4;
%% real simulation
dt=0.008;
t=0:dt:84000;
%
t=t';
Nnt=length(t);
% format shortG
Ifunc=@(t)10; 
Area=100;
%T_MC=nan(Nsim,1);
T_14D=nan(Nsim,1);
T_Orio=nan(Nsim,1);
T_SS=nan(Nsim,1);
T_Dangerfield=nan(Nsim,1);
T_Fox97=nan(Nsim,1);
T_Fox94=nan(Nsim,1);

seed=nan(Nsim,1);
% load('seed.mat')

%ISI_data_MC=cell(Nsim,1);
ISI_data_Fox94=cell(Nsim,1);
ISI_data_Fox97=cell(Nsim,1);
ISI_data_14D=cell(Nsim,1);
ISI_data_Orio=cell(Nsim,1);
ISI_data_D=cell(Nsim,1);
ISI_data_SS=cell(Nsim,1);

%T_MC_ind=zeros(Nsim,1,'uint64');
T_Fox97_ind=zeros(Nsim,1,'uint64');
T_Fox94_ind=zeros(Nsim,1,'uint64');
T_14D_ind=zeros(Nsim,1,'uint64');
T_Orio_ind=zeros(Nsim,1,'uint64');
T_D_ind=zeros(Nsim,1,'uint64');
T_SS_ind=zeros(Nsim,1,'uint64');

parfor i=1:Nsim
  
    seed(i)=sum(clock)*10^6+i;
    rng(seed(i))
    fprintf('MC - %d \n', i); 
%     rng(seed(i))
 %   [Y_MC,T_MC_run] = MC(t, Ifunc, Area);
  %  T_MC(i)=T_MC_run;
   % V_MC=Y_MC(:,2);

   
    % Fox and Lu 1997
    fprintf('Fox and Lu 1997- %d \n', i);
%     rng(seed(i))
    [Y_Fox97,T_Fox97_run] = Fox97(t, Ifunc, Area);
    T_Fox97(i)=T_Fox97_run;
    V_Fox97=Y_Fox97(:,2);

    % Fox and Lu 1994
    fprintf('Fox and Lu 1994- %d \n', i);
%     rng(seed(i))
    [Y_Fox94,T_fox94_run] = FoxandLu94(t, Ifunc, Area);
    T_Fox94(i)=T_fox94_run;
    V_Fox94=Y_Fox94(:,2);

    % 14D HH
    fprintf('14D HH- %d \n', i);
%     rng(seed(i))
    [V_full,T_14D_run] =HHSS14D(t, Ifunc, Area);
    T_14D(i)=T_14D_run;

    % Orio and Soudry
    fprintf('Orio and Soudry- %d \n', i);
%     rng(seed(i))
    [V_Orio,T_Orio_run] = OrioSoudry(t, Ifunc, Area);
    T_Orio(i)=T_Orio_run;

     % Dangerfield
    fprintf('Dangerfield- %d \n', i);
%     rng(seed(i))
    [V_Dangerfield,T_D_run] = Dangerfield(t, Ifunc, Area);
    T_Dangerfield(i)=T_D_run;
    
    
    fprintf('Stochastic Shielding New- %d \n', i);
%     rng(seed(i))
    [V_SS,T_SS_run] = SS_Na_K(t, Ifunc, Area);
    T_SS(i)=T_SS_run;
    
    
     %calculate all ISI's
    %[ISI_MC] = eva_ISI(V_MC,t);
    [ISI_Fox94] = eva_ISI(V_Fox94,t);
    [ISI_Fox97] = eva_ISI(V_Fox97,t);
    [ISI_14D] = eva_ISI(V_full,t);
    [ISI_Orio] = eva_ISI(V_Orio,t);
    [ISI_Dangerfield] = eva_ISI(V_Dangerfield,t);
    [ISI_SS] = eva_ISI(V_SS,t);
   % ISI_data_MC{i}=ISI_MC;
    ISI_data_Fox94{i}=ISI_Fox94;
    ISI_data_Fox97{i}=ISI_Fox97;
    ISI_data_14D{i}=ISI_14D;
    ISI_data_Orio{i}=ISI_Orio;
    ISI_data_D{i}=ISI_Dangerfield;
    ISI_data_SS{i}=ISI_SS;

%     disp(i)

end

%% save data
%save T_MC T_MC
save T_14D T_14D
save T_Orio T_Orio
save T_SS T_SS
save T_Dangerfield T_Dangerfield
save T_Fox97 T_Fox97
save T_Fox94 T_Fox94

%save ISI_data_MC ISI_data_MC
save ISI_data_Fox94 ISI_data_Fox94
save ISI_data_Fox97 ISI_data_Fox97
save ISI_data_14D ISI_data_14D 
save ISI_data_Orio ISI_data_Orio
save ISI_data_D ISI_data_D 
save ISI_data_SS ISI_data_SS 

save seed seed
