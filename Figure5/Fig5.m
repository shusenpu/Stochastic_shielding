
% Written by Shusen Pu and Peter Thomas
% June, 2020


% this code is to load data from HPC and make plots
% to generate the data, one need to run gege_data_SS.m first
% after saving all data, one can combine data and make the plot
% the following code is to load data from 100 separate files and make plots


%%


% we use 100 here as an example 
Nsim = 100;    % number of .mat-Files    

K_data=nan(8,Nsim*16);
Na_data=nan(20,Nsim*16);

%%
ii=1;

%code to collect all data
for i=1:Nsim
    file=['trial_',num2str(i)];   
    Files=dir(fullfile(file,'*.mat'));
    ll=length(Files);
    
    if isempty(Files)
        sprintf('%s has no data ii= %d ',file,i)
    else  
       
        cd(file);
        load('K_SS');
        load('Na_SS');
        K_data(:,((i-1)*16+1):i*16)=K_SS;
        Na_data(:,((i-1)*16+1):i*16)=Na_SS;
       
    end    
    
        ii=ii+1;
        cd('..');
end
  
 %% mean and std
 mean_Na=mean(Na_data,2);
 mean_K=mean(K_data,2);
 
 std_Na=nan(20,1);
 std_K=nan(8,1);
 for i=1:20
     std_Na(i)=std(Na_data(i,:));
 end
 
 for i=1:8
     std_K(i)=std(K_data(i,:));
 end
 
 %%
 
 % 95% confidence interval
 err95_K=1.96*std_K./sqrt(1600);
 err95_Na=1.96*std_Na./sqrt(400);

 % 99% confidence interval
 err99_K=2.576*std_K./sqrt(400);
 err99_Na=2.576*std_Na./sqrt(400);

 %% 
 %plot for ISI
x=1:20;
figure
bar(mean_Na)
set(gca,'Fontsize',16)
hold on
errorbar(x,mean_Na,err95_Na,'vertical','bs','MarkerSize',2,'MarkerEdgeColor','black','MarkerFaceColor','white')
xlabel('Edge Number (Na^+)')
ylabel('Var(ISI)')
grid on

%%
x=1:8;
figure
% bar(mean_K)
bar(mean_K,'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5)
set(gca,'Fontsize',16)
hold on
errorbar(x,mean_K,err95_K,'vertical','bs','MarkerSize',2,'MarkerEdgeColor','black','MarkerFaceColor','white')
xlabel('Edge Number (K^+)')
ylabel('Var(ISI)')
grid on

%%
yneg=log(mean_K)-log(mean_K-1.96*std_K./sqrt(1600));
ypos=log(mean_K+1.96*std_K./sqrt(1600))-log(mean_K);
x=1:8;
% log_err95_K=log(1.96*std_K./sqrt(400));
figure
% bar(log(mean_K))
%bar(log(mean_K),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5)
bar(log(mean_K),'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor',[0.3010 0.7450 0.9330],'LineWidth',1.5)

set(gca,'Fontsize',16)
hold on
errorbar(x,log(mean_K),yneg,ypos,'bs','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','black','MarkerFaceColor','white')
xlabel('Edge Number (K^+)')
ylabel('$\log(\sigma^2$(ISI))','Interpreter','latex')

%%
yneg=log(mean_K)-log(mean_K-1.96*std_K./sqrt(1600));
ypos=log(mean_K+1.96*std_K./sqrt(1600))-log(mean_K);
x=1:8;
% log_err95_K=log(1.96*std_K./sqrt(400));
figure
% bar(log(mean_K))
%bar(log(mean_K),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5)
bar(log(mean_K),'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor',[0.3010 0.7450 0.9330],'LineWidth',1.5)

set(gca,'Fontsize',16)
hold on
errorbar(x,log(mean_K),yneg,ypos,'ms','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','magenta','MarkerFaceColor','white')
xlabel('Edge Number (K^+)')
ylabel('$\log(\sigma^2$(ISI))','Interpreter','latex')
grid on

%%
yneg1=log(mean_Na)-log(mean_Na-1.96*std_Na./sqrt(1600));
ypos1=log(mean_Na+1.96*std_Na./sqrt(1600))-log(mean_Na);
x=1:20;
% log_err95_K=log(1.96*std_K./sqrt(400));
figure
% bar(log(mean_K))
%bar(log(mean_K),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5)
bar(log(mean_Na),'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor',[0.3010 0.7450 0.9330],'LineWidth',1.5)

set(gca,'Fontsize',16)
hold on
errorbar(x,log(mean_Na),yneg1,ypos1,'ms','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','magenta','MarkerFaceColor','white')
xlabel('Edge Number (Na^+)')
ylabel('$\log(\sigma^2$(ISI))','Interpreter','latex')
grid on

%%
x=1:20;
err95_Na1=1.96*log(std_Na)./sqrt(400);
figure
bar(log(mean_Na))
set(gca,'Fontsize',16)
hold on
errorbar(x,log(mean_Na),err95_Na1,'vertical','bs','MarkerSize',2,'MarkerEdgeColor','black','MarkerFaceColor','white')
xlabel('Edge Number (Na^+)')
ylabel('$\log(\sigma^2$(ISI))','Interpreter','latex')
grid on


%%
figure
bar(mean_K)
set(gca,'Fontsize',16)
% xlabel('Edge Number (K^+)')
xlabel('Edge Number')
ylabel('Var(ISI)')
grid on


 
 
 
 
