
load('fig5data.mat')

%% K gates
% yneg=log(mean_K)-log(mean_K-1.96*std_K./sqrt(400));
% ypos=log(mean_K+1.96*std_K./sqrt(400))-log(mean_K);
% x=1:8;
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


%% Na gates
% yneg1=log(mean_Na)-log(mean_Na-1.96*std_Na./sqrt(400));
% ypos1=log(mean_Na+1.96*std_Na./sqrt(400))-log(mean_Na);
% x1=1:20;
% log_err95_K=log(1.96*std_K./sqrt(400));
figure
% bar(log(mean_K))
%bar(log(mean_K),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5)
bar(log(mean_Na),'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor',[0.3010 0.7450 0.9330],'LineWidth',1.5)

set(gca,'Fontsize',16)
hold on
errorbar(x1,log(mean_Na),yneg1,ypos1,'ms','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','magenta','MarkerFaceColor','white')
xlabel('Edge Number (Na^+)')
ylabel('$\log(\sigma^2$(ISI))','Interpreter','latex')
grid on