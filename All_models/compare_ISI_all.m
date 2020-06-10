
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Shusen Pu and Peter Thomas
% June, 2020
% pathwise equivalency including either Na and K gates
% sample code to to calculate the L1-Wasserstein distance

n_MC=length(sort_MC);

S_14D=sort(ISI_New14D_V);
n_14D=length(S_14D);


%%

% pick_14D = sort(randperm(n_14D,n_MC));
% L_14D=mean(abs(S_14D(pick_14D)-sort_MC));
% 
% pick_Dang = sort(randperm(n_Dang,n_MC));
% L_Dang=mean(abs(S_Dang(pick_Dang)-sort_MC));
% 
% 
% pick_Fox94 = sort(randperm(n_Fox94,n_MC));
% L_Fox94=mean(abs(S_Fox94(pick_Fox94)-sort_MC));
% 
% 
% pick_Fox97 = sort(randperm(n_Fox97,n_MC));
% L_Fox97=mean(abs(S_Fox97(pick_Fox97)-sort_MC));
% 
% pick_Orio = sort(randperm(n_Orio,n_MC));
% L_Orio=mean(abs(S_Orio(pick_Orio)-sort_MC));


% S_Dang=sort(ISI_D_All);

%% repeatly draw samples for the difference
tic
Nsimj=100;
L_14D=nan(Nsimj,1);

for simj=1:Nsimj

pick_14D = sort(randperm(n_MC,n_14D));
L_14D(simj)=mean(abs(S_14D-sort_MC(pick_14D)));

end
trun=toc

%% align all time data

%% view time data

figure
hist(T_New14D_V,500)
[ST_Dang,indx]=sort(T_Dang);
t1ff=diff(ST_Dang);
figure
plot(t1ff)
tt1=find(t1ff>20);

ind_jump=tt1(2:3);
ind_major=indx((ind_jump(1)+1):ind_jump(2));

%%
% T_Dang=T_Dang(ind_major);
% T_14D=T_14D(ind_major);
% T_Fox94=T_Fox94(ind_major);
% T_Fox97=T_Fox97(ind_major);
% T_Orio=T_Orio(ind_major);

%%
% std of errors
std_Q_D=std(L_Dang);
std_Q_14D=std(L_14D);
std_Q_Fox94=std(L_Fox94);
std_Q_Fox97=std(L_Fox97);
std_Q_Orio=std(L_Orio);


std_T_D=std(T_Dang);
std_T_14D=std(T_14D);
std_T_Fox94=std(T_Fox94);
std_T_Fox97=std(T_Fox97);
std_T_Orio=std(T_Orio);


x_Fox94=(mean(T_Fox94));
y_Fox94=(mean(L_Fox97));
yp_Fox94=(1.96*std_Q_Fox94/sqrt(length(L_Fox97)));
xp_Fox94=(1.96*std_T_Fox94/sqrt(length(L_Fox97)));

x_Fox97=(mean(T_Fox97));
y_Fox97=(mean(L_Fox97));
yp_Fox97=(1.96*std_Q_Fox97/sqrt(length(L_Fox97)));
xp_Fox97=(1.96*std_T_Fox97/sqrt(length(L_Fox97)));

x_D=(mean(T_Dang));
y_D=(mean(L_Dang));
yp_D=(1.96*std_Q_D/sqrt(length(L_Dang)));
xp_D=(1.96*std_T_D/sqrt(length(L_Dang)));

x_14D=(mean(T_14D));
y_14D=(mean(L_14D));
yp_14D=(1.96*std_Q_14D/sqrt(length(L_14D)));
xp_14D=(1.96*std_T_14D/sqrt(length(L_14D)));

x_Orio=(mean(T_Orio));
y_Orio=(mean(L_Orio));
yp_Orio=(1.96*std_Q_Orio/sqrt(length(L_Orio)));
xp_Orio=(1.96*std_T_Orio/sqrt(length(L_Orio)));


%%
figure
% errorbar(x_IG,y_IG,yp_IG,yp_IG,xp_IG,xp_IG,'r-+','MarkerSize',12,'MarkerEdgeColor','red','MarkerFaceColor','red')

errorbar(x_Fox94,y_Fox94,yp_Fox94,yp_Fox94,xp_Fox94,xp_Fox94,'g-o','MarkerSize',12,'MarkerEdgeColor','green','MarkerFaceColor','green')
hold on
errorbar(x_Fox97,y_Fox97,yp_Fox97,yp_Fox97,xp_Fox97,xp_Fox97,'k-d','MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','black')

errorbar(x_Orio,y_Orio,yp_Orio,yp_Orio,xp_Orio,xp_Orio,'c-s','MarkerSize',12,'MarkerEdgeColor','cyan','MarkerFaceColor','cyan')

errorbar(x_14D,y_14D,yp_14D,yp_14D,xp_14D,xp_14D,'b-*','MarkerSize',12,'MarkerEdgeColor','blue','MarkerFaceColor','blue')

%errorbar(x_SS,y_SS,yp_SS,yp_SS,xp_SS,xp_SS,'m-v','MarkerSize',12,'MarkerEdgeColor','magenta','MarkerFaceColor','magenta')

errorbar(x_D,y_D,yp_D,yp_D,xp_D,xp_D,'r-^','MarkerSize',12,'MarkerEdgeColor','red','MarkerFaceColor','red')

% errorbar(x,y,yp,yp,xp,xp,'-o','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor','red')

set(gca,'Fontsize',16)
grid on
legend('Fox94','Fox97','Orio','14D','Dangerfield')
title('Now')


%%
