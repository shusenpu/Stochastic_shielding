
%simulation time
t=0:0.001:(0.001*100000);

%parameters
Ifunc=@(t)10; 
Area=100;

Nt=length(t); 

% mask will shield selected edges
Na_mask=ones(20,1);
K_mask=zeros(8,1);

% Fox and Lu noise
KNoise1=randn(5,Nt+1);
NaNoise1=randn(8,Nt+1);


%[Y3] = PWC(t, Ifunc, Area, Na_mask,K_mask,NaNoise1,KNoise1);

% from Fox and Lu's method to the 14-by-28D model
[Y0] = PWC_H_exact(t, Ifunc, Area, Na_mask,K_mask,NaNoise1,KNoise1);

%% from 14D HH to Fox and Lu

% 14-by-28 D stochastic HH model noise
KNoise2=0; % in this simulation K gates are not included by Stochastic Shielding method
NaNoise2=randn(20,Nt+1);
[Y] = PWC_exact(t, Ifunc, Area, Na_mask,K_mask,NaNoise2,KNoise2);

%% from Fox and Lu to 14D
figure

subplot(2,3,1)
plot(t,Y0(:,2),'b--','LineWidth',2)
hold on
plot(t,Y0(:,16),'k-','LineWidth',1)
%legend('HH (Na)', 'Fox and Lu (Na)')
%title('From 14D to Fox and Lu')
grid on
set(gca,'Fontsize',16)
ylabel('V')

subplot(2,3,4)
plot(t,Y0(:,2)-Y0(:,16))
grid on
%ylabel('mv')
ylabel('voltage difference')
%xlabel('Time (ms)')
set(gca,'Fontsize',16)
xlabel('Time (ms)')

subplot(2,3,2)
plot(t,Y0(:,10),'b--','LineWidth',2)
hold on
plot(t,Y0(:,24),'k-','LineWidth',1)
%legend('14D HH', 'Fox and Lu')
% title('B: From Fox and Lu to 14D HH')
grid on
set(gca,'Fontsize',16)
ylabel('M_{31}')

subplot(2,3,5)
plot(t,Y0(:,10)-Y0(:,24))
grid on
%ylabel('mv')
ylabel('difference in M_{31}')
xlabel('Time (ms)')
set(gca,'Fontsize',16)

subplot(2,3,3)
plot(t,Y0(:,15),'b--','LineWidth',2)
hold on
plot(t,Y0(:,29),'k-','LineWidth',1)
legend('Fox and Lu', '14D')
%title('From 14D to Fox and Lu')
grid on
set(gca,'Fontsize',16)
ylabel('N_{4}')

subplot(2,3,6)
plot(t,Y0(:,15)-Y0(:,29))
grid on
%ylabel('mv')
ylabel('difference in N_{4}')
%xlabel('Time (ms)')
set(gca,'Fontsize',16)
xlabel('Time (ms)')

%% from 14D HH to Fox and Lu

figure

subplot(2,3,1)
plot(t,Y(:,2),'b--','LineWidth',2)
hold on
plot(t,Y(:,16),'k-','LineWidth',1)
%legend('HH (Na)', 'Fox and Lu (Na)')
%title('From 14D to Fox and Lu')
%title('A: From 14D to Fox and Lu')
grid on
set(gca,'Fontsize',16)
ylabel('V')

subplot(2,3,4)
plot(t,Y(:,2)-Y(:,16))
grid on
%ylabel('mv')
ylabel('voltage difference')
%xlabel('Time (ms)')
set(gca,'Fontsize',16)
xlabel('Time (ms)')
subplot(2,3,2)
plot(t,Y(:,10),'b--','LineWidth',2)
hold on
plot(t,Y(:,24),'k-','LineWidth',1)
%legend('14D HH', 'Fox and Lu')
grid on
set(gca,'Fontsize',16)
ylabel('M_{31}')

subplot(2,3,5)
plot(t,Y(:,10)-Y(:,24))
grid on
%ylabel('mv')
ylabel('difference in M_{31}')
xlabel('Time (ms)')
set(gca,'Fontsize',16)

subplot(2,3,3)
plot(t,Y(:,15),'b--','LineWidth',2)
hold on
plot(t,Y(:,29),'k-','LineWidth',1)
legend('14D', 'Fox and Lu')
%title('From 14D to Fox and Lu')
grid on
set(gca,'Fontsize',16)
ylabel('N_{4}')

subplot(2,3,6)
plot(t,Y(:,15)-Y(:,29))
grid on
%ylabel('mv')
ylabel('difference in N_{4}')
%xlabel('Time (ms)')
set(gca,'Fontsize',16)
xlabel('Time (ms)')


