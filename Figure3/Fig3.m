
% Written by Shusen Pu and Peter Thomas
% June, 2020
% convergece of 14D HH model


%% code to verify that the first non-zero eigenvector in parallel to dn,dh
%the code also finds the second lagest eigenvalue of both A_Na and A_K
V2=-125:0.01:40;
% [EV_Na2,EV_K2,CoeffK,CoeffNa,CoeffNah,ErrorK,ErrorNa] = HHeigens1(V2);
[EV_Na2,EV_K2,CoeffK,CoeffNa,CoeffNah,ErrorK,ErrorNa,GL] = HHeigens_new(V2);

figure
subplot(2,1,1)
plot(V2,ErrorK)
ylabel('L2 difference for K')
set(gca,'Fontsize',16)
grid on
subplot(2,1,2)
plot(V2,ErrorNa)
ylabel('L2 difference for Na')
set(gca,'Fontsize',16)
grid on

figure
subplot(2,1,1)
plot(V2,EV_Na2)
ylabel('2nd largest eigenvalue of Na')
set(gca,'Fontsize',16)
grid on
subplot(2,1,2)
plot(V2,EV_K2)
ylabel('2nd largest eigenvalue of K')
set(gca,'Fontsize',16)
grid on
%% find the slowest decay value on plots
max1=max(EV_Na2);max2=max(EV_K2);
lambda_decay=max(max1,max2);

%% plot the convergence
% converge_100_HR

% from Shusen: convergence of 14D trajectories to 4D manifold for HH
% sometime in 2019?
dt=0.0001;
t=0:dt:20;
Ifunc=@(t)10; 
Area=100;
Nt=length(t);

%%
% figure(1)
% hold on

%% random initial conditions
for j=1:30
%from Lemma 2 & 3 in Appendix B
V_min=-77;
V_max=10/0.3+50;
V0 = (V_max-V_min).*rand + V_min;

m00=rand;
m01=(1-m00)*rand; % make sure that m01+m00<=1, same below
m02=(1-m00-m01)*rand;
m03=(1-m00-m01-m02)*rand;
m10=(1-m00-m01-m02-m03)*rand;
m11=(1-m00-m01-m02-m03-m10)*rand;
m12=(1-m00-m01-m02-m03-m10-m11)*rand;
m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1

n0=rand;
n1=(1-n0)*rand;
n2=(1-n0-n1)*rand;
n3=(1-n0-n1-n2)*rand;
n4=1-n0-n1-n2-n3; %normalize to 1

Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
Y = HH14D(t, Ifunc, Area, Y0);

%% Convergence rate

M00=Y(:,3);
M01=Y(:,4);
M02=Y(:,5);
M03=Y(:,6);
M10=Y(:,7);
M11=Y(:,8);
M12=Y(:,9);
M13=Y(:,10);

N0=Y(:,11);
N1=Y(:,12);
N2=Y(:,13);
N3=Y(:,14);
N4=Y(:,15);
% R: 14D to 4D
m=(M11+M01)/3+2*(M12+M02)/3+M13+M03;
h=M10+M11+M12+M13;
n=N1/4+N2/2+3*N3/4+N4;

NaBar = [(1-m).^3.*(1-h) , ...
        3*(1-m).^2.*m.*(1-h) , ...
        3*(1-m).*m.^2.*(1-h) , ...
        m.^3.*(1-h) , ...
        (1-m).^3.*h , ...
        3*(1-m).^2.*m.*h , ...
        3*(1-m).*m.^2.*h , ...
        m.^3.*h];

KBar  = [(1-n).^4 , 4.*n.*(1-n).^3 , 6*n.^2.*(1-n).^2  , 4*n.^3.*(1-n)  , n.^4];

HR_Y=[NaBar,KBar];

Diff_HR=Y(:,3:15)-HR_Y;

Diff=nan(length(t),1);
for i=1:length(t)
    Diff(i)=norm(Diff_HR(i,:));
end

disp(j)

left=log(1./Diff(2:end))./t(2:end)';
left1=log(Diff(1)./Diff(2:end))./t(2:end)';

figure(111)
% hold all

subplot(2,1,1)
plot(t,Diff,'--')
hold on

subplot(2,1,2)
plot(t,log(Diff),'--')
hold on

% figure(1)
% plot(t(2:end),left)
% 
% figure(2)
% plot(t(2:end),left1)
end

%% maxmium distance

disp('sample max distance traces')

for jj=1:10
    V0 = (V_max-V_min).*rand + V_min;
    m00=0.5;
    m01=0;
    m02=0;
    m03=0.5;
    m10=0;
    m11=0;
    m12=0;
    m13=0;

    n0=0.5;
    n1=0;
    n2=0;
    n3=0;
    n4=0.5;
    % R: 14D to 4D
    % m=(M11+M01)/3+2*(M12+M02)/3+M13+M03;
    % h=M10+M11+M12+M13;
    % n=N1/4+N2/2+3*N3/4+N4;

    Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
    Y = HH14D(t, Ifunc, Area, Y0);
    % Convergence rate

    M00=Y(:,3);
    M01=Y(:,4);
    M02=Y(:,5);
    M03=Y(:,6);
    M10=Y(:,7);
    M11=Y(:,8);
    M12=Y(:,9);
    M13=Y(:,10);

    N0=Y(:,11);
    N1=Y(:,12);
    N2=Y(:,13);
    N3=Y(:,14);
    N4=Y(:,15);
    % R: 14D to 4D
    m=(M11+M01)/3+2*(M12+M02)/3+M13+M03;
    h=M10+M11+M12+M13;
    n=N1/4+N2/2+3*N3/4+N4;

    NaBar = [(1-m).^3.*(1-h) , ...
            3*(1-m).^2.*m.*(1-h) , ...
            3*(1-m).*m.^2.*(1-h) , ...
            m.^3.*(1-h) , ...
            (1-m).^3.*h , ...
            3*(1-m).^2.*m.*h , ...
            3*(1-m).*m.^2.*h , ...
            m.^3.*h];

    KBar  = [(1-n).^4 , 4.*n.*(1-n).^3 , 6*n.^2.*(1-n).^2  , 4*n.^3.*(1-n)  , n.^4];

    HR_Y=[NaBar,KBar];

    Diff_HR=Y(:,3:15)-HR_Y;

    Diff=nan(length(t),1);
    for i=1:length(t)
        Diff(i)=norm(Diff_HR(i,:));
    end



    left=log(1./Diff(2:end))./t(2:end)';
    left1=log(Diff(1)./Diff(2:end))./t(2:end)';

    figure(111)
    subplot(2,1,1)
    plot(t,Diff,'--','LineWidth',0.5)

    subplot(2,1,2)
    plot(t,log(Diff),'--','LineWidth',0.5)
    
end
%% plot the maximum

subplot(2,1,1)
plot(t,(Diff(1)).*exp(lambda_decay*t),'r-','LineWidth',2)
set(gca,'Fontsize',16)
grid on

subplot(2,1,2)
plot(t,log(Diff(1))+lambda_decay*t,'r-','LineWidth',2)
set(gca,'Fontsize',16)
grid on


