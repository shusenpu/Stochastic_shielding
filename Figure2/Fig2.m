
% Time
dt=0.001;
t0=0;
tmax=100;
t=t0:dt:tmax;
nt=length(t)-1;
NaNoise1=randn(20,nt);
KNoise1=randn(8,nt);
Area=100;
%Area=1;
Ifunc=@(t)10;

% Choose which variables to plot:
plotidx1=2; % voltage
%plotidx2=10; % m13
plotidx2=11; % n0
plotidx3=12; % n1 % maybe this one?
%plotidx2=13; % n2
%plotidx2=14; % n3
% plotidx3=15; % n4
p=[plotidx1,plotidx2,plotidx3];
%Y(:,2) %voltage
%Y(:,3) %m00
%Y(:,4) %m01
%Y(:,5) %m02
%Y(:,6) %m03
%Y(:,7) %m10
%Y(:,8) %m11
%Y(:,9) %m12
%Y(:,10) %m13
%Y(:,11) %n0
%Y(:,12) %n1
%Y(:,13) %n2
%Y(:,14) %n3
%Y(:,15) %n4

figure(1)
[vv,nn]=meshgrid(-80:5:50,0:.01:1);
if p==[2,14,15]
    h=mesh(vv,4*nn.^3.*(1-nn),nn.^4);
elseif p==[2,13,15]
    h=mesh(vv,6*nn.^2.*(1-nn).^2,nn.^4);
elseif p==[2,12,15]
    h=mesh(vv,4*nn.*(1-nn).^3,nn.^4);
elseif p==[2,11,15]
    h=mesh(vv,(1-nn).^4,nn.^4);
elseif p==[2,11,12]
    h=mesh(vv,(1-nn).^4,4*(1-nn).^3.*nn);
end
set(h,'FaceAlpha',0.1)
set(h,'EdgeColor','k')
set(h,'FaceColor','k')
hold on

for ic=1:17
    
    if ic==1 % NOT on MN manifold
        V0=20;
        m00=.5;
        m01=0;
        m02=0;
        m03=0;
        m10=0;
        m11=0;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=.5;
        n1=0;
        n2=0;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
       plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'g-','LineWidth',2)
%         plot3(Y(5:end,plotidx1),Y(5:end,plotidx2),Y(5:end,plotidx3),'c-','LineWidth',2)
        hold on
        plot3(Y(1,plotidx1),Y(1,plotidx2),Y(1,plotidx3),'g*','MarkerSize',8,'LineWidth',4)
    end
    
    if ic==2 % NOT on MN manifold. Stochastic
        V0=20;
        m00=.5;
        m01=0;
        m02=0;
        m03=0;
        m10=0;
        m11=0;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=.5;
        n1=0;
        n2=0;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HHSto15D(t, Ifunc, Area, NaNoise1, KNoise1,Y0);
        %[Y] = HH14D(t, Ifunc, Area, Y0);
        plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'r-','LineWidth',2)
        %plot3(Y(5:end,plotidx1),Y(5:end,plotidx2),Y(5:end,plotidx3),'r-','LineWidth',1.5)
        hold on
        plot3(Y(1,plotidx1),Y(1,plotidx2),Y(1,plotidx3),'r+','MarkerSize',6,'LineWidth',4)
    end
    
    if ic==3 % on MN manifold
        V0=50;
        m00=0;
        m01=0;
        m02=0;
        m03=0;
        m10=0.5;
        m11=0.5;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=0;
        n1=0.2;
        n2=0.5;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
       % plot3(Y(45:end,plotidx1),Y(45:end,plotidx2),Y(45:end,plotidx3),'b-','LineWidth',2)
        hold on
    end
    
    if ic==1+8 % on MN manifold
        V0=50;
        m00=0;
        m01=0;
        m02=0;
        m03=0;
        m10=0;
        m11=0;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=1;
        n1=0;
        n2=0;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
        plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'b-','LineWidth',2)
        hold on
        plot3(Y(1,plotidx1),Y(1,plotidx2),Y(1,plotidx3),'b+','MarkerSize',16)
     
    end
    if ic==2+8 % on MN manifold Keep this one
        V0=-70;
        m00=1;
        m01=0;
        m02=0;
        m03=0;
        m10=0;
        m11=0;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=1;
        n1=0;
        n2=0;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
        plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'b-','LineWidth',2)
        hold on
        plot3(Y(1,plotidx1),Y(1,plotidx2),Y(1,plotidx3),'b+','MarkerSize',8,'LineWidth',2)
     
    end
    if ic==3+8 % on MN manifold OK
        V0=20;
        m00=0;
        m01=0;
        m02=0;
        m03=0;
        m10=0;
        m11=0;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=0;
        n1=0;
        n2=0;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
        %plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'b-','LineWidth',2)
        hold on
    end
    if ic==6+8 % on MN manifold OK
        V0=30;
        m00=1;
        m01=0;
        m02=0;
        m03=0;
        m10=0;
        m11=0;
        m12=0;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=1;
        n1=0;
        n2=0;
        n3=0;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
        plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'b-','LineWidth',2)
        hold on
        plot3(Y(1,plotidx1),Y(1,plotidx2),Y(1,plotidx3),'b+','MarkerSize',8,'LineWidth',2)
     
    end
    if ic==17 % Initial condition on LC
        V0 = -61.897274;
        m00=0.4329406;
        m01=0.10034765;
        m02=0.0077529117;
        m03=0.00019966466;
        m10=0.36696321;
        m11=0.085055299;
        m12=0.0065714167;
        % m13=0.00016923703;
        m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
        n0=0.13761529;
        n1=0.35331318;
        n2=0.34016079;
        n3=0.14555468;
        %n4=0.023356047;
        n4=1-n0-n1-n2-n3; %normalize to 1
        Y0=[t(1),V0,m00,m01,m02,m03,m10,m11,m12,m13,n0,n1,n2,n3,n4];
        % Find trajectory
        [Y] = HH14D(t, Ifunc, Area, Y0);
        plot3(Y(:,plotidx1),Y(:,plotidx2),Y(:,plotidx3),'k--','LineWidth',3)
        plot3(Y(1,plotidx1),Y(1,plotidx2),Y(1,plotidx3),'k+','MarkerSize',8,'LineWidth',2)
     
    end
end
rotate3d on
grid on
set(gca,'FontSize',20)
xlabel('Voltage (mV)')
ylabel('n_0')
zlabel('n_1')
% alternative viewpoints
%view(-39.1000,13.2000)
%view(-132,-4.4)
% view(33.6,-23.2)
% view(67.8589,54.1642)
% view(-117.8189,2.5802)
% Adjust limits to maximize aesthetic impact
% ylim([0,.2])
% zlim([0,.5])
% view(-145.359,-54.8)
% view(-55.8823,39.3426)
% view( -148.2093,-44.9584)
view(-161.4713,-48.7695)

ylim([0, 0.5])