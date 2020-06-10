function [Y,T_run] = HHSS14D_skip(t, Ifunc, Area)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM_HH.m
% Written by Shusen Pu and Peter Thomas
% Oct, 2018

%%% Inputs
% t is vector of time values (ms)
% Ifunc is a function @(t)f(t) that returns stimulus value as a function of time (in ms)
% Area: Membrane Area (mu m^2)

% Noise Model (String), possible values: 
%   None: 'ODE'
%   EM: 'Euler Maruyama', must also have value for Area

%%% Outputs
% Y(:,1) : t
% Y(:,2) : V
% Y(:,3:10) : fractions in Na states m00-m13 
% Y(:,11:15) : fraction in K states n0-n4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize quantities needed to run solver 
t_ind=tic;
% time step size
dt = t(2)-t(1);

% Number of time steps
nt = length(t);  % total
nt1 = nt-1;  % at which to solve

NaNoise1=randn(20,nt);
KNoise1=randn(8,nt);

% Initial Values
% from XPP code on the limit cycle
t0 = t(1);
V0 = -61.897274; 
m00=0.4329406;m01=0.10034765;m02=0.0077529117;m03=0.00019966466;
m10=0.36696321;m11=0.085055299;m12=0.0065714167;
% m13=0.00016923703;
m13=1-m00-m01-m02-m03-m10-m11-m12; %normalize to 1
n0=0.13761529;n1=0.35331318;n2=0.34016079;n3=0.14555468;
%n4=0.023356047;
n4=1-n0-n1-n2-n3; %normalize to 1

 
% Initialize Output
Y = zeros(nt,14); 
Y(1,1) = V0;

Na_gates =[m00;m01;m02;m03;m10;m11;m12;m13];
K_gates=[n0,n1,n2,n3,n4]';
Y(1,1) = V0;

Na_gates =Na_gates(:);
K_gates=K_gates(:);

Y(1,2:9)=Na_gates';
Y(1,10:14)=K_gates';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Values

% Number of Channels
NNa = round(Area*60); % Na
NK = round(Area*18); % K

% Capacitance
C = 1; % muF /cm^2

% Na Current
gNa = 120; % mS/cm^2
ENa = 50; % mV

% K Current
gK = 36; % mS/cm^2
EK = -77; % mV

% Passive Leak
gL = 0.3; % mS / cm^2
EL = -54.4; % mV



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Which Noise Model and Do Some Necessary Setup


% Conductance Noise Model using Euler Maruyama (Thomas and Pu)
    
    % Drift Na
    ANa = @(V) ...
         [ -3*alpham(V)-alphah(V)       , betam(V)                , 0              , 0                      ,  betah(V)                 , 0                       , 0          , 0 ;
            3*alpham(V)               ,-2*alpham(V)-betam(V)-alphah(V), 2*betam(V)        , 0                      ,  0                     , betah(V)                   , 0          , 0 ;
            0                      , 2*alpham(V)             , -alpham(V)-2*betam(V)-alphah(V),  3*betam(V)        ,  0                     ,  0                      , betah(V)      , 0 ;
            0                      , 0                    , alpham(V)         , -3*betam(V)-alphah(V)        ,  0                     ,  0                      , 0          , betah(V)    ;
            alphah(V)                 , 0                    , 0              , 0                      ,  -3*alpham(V) - betah(V)     , betam(V)                   , 0          , 0 ;
            0                      , alphah(V)               , 0              , 0                      ,  3*alpham(V)              ,  -2*alpham(V)-betam(V)-betah(V)  ,   2*betam(V)  , 0 ;
            0                      , 0                    , alphah(V)         , 0                      ,  0                     ,  2*alpham(V)               ,   -alpham(V)-2*betam(V)-betah(V) , 3*betam(V)  ;
            0                      , 0                    , 0              , alphah(V)                 ,           0            ,  0                      ,  alpham(V)    , -3*betam(V)-betah(V)];

    % Drift K
    AK = @(V) ...
        [-4*alphan(V), betan(V)             , 0                , 0                  , 0
          4*alphan(V), -3*alphan(V)-betan(V),  2*betan(V)               , 0,                   0;
          0,        3*alphan(V),        -2*alphan(V)-2*betan(V), 3*betan(V),          0;
          0,        0,               2*alphan(V),          -alphan(V)-3*betan(V), 4*betan(V);
          0,        0,               0,                 alphan(V),          -4*betan(V)];
    
    
    %Noise Term: adding white noise to each transaction
    %defined in two functions namely, DNafull and DKfull 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HERE IS THE SOLVER            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% USING EULER FOR ODEs,         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Milstein FOR SDEs,            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:nt

  % Input Current
  I = Ifunc(t(i-1));
  % update the gates
   Na_gates1 = Na_gates + dt*ANa(V0)*Na_gates+ sqrt(dt)*DNafull(V0,Na_gates,NNa,NaNoise1(:,i-1));
   
   %check for out of bounds
   j1=0;j2=0;
   while sum(Na_gates1<0)>0
       new_ran=randn(20,1);
       j1=j1+1;
       Na_gates1 = Na_gates + dt*ANa(V0)*Na_gates+ sqrt(dt)*DNafull(V0,Na_gates,NNa,new_ran);
       disp(i)
       disp(j1)
       if j1>50
           keyboard
       end
   end
   
    if ~isreal(K_gates)
        keyboard
    end
    
   K_gates1 =  K_gates + dt*AK(V0) *K_gates + sqrt(dt)*DKfull(V0,K_gates,NK,KNoise1(:,i-1)) ;
   
   if ~isreal(K_gates1)
        keyboard
    end
   
   while sum(K_gates1<0)>0
       new_ran1=randn(8,1);
       j2=j2+1;
       disp(i)
       disp(j2)
       K_gates1 = K_gates + dt*AK(V0) *K_gates + sqrt(dt)*DKfull(V0,K_gates,NK,new_ran1);
       if j2>50
           keyboard
       end
   end
    
  % Compute Fraction of open channels
  NaFraction= Na_gates1(8);
  KFraction=K_gates1(5);
      
  % Update Voltage
  Vrhs = (-gNa*(NaFraction)*(V0 - ENa)-gK*(KFraction)*(V0 - EK) - gL*(V0-EL) + I)/C;
  V0 = V0 + dt*Vrhs ;   % VNoise is non-zero for Current Noise Model
  
  % update gates
  Na_gates=Na_gates1;
  K_gates= K_gates1;
%   
%   if i==437796
%       keyboard
%   end
%   
  
  if ~isreal(V0)
      
     keyboard
     
  end
      
  
  % Save Outputs    
  Y(i,1) = V0;
  Na_gates =Na_gates(:);
  K_gates=K_gates(:);

  Y(i,2:9)=Na_gates';
  Y(i,10:14)=K_gates';
end % End loop over time for SDE solver

T_run=toc(t_ind);
end




 % End Function Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF SOLVER                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define functions used above

% subunit kinetics (Hodgkin and Huxley parameters, modified from XPP code by PJT OCt, 2018)
function out = alpham(V)
  out =  .1*(V+40)/(1-exp(-(V+40)/10));
  %0.1 * (25-V)/ (exp((25-V)/10)-1);
end

function out = betam(V)
  out = 4*exp(-(V+65)/18);
  %4 * exp(-V/18);
end

function out = alphah(V)
  out =  .07*exp(-(V+65)/20);
  %0.07 * exp(-V/20);
end

function out = betah(V)
  out = 1/(1+exp(-(V+35)/10));
  %1/ (exp((30-V)/10)+1);
end

function out = alphan(V)
  out = .01*(V+55)/(1-exp(-(V+55)/10));
  %0.01 * (10-V) / (exp((10-V)/10)-1);
end

function out = betan(V)
  out = .125*exp(-(V+65)/80);
  %0.125 * exp(-V/80);
end

% Full Model (Thomas and Pu)

% Diffusion matrix for Na (Thomas and Pu)
function   D = DNafull(V,Y,N,M) 
%V is voltage, Y stands for probablity at each state (we use Y-Ybar)
%N is the dimension or number of Na Channels, M is randn 20 by 1
D = zeros(8,1);
% Y=Y.*(Y>0); %to exclude negative numbers added by SP&PJT 10/08/2018
% Y=abs(Y);
y00 = Y(1);
y10 = Y(2);
y20 = Y(3);
y30 = Y(4);
y01 = Y(5);
y11 = Y(6);
y21 = Y(7);
y31 = Y(8);

% noise for these 20 transactions
N20=zeros(20,1);
%save all these alpha and beta values for computational efficiency
ah=alphah(V);
bh=betah(V);
am=alpham(V);
bm=betam(V);
N20(1)=sqrt(ah*y00)*M(1); %00--01
N20(2)=sqrt(bh*y01)*M(2); %01--00
N20(3)=sqrt(3*am*y00)*M(3); %00--10
N20(4)=sqrt(bm*y10)*M(4); %10--00
N20(5)=sqrt(ah*y10)*M(5);  %10--11
N20(6)=sqrt(bh*y11)*M(6);  %11--10
N20(7)=sqrt(2*am*y10)*M(7);%10--20
N20(8)=sqrt(2*bm*y20)*M(8); %20--10
N20(9)=sqrt(ah*y20)*M(9); %20--21
N20(10)=sqrt(bh*y21)*M(10); %21--20
N20(11)=sqrt(am*y20)*M(11); %20--30
N20(12)=sqrt(3*bm*y30)*M(12);%30--20
N20(13)=sqrt(ah*y30)*M(13); %30--31
N20(14)=sqrt(bh*y31)*M(14); %31--30
N20(15)=sqrt(3*am*y01)*M(15);%01--11
N20(16)=sqrt(bm*y11)*M(16); %11--01
N20(17)=sqrt(2*am*y11)*M(17);%11--21
N20(18)=sqrt(2*bm*y21)*M(18);%21--11
N20(19)=sqrt(am*y21)*M(19);%21--31
N20(20)=sqrt(3*bm*y31)*M(20);%31--21

D(1) = -N20(1)-N20(3)+N20(2)+N20(4) ;
D(2) = -N20(4)-N20(5)-N20(7)+N20(3)+N20(6)+N20(8);
D(3) = -N20(8)-N20(9)-N20(11)+N20(7)+N20(10)+N20(12);
D(4) = -N20(12)-N20(13)+N20(11)+N20(14);
D(5) = -N20(2)-N20(15)+N20(1)+N20(16);
D(6) = -N20(6)-N20(16)-N20(17)+N20(5)+N20(15)+N20(18);
D(7) = -N20(10)-N20(18)-N20(19)+N20(9)+N20(17)+N20(20);
D(8) = -N20(14)-N20(20)+N20(13)+N20(19);

D = D/sqrt(N);
end

% Diffusion matrix for K (Thomas and Pu)
function   D = DKfull(V,Y,N,M) 
%V is voltage, Y stands for probablity at each state (we use Y-Ybar)
%N is the dimension or number of K Channels, M is randn 8 by 1
D = zeros(5,1);

% noise for these 8 transactions
N8=zeros(8,1);
%save all these alpha and beta values for computational efficiency
an=alphan(V);
bn=betan(V);
% Y=Y.*(Y>0); %to exclude negative numbers added by SP&PJT 10/08/2018
% Y=abs(Y);

N8(1)=sqrt(4*an*Y(1))*M(1); %0--1
N8(2)=sqrt(bn*Y(2))*M(2); %1--0
N8(3)=sqrt(3*an*Y(2))*M(3); %1--2
N8(4)=sqrt(2*bn*Y(3))*M(4); %2--1
N8(5)=sqrt(2*an*Y(3))*M(5); %2--3
N8(6)=sqrt(3*bn*Y(4))*M(6); %3--2
N8(7)=sqrt(an*Y(4))*M(7); %3--4
N8(8)=sqrt(4*bn*Y(5))*M(8); %4--3


D(1) = -N8(1)+N8(2);
D(2) = -N8(2)-N8(3)+N8(1)+N8(4);
D(3) = -N8(4)-N8(5)+N8(3)+N8(6);
D(4) = -N8(6)-N8(7)+N8(5)+N8(8);
D(5) = -N8(8)+N8(7);

D = D/sqrt(N);
end
