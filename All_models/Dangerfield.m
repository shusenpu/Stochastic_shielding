function [Y,T_run] = Dangerfield(t, Ifunc, Area)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dangerfield.m
% Written by Shusen Pu
% Sept, 2019

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
cont=0;
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
Y = zeros(nt,1); 
Y(1) = V0;

Na_gates =[m00;m01;m02;m03;m10;m11;m12;m13];
K_gates=[n0,n1,n2,n3,n4]';
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
   Na_gates = Na_gates + dt*ANa(V0)*Na_gates+ sqrt(dt)*DNaOrio(V0,Na_gates,NNa,NaNoise1(:,i-1));
   K_gates =  K_gates + dt*AK(V0) *K_gates + sqrt(dt)*DKOrio(V0,K_gates,NK,KNoise1(:,i-1)) ;
 
   % Check to see if point lies in correct domain
    if min(Na_gates)<0  || max(Na_gates)>1         
        % Project into the correct domain
        Na_gates = projsplx(Na_gates);
        cont=cont+1;
    else       
        % Point lies in correct domain so no action needed           
    end
    
     % Check to see if point lies in correct domain
    if min(K_gates)<0  || max(K_gates)>1         
        % Project into the correct domain
        K_gates = projsplx(K_gates);
        cont=cont+1;
    else       
        % Point lies in correct domain so no action needed           
    end
    
    
    
  % Compute Fraction of open channels
  NaFraction= Na_gates(8);
  KFraction=K_gates(5);
      
  % Update Voltage
  Vrhs = (-gNa*(NaFraction)*(V0 - ENa)-gK*(KFraction)*(V0 - EK) - gL*(V0-EL) + I)/C;
  V0 = V0 + dt*Vrhs ;   % VNoise is non-zero for Current Noise Model
  
  % Save Outputs    
  Y(i) = V0;
  
end  % End loop over time for SDE solver

T_run=toc(t_ind);

end % End Function Definition
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

% Diffusion matrix for Na (Thomas and Pu) for Dangerfield 2012 methods
function   D = DNaOrio(V,Y,N,M) 
%V is voltage, Y stands for probablity at each state (we use Y-Ybar)
%N is the dimension or number of Na Channels, M is randn 20 by 1
D = zeros(8,1);
% Y=Y.*(Y>0); %to exclude negative numbers added by SP&PJT 10/08/2018
y00 = Y(1);
y10 = Y(2);
y20 = Y(3);
y30 = Y(4);
y01 = Y(5);
y11 = Y(6);
y21 = Y(7);
y31 = Y(8);

% noise for these 10 paris of transactions 
N10=zeros(10,1);
%save all these alpha and beta values for computational efficiency
ah=alphah(V);
bh=betah(V);
am=alpham(V);
bm=betam(V);


N10(1)=sqrt((3*ah*y00+bh*y01))*M(1); %00--01
N10(2)=sqrt((3*am*y00+bm*y10))*M(2); %00--10
N10(3)=sqrt((ah*y10+bh*y11))*M(3);  %10--11
N10(4)=sqrt((2*am*y10+2*bm*y20))*M(4);%10--20
N10(5)=sqrt((ah*y20+bh*y21))*M(5); %20--21
N10(6)=sqrt((am*y20+3*bm*y30))*M(6); %20--30
N10(7)=sqrt((ah*y30+bh*y31))*M(7); %30--31
N10(8)=sqrt((3*am*y01+bm*y11))*M(8);%01--11
N10(9)=sqrt((2*am*y11+2*bm*y21))*M(9);%11--21
N10(10)=sqrt((am*y21+3*bm*y31))*M(10);%21--31

D(1) = N10(1)+N10(2);
D(2) = -N10(2)+N10(3)+N10(4);
D(3) = -N10(4)+N10(5)+N10(6);
D(4) = -N10(6)+N10(7);
D(5) = -N10(1)+N10(8);
D(6) = -N10(3)-N10(8)+N10(9);
D(7) = -N10(5)-N10(9)+N10(10);
D(8) = -N10(7)-N10(10);

D = D/sqrt(N);
end

% Diffusion matrix for K (Orio and Soudry)
function   D = DKOrio(V,Y,N,M) 
%V is voltage, Y stands for probablity at each state
%N is the dimension or number of K Channels, M is randn 4 by 1
D = zeros(5,1);

% noise for these 8 transactions
N4=zeros(4,1);
%save all these alpha and beta values for computational efficiency
an=alphan(V);
bn=betan(V);

N4(1)=sqrt((4*an*Y(1)+bn*Y(2)))*M(1); %0--1
N4(2)=sqrt((3*an*Y(2)+2*bn*Y(3)))*M(2); %1--2
N4(3)=sqrt((2*an*Y(3)+3*bn*Y(4)))*M(3); %2--3
N4(4)=sqrt((an*Y(4)+4*bn*Y(5)))*M(4); %3--4

D(1) = N4(1);
D(2) = -N4(1)+N4(2);
D(3) = -N4(2)+N4(3);
D(4) = -N4(3)+N4(4);
D(5) = -N4(4);

D = D/sqrt(N);
end


function x = projsplx(y)
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}

% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.

m = length(y); bget = false;

s = sort(y,'descend'); tmpsum = 0;

for ii = 1:m-1
    tmpsum = tmpsum + s(ii);
    tmax = (tmpsum - 1)/ii;
    if tmax >= s(ii+1)
        bget = true;
        break;
    end
end
    
if ~bget
    tmax = (tmpsum + s(m) -1)/m;
end

x = max(y-tmax,0);

end
