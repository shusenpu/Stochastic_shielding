function [Y,T_run] = Goldwyn(t, Ifunc, Area)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Shusen Pu
% Aug 8, 2019

%%% Inputs
% t is vector of time values (ms)
% Ifunc is a function @(t)f(t) that returns stimulus value as a function of time (in ms)
% Area: Membrane Area (mu m^2)

%%% Outputs
% Y(:,1) : t
% Y(:,2) : V
% Y(:,3) : fraction open Na channels
% Y(:,4) : fraction open K channels
% Y(:,5) : m
% Y(:,6) : h
% Y(:,7) : n


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize quantities needed to run solver 
% time step size
t_ind=tic;
dt = t(2)-t(1);

% Number of time steps
nt = length(t);  % total
nt1 = nt-1;  % at which to solve

% Initial Values
t0 = t(1);
V0 = -61.897274; 
m0 = alpham(V0) / (alpham(V0) + betam(V0)); % m
h0 = alphah(V0) / (alphah(V0) + betah(V0)); % h
n0 = alphan(V0) / (alphan(V0) + betan(V0)); % n
NaFraction = m0^3*h0;
KFraction = n0^4;
 
% Initialize Output
Y = zeros(nt,7); 
Y(1,:) = [t0, V0, m0^3*h0, n0^4, m0, h0, n0];


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

% Conductance Noise (FL Channel Model)

    NaHat = zeros(8,1);  %Initial values set to 0
    KHat = zeros(5,1); %Initial values set to 0
    NaNoise = randn(8, nt1);
    KNoise = randn(5, nt1);
    
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
    
    
    % Diffusion Na : Defined in a afunction below
    
    
    % Diffusion K
    DK = @(V,X) (1/(NK)) * ...
        [   (4*alphan(V)*X(1) + betan(V)*X(2)) ,  -(4*alphan(V)*X(1) + betan(V)*X(2))                                    ,   0                                                                  , 0                                                                , 0 ;
           -(4*alphan(V)*X(1) + betan(V)*X(2)),  (4.*alphan(V)*X(1) + (3*alphan(V)+ betan(V))*X(2) + 2.*betan(V)*X(3))  ,  -(2*betan(V)*X(3) + 3*alphan(V)*X(2) )                               , 0                                                                , 0 ; 
           0                                  ,  -(2*betan(V)*X(3) + 3*alphan(V)*X(2))                                   , (3*alphan(V)*X(2) + (2*alphan(V)+2*betan(V))*X(3) + 3*betan(V)*X(4)) , -(3*betan(V)*X(4) + 2*alphan(V)*X(3))                            , 0 ;
           0                                  ,  0                                                                      , -(3*betan(V)*X(4) + 2*alphan(V)*X(3))                                , (2*alphan(V)*X(3) + (alphan(V)+3*betan(V))*X(4) +4*betan(V)*X(5)), -(4*betan(V)*X(5) + alphan(V)*X(4)) ;
           0                                  ,  0                                                                      , 0                                                                    , -(4*betan(V)*X(5) + alphan(V)*X(4))                              , (alphan(V)*X(4) + 4*betan(V)*X(5)) ];


    % Take Matrix square roots numerically using SVD
    SNa = @(V,Y,NNa) mysqrtm(DNa(V,Y,NNa));
    SK = @(V,X) mysqrtm(DK(V,X)); 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HERE IS THE SOLVER            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% USING EULER FOR ODEs,         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EULER-MARUYAMA FOR SDEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:nt
    I = Ifunc(t(i-1));    
    
    % Update subunits
  % Noise terms are non-zero for Subunit Noise model
  m = m0 + dt*(alpham(V0)*(1-m0) - betam(V0)*m0) ;  % shifted to i-1 in function
  h = h0 + dt*(alphah(V0)*(1-h0) - betah(V0)*h0) ;
  n = n0 + dt*(alphan(V0)*(1-n0) - betan(V0)*n0) ;
    
    
  NaBar = [(1-m0)^3*(1-h0) , 3*(1-m0)^2*m0*(1-h0) , 3*(1-m0)*m0^2*(1-h0) , m0^3*(1-h0) , (1-m0)^3*h0 , 3*(1-m0)^2*m0*h0 , 3*(1-m0)*m0^2*h0 , m0^3*h0];
  KBar  = [(1-n0)^4 , 4*n0*(1-n0)^3 , 6*n0^2*(1-n0)^2  , 4*n0^3*(1-n0)  , n0^4];
  NaHat = NaHat + dt*ANa(V0)*NaHat + sqrt(dt)*SNa(V0,NaBar,NNa)*NaNoise(:,i-1);
  KHat =  KHat  + dt*AK(V0) *KHat  + sqrt(dt)*SK(V0,KBar)*KNoise(:,i-1);
  NaFluctuation = NaHat(end) ;
  KFluctuation  = KHat(end) ;
        
    % Note: Impose bounds on fractions to avoid <0 or >1 in dV/dt equation, this doesn't directly alter the dynamics of the subunits or channels
    NaFraction = max(0, min(1, m0^3*h0 + NaFluctuation));  % Fluctuations are non-zero for Conductance Noise Models
    KFraction = max(0, min(1, n0^4 + KFluctuation));
  
  % Update Voltage
  Vrhs = (-gNa*(NaFraction)*(V0 - ENa)-gK*(KFraction)*(V0 - EK) - gL*(V0-EL) + I)/C;
  V = V0 + dt*Vrhs;   % VNoise is non-zero for Current Noise Model
  
  % Save Outputs  
  Y(i,1) = t(i);
  Y(i,2) = V;
  Y(i,3) = NaFraction;
  Y(i,4) = KFraction;
  Y(i,5) = m;
  Y(i,6) = h;
  Y(i,7) = n;

  % Keep "old values" to use in next Euler time step
  V0 = V;
  m0 = m;
  h0 = h;
  n0 = n;
  
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

% Computing matrix square roots with SVD
function S = mysqrtm(D)
 [u,s,v] = svd(D);
 S = u*sqrt(s)*v';
end

% Diffusion matrix for Na
function   D = DNa(V,Y,N)
D = zeros(8,8);
y00 = Y(1);
y10 = Y(2);
y20 = Y(3);
y30 = Y(4);
y01 = Y(5);
y11 = Y(6);
y21 = Y(7);
y31 = Y(8);

D(1,1) = (  (3*alpham(V) + alphah(V))*y00 + betam(V)*y10 + betah(V)*y01 )  ;
D(1,2) = (-3*alpham(V)*y00 - betam(V)*y10); 
D(1,3) = 0;
D(1,4) = 0;
D(1,5) = - (alphah(V)*y00 + betah(V)*y01);
D(1,6) = 0;
D(1,7) = 0;
D(1,8) = 0;

D(2,1) = D(1,2);
D(2,2) = ((betam(V)+2*alpham(V))*y10 + 2*betam(V)*y20 + 3*alpham(V)*y00 + alphah(V)*y10 + betah(V)*y11) ;
D(2,3) = -(2*alpham(V)*y10 + 2*betam(V)*y20); 
D(2,4) = 0;
D(2,5) = 0;
D(2,6) = -(alphah(V)*y10 + betah(V)*y11);
D(2,7) = 0;
D(2,8) = 0;

D(3,1) = D(1,3);
D(3,2) = D(2,3);
D(3,3) = ((2*betam(V)+alpham(V))*y20 + 3*betam(V)*y30 + 2*alpham(V)*y10 + alphah(V)*y20 + betah(V)*y21) ;
D(3,4) = -(alpham(V)*y20+3*betam(V)*y30); 
D(3,5) = 0;
D(3,6) = 0;
D(3,7) = -(alphah(V)*y20+betah(V)*y21);
D(3,8) = 0;

D(4,1) = D(1,4);
D(4,2) = D(2,4);
D(4,3) = D(3,4);
D(4,4) = (3*betam(V)*y30 + alpham(V)*y20 + alphah(V)*y30 + betah(V)*y31); 
D(4,5) = 0;
D(4,6) = 0;
D(4,7) = 0;
D(4,8) = -(alphah(V)*y30 + betah(V)*y31); 

D(5,1) = D(1,5);
D(5,2) = D(2,5);
D(5,3) = D(3,5);
D(5,4) = D(4,5);
D(5,5) = (3*alpham(V)*y01 + betam(V)*y11 + betah(V)*y01 + alphah(V)*y00) ;
D(5,6) = -(3*alpham(V)*y01 + betam(V)*y11) ;
D(5,7) = 0;
D(5,8) = 0;

D(6,1) = D(1,6);
D(6,2) = D(2,6);
D(6,3) = D(3,6);
D(6,4) = D(4,6);
D(6,5) = D(5,6);
D(6,6) = ((betam(V)+2*alpham(V))*y11 + 2*betam(V)*y21 + 3*alpham(V)*y01 + betah(V)*y11 + alphah(V)*y10);
D(6,7) = -(2*alpham(V)*y11+2*betam(V)*y21);
D(6,8) = 0;

D(7,1) = D(1,7);
D(7,2) = D(2,7);
D(7,3) = D(3,7);
D(7,4) = D(4,7);
D(7,5) = D(5,7);
D(7,6) = D(6,7);
D(7,7) = ((2*betam(V)+alpham(V))*y21+3*betam(V)*y31+2*alpham(V)*y11+betah(V)*y21+alphah(V)*y20);
D(7,8) = -(alpham(V)*y21+3*betam(V)*y31);

D(8,1) = D(1,8);
D(8,2) = D(2,8);
D(8,3) = D(3,8);
D(8,4) = D(4,8);
D(8,5) = D(5,8);
D(8,6) = D(6,8);
D(8,7) = D(7,8);
D(8,8) = (3*betam(V)*y31 + alpham(V)*y21 + betah(V)*y31 + alphah(V)*y30);

D = D/(N);
end
