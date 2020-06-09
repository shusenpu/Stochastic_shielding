function [Y] = HH14D(t, Ifunc, Area, Y0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HH14D.m
% Written by Shusen Pu and Peter Thomas
% June, 2020


%%% Inputs
% t is vector of time values (ms)
% Ifunc is a function @(t)f(t) that returns stimulus value as a function of time (in ms)
% Area: Membrane Area (mu m^2)

%%% Outputs
% Y(:,1) : t
% Y(:,2) : V
% Y(:,3:10) : fractions in Na states m00-m13 
% Y(:,11:15) : fraction in K states n0-n4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize quantities needed to run solver 

% time step size
dt = t(2)-t(1);

% Number of time steps
nt = length(t);  % total
nt1 = nt-1;  % at which to solve

% Initial Values
% from XPP code on the limit cycle
t0 = t(1);

% Initialize Output
Y = zeros(nt,15); 
Y(1,:) = Y0;

Na_gates =Y0(3:10)';
K_gates=Y0(11:15)';
V0 =Y0(2);
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
    
    
for i=2:nt

  % Input Current
  I = Ifunc(t(i-1));
  % update the gates

  Na_gates = Na_gates + dt*ANa(V0)*Na_gates;
  K_gates =  K_gates + dt*AK(V0) *K_gates;
    
    
  % Compute Fraction of open channels
  NaFraction= Na_gates(8);
  KFraction=K_gates(5);
      
  % Update Voltage
  Vrhs = (-gNa*(NaFraction)*(V0 - ENa)-gK*(KFraction)*(V0 - EK) - gL*(V0-EL) + I)/C;
  V = V0 + dt*Vrhs ;   % VNoise is non-zero for Current Noise Model
  
  % Save Outputs    
  Y(i,1) = t(i);
  Y(i,2) = V;
  Y(i,3:10) = Na_gates';
  Y(i,11:15) = K_gates';
  % Keep "old values" to use in next Euler time step
  V0 = V;
end

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

