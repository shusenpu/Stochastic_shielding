function [Y,T_run] = FoxLu_StoHH(t, Ifunc, Area,seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Shusen Pu


%%% Outputs
% Y(:,1) : t
% Y(:,2) : V
% Y(:,3) : fraction open Na channels
% Y(:,4) : fraction open K channels
% Y(:,5) : m
% Y(:,6) : h
% Y(:,7) : n
    
rng(seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize quantities needed to run solver 

t_ind=tic;

% time step size
dt = t(2)-t(1);

% Number of time steps
nt = length(t);  % total
nt1 = nt-1;  % at which to solve

% Initial Values
t0 = t(1);
V0 = 0; 
m0 = alpham(V0) / (alpham(V0) + betam(V0)); % m
h0 = alphah(V0) / (alphah(V0) + betah(V0)); % h
n0 = alphan(V0) / (alphan(V0) + betan(V0)); % n
 
% Initialize Output
Y = zeros(nt,4); 
Y(1,1) = V0;
Y(1,2) = m0;
Y(1,3) = h0;
Y(1,4) = n0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Values

% Number of Channels
NNa = round(Area*60); % Na
NK = round(Area*18); % K

% Capacitance
C = 1; % muF /cm^2

% Na Current
gNa = 120; % mS/cm^2
ENa = 120; % mV

% K Current
gK = 36; % mS/cm^2
EK = -12; % mV

% Passive Leak
gL = 0.3; % mS / cm^2
EL = 10.6; % mV


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Which Noise Model and Do Some Necessary Setup
    mNoiseVec = randn(nt1,1);
    % Imposing bounds on argument of sqrt functions, not directly altering dynamics of the subunits
   % mNoise = @(V,m,i) sqrt((alpham(V)*(1-m) + betam(V)*m)/NNa) * mNoiseVec(i-1);
    mNoise = @(V,m,i) sqrt(2)*sqrt((alpham(V)*(1-m) + betam(V)*m)/NNa) * mNoiseVec(i-1);
    hNoiseVec = randn(nt1,1);
    %hNoise = @(V,h,i) sqrt((alphah(V)*(1-h) + betah(V)*h)/NNa) * hNoiseVec(i-1);
    hNoise = @(V,h,i) sqrt(2)*sqrt((alphah(V)*(1-h) + betah(V)*h)/NNa) * hNoiseVec(i-1);
    nNoiseVec = randn(nt1,1);
%     nNoise = @(V,n,i) sqrt((alphan(V)*(1-n) + betan(V)*n)/NK)  * nNoiseVec(i-1);
    nNoise = @(V,n,i) sqrt(2)*sqrt((alphan(V)*(1-n) + betan(V)*n)/NK)  * nNoiseVec(i-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HERE IS THE SOLVER            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% USING EULER FOR ODEs,         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EULER-MARUYAMA FOR SDEs, and  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GILLESPIE FOR MARKOV CHAIN    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:nt

  % Input Current
  I = Ifunc(t(i-1));
    
  % Update subunits
  % Noise terms are non-zero for Subunit Noise model
  m = m0 + dt*(alpham(V0)*(1-m0) - betam(V0)*m0) + mNoise(V0,m0,i)*sqrt(dt);  % shifted to i-1 in function
  h = h0 + dt*(alphah(V0)*(1-h0) - betah(V0)*h0) + hNoise(V0,h0,i)*sqrt(dt);
  n = n0 + dt*(alphan(V0)*(1-n0) - betan(V0)*n0) + nNoise(V0,n0,i)*sqrt(dt);
  
  % skip that method (from Fox and Lu 1994    )
  while (m-1)*(m-0)>0
      m = m0 + dt*(alpham(V0)*(1-m0) - betam(V0)*m0) + sqrt(2)*sqrt((alpham(V)*(1-m) + betam(V)*m)/NNa) * randn*sqrt(dt);     
  end
  
  while (h-1)*(h-0)>0
     h = h0 + dt*(alphah(V0)*(1-h0) - betah(V0)*h0) + sqrt(2)*sqrt((alphah(V)*(1-h) + betah(V)*h)/NNa) * randn*sqrt(dt);     
  end
  
  while (n-1)*(n-0)>0
     n = n0 + dt*(alphan(V0)*(1-n0) - betan(V0)*n0) + sqrt(2)*sqrt((alphan(V)*(1-n) + betan(V)*n)/NK)  * randn*sqrt(dt);   
     disp('n goes out of bounds') 
  end

  % Enforce boundary conditions (only necessary for subunit noise model
%   m = max(0,min(1,m));
%   h = max(0,min(1,h));
%   n = max(0,min(1,n));
 
  % Note: Impose bounds on fractions to avoid <0 or >1 in dV/dt equation, this doesn't directly alter the dynamics of the subunits or channels
  NaFraction =  m0^3*h0;  % Fluctuations are non-zero for Conductance Noise Models
  KFraction = n0^4;

  
  % Update Voltage
  Vrhs = (-gNa*(NaFraction)*(V0 - ENa)-gK*(KFraction)*(V0 - EK) - gL*(V0-EL) + I)/C;
  V = V0 + dt*Vrhs ;   % VNoise is non-zero for Current Noise Model
  
  % Save Outputs  
  Y(i,1) = V;
  Y(i,2) = m;
  Y(i,3) = h;
  Y(i,4) = n;

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

% subunit kinetics (Hodgkin and Huxley parameters)
function out = alpham(V)
  out =  0.1 * (25-V)/ (exp((25-V)/10)-1);
end

function out = betam(V)
  out = 4 * exp(-V/18);
end

function out = alphah(V)
  out =  0.07 * exp(-V/20);
end

function out = betah(V)
  out = 1/ (exp((30-V)/10)+1);
end

function out = alphan(V)
  out = 0.01 * (10-V) / (exp((10-V)/10)-1);
end

function out = betan(V)
  out = 0.125 * exp(-V/80);
end
