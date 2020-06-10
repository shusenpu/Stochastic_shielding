function [Y,T_run] = MC(t, Ifunc, Area)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Shusen Pu
% Aug 08, 2019

%%% Inputs
% t is vector of time values (ms)
% Ifunc is a function @(t)f(t) that returns stimulus value as a function of time (in ms)
% SigmaIn: st. dev. of current noise
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
t_ind=tic;
% time step size
dt = t(2)-t(1);
% Number of time steps
nt = length(t);  % total
nt1 = nt-1;  % at which to solve

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
 
% Initialize Output
Y = zeros(nt,15); 
Y(1,:) = [t0, V0, Na_gates',K_gates'];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parameter Values
% 
% % Number of Channels
% NNa = round(Area*60); % Na
% NK = round(Area*18); % K
% 
% % Capacitance
% C = 1; % muF /cm^2
% 
% % Na Current
% gNa = 120; % mS/cm^2
% ENa = 120; % mV
% 
% % K Current
% gK = 36; % mS/cm^2
% EK = -12; % mV
% 
% % Passive Leak
% gL = 0.3; % mS / cm^2
% EL = 10.6; % mV


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


%states occupancy

MCNa=zeros(1,8);
MCNa(1) = floor(NNa*m00);
MCNa(2) = floor(NNa*m01); 
MCNa(3) = floor(NNa*m02);
MCNa(4) = floor(NNa*m03); 
MCNa(5) = floor(NNa*m10);
MCNa(6) = floor(NNa*m11);
MCNa(7) = floor(NNa*m12); 
MCNa(8) = NNa - sum(sum(MCNa));

MCK=zeros(1,5);
MCK(1:4) = floor(NK*[n0 n1 n2 n3]);
MCK(5) = NK-sum(sum(MCK));

%save MCNA(4,2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HERE IS THE SOLVER            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GILLESPIE FOR MARKOV CHAIN    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:nt

  % Input Current
  I = Ifunc(t(i-1));
    
  
  % Compute Fraction of open channels

    [MCNa, MCK]= MarkovChainFraction(V0, MCNa, MCK, t0,dt);
    NaFraction = MCNa / NNa;
    KFraction = MCK / NK;
 
  
  % Update Voltage
  Vrhs = (-gNa*(NaFraction(8))*(V0 - ENa)-gK*(KFraction(5))*(V0 - EK) - gL*(V0-EL) + I)/C;
  V = V0 + dt*Vrhs;   % VNoise is non-zero for Current Noise Model
  
  % Save Outputs  
  Y(i,1) = t(i);
  Y(i,2) = V;
  Y(i,3:10) = NaFraction;
  Y(i,11:15) = KFraction;

  % Keep "old values" to use in next Euler time step
  V0 = V;

  
end  % End loop over time for SDE solver

T_run=toc(t_ind);
end % End Function Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF SOLVER                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define functions used above

% % subunit kinetics (Hodgkin and Huxley parameters)
% function out = alpham(V)
%   out =  0.1 * (25-V)/ (exp((25-V)/10)-1);
% end
% 
% function out = betam(V)
%   out = 4 * exp(-V/18);
% end
% 
% function out = alphah(V)
%   out =  0.07 * exp(-V/20);
% end
% 
% function out = betah(V)
%   out = 1/ (exp((30-V)/10)+1);
% end
% 
% function out = alphan(V)
%   out = 0.01 * (10-V) / (exp((10-V)/10)-1);
% end
% 
% function out = betan(V)
%   out = 0.125 * exp(-V/80);
% end

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

% Markov chain
function [NaStateOut, KStateOut]= MarkovChainFraction(V, NaStateIn, KStateIn, t,dt)

tswitch = t;
Nastate = NaStateIn;
Kstate = KStateIn;
  % Update Channel States
  while (tswitch < (t+dt)) 

% MCNa(1) = floor(NNa*m00);
% MCNa(2) = floor(NNa*m01); 
% MCNa(3) = floor(NNa*m02);
% MCNa(4) = floor(NNa*m03); 
% MCNa(5) = floor(NNa*m10);
% MCNa(6) = floor(NNa*m11);
% MCNa(7) = floor(NNa*m12); 
% MCNa(8) = NNa - sum(sum(MCNa));

%  1 MCNa(1,1) = floor(NNa*(1-m0)^3*(1-h0));
%  2  MCNa(2,1) = floor(NNa*3*(1-m0)^2*m0*(1-h0)); 
%  3 MCNa(3,1) = floor(NNa*3*(1-m0)^1*m0^2*(1-h0));
%  4 MCNa(4,1) = floor(NNa*(1-m0)*m0^3*(1-h0)); 
%  5 MCNa(1,2) = floor(NNa*(1-m0)^3*(h0));
%  6 MCNa(2,2) = floor(NNa*3*(1-m0)^2*m0*(h0));
%  7 MCNa(3,2) = floor(NNa*3*(1-m0)^1*m0^2*(h0)); 
%  8 MCNa(4,2) = NNa - sum(sum(MCNa));

    % Determine which state switches by partitioning total rate into its 28 components
    rate(1) = 3.*alpham(V) * Nastate(1);
    rate(2) = rate(1) + 2.*alpham(V) * Nastate(2);
    rate(3) = rate(2) + 1.*alpham(V) * Nastate(3);
    rate(4) = rate(3) + 3.*betam(V) * Nastate(4);
    rate(5) = rate(4) + 2.*betam(V) * Nastate(3);
    rate(6) = rate(5) + 1.*betam(V) * Nastate(2);
    rate(7) = rate(6) + alphah(V) * Nastate(1);
    rate(8) = rate(7) + alphah(V) * Nastate(2);
    rate(9) = rate(8) + alphah(V) * Nastate(3);
    rate(10) = rate(9) + alphah(V) * Nastate(4);
    rate(11) = rate(10) + betah(V) * Nastate(5);
    rate(12) = rate(11) + betah(V) * Nastate(6);
    rate(13) = rate(12) + betah(V) * Nastate(7);
    rate(14) = rate(13) + betah(V) * Nastate(8);
    rate(15) = rate(14) + 3.*alpham(V) * Nastate(5);
    rate(16) = rate(15) + 2.*alpham(V) * Nastate(6);
    rate(17) = rate(16) + 1.*alpham(V) * Nastate(7);
    rate(18) = rate(17) + 3.*betam(V) * Nastate(8);
    rate(19) = rate(18) + 2.*betam(V) * Nastate(7);
    rate(20) = rate(19) + 1.*betam(V) * Nastate(6);
    rate(21) = rate(20) + 4.*alphan(V) * Kstate(1);
    rate(22) = rate(21) + 3.*alphan(V) * Kstate(2);
    rate(23) = rate(22) + 2.*alphan(V) * Kstate(3);
    rate(24) = rate(23) + 1.*alphan(V) * Kstate(4);
    rate(25) = rate(24) + 4.*betan(V) * Kstate(5);
    rate(26) = rate(25) + 3.*betan(V) * Kstate(4);
    rate(27) = rate(26) + 2.*betan(V) * Kstate(3);
    rate(28) = rate(27) + 1.*betan(V) * Kstate(2);

    % Total Transition Rate
    totalrate = rate(28);

    % Exponential Waiting Time Distribution
    tupdate = -log(rand()) / totalrate;

    % Time of Next Switching Event (Exp Rand Var)
    tswitch = tswitch + tupdate;

    if (tswitch < (t+dt)) 

      % Scaled Uniform RV to determine which state to switch
      r = totalrate*rand();

      if (r < rate(1)) 
        Nastate(1) = Nastate(1)-1;
        Nastate(2) = Nastate(2)+1 ;
      elseif (r < rate(2)) 
       Nastate(2) = Nastate(2)-1;
       Nastate(3) = Nastate(3)+1 ;
      elseif (r < rate(3)) 
       Nastate(3) = Nastate(3)-1;
       Nastate(4) = Nastate(4)+1 ;
      elseif (r < rate(4)) 
       Nastate(4) = Nastate(4)-1;
       Nastate(3) = Nastate(3)+1 ; 
      elseif (r < rate(5)) 
       Nastate(3) = Nastate(3)-1;
       Nastate(2) = Nastate(2)+1;
      elseif (r < rate(6)) 
       Nastate(2) = Nastate(2)-1;
       Nastate(1) = Nastate(1)+1;
      elseif (r < rate(7)) 
       Nastate(1) = Nastate(1)-1;
       Nastate(5) = Nastate(5)+1;
      elseif (r < rate(8)) 
       Nastate(2) = Nastate(2)-1;
       Nastate(6) = Nastate(6)+1;
      elseif (r < rate(9)) 
       Nastate(3) = Nastate(3)-1;
       Nastate(7) = Nastate(7)+1;
      elseif (r < rate(10)) 
       Nastate(4) = Nastate(4)-1;
       Nastate(8) = Nastate(8)+1;
      elseif (r < rate(11)) 
       Nastate(5) = Nastate(5)-1;
       Nastate(1) = Nastate(1)+1;
      elseif (r < rate(12)) 
       Nastate(6) = Nastate(6)-1;
       Nastate(2) = Nastate(2)+1;
      elseif (r < rate(13)) 
       Nastate(7) = Nastate(7)-1;
       Nastate(3) = Nastate(3)+1;
      elseif (r < rate(14)) 
       Nastate(8) = Nastate(8)-1;
       Nastate(4) = Nastate(4)+1;
      elseif (r < rate(15)) 
       Nastate(5) = Nastate(5)-1;
       Nastate(6) = Nastate(6)+1;
      elseif (r < rate(16)) 
       Nastate(6) = Nastate(6)-1;
       Nastate(7) = Nastate(7)+1;
      elseif (r < rate(17)) 
       Nastate(7) = Nastate(7)-1;
       Nastate(8) = Nastate(8)+1;
      elseif (r < rate(18)) 
       Nastate(8) = Nastate(8)-1;
       Nastate(7) = Nastate(7)+1;
      elseif (r < rate(19)) 
       Nastate(7) = Nastate(7)-1;
       Nastate(6) = Nastate(6)+1;
      elseif (r < rate(20)) 
       Nastate(6) = Nastate(6)-1;
       Nastate(5) = Nastate(5)+1;
      elseif (r < rate(21)) 
       Kstate(1) = Kstate(1)-1;
       Kstate(2) = Kstate(2)+1;
      elseif (r < rate(22)) 
       Kstate(2) = Kstate(2)-1;
       Kstate(3) = Kstate(3)+1;
      elseif (r < rate(23)) 
       Kstate(3) = Kstate(3)-1;
       Kstate(4) = Kstate(4)+1;
      elseif (r < rate(24)) 
       Kstate(4) = Kstate(4)-1;
       Kstate(5) = Kstate(5)+1;
      elseif (r < rate(25)) 
       Kstate(5) = Kstate(5)-1;
       Kstate(4) = Kstate(4)+1;
      elseif (r < rate(26)) 
       Kstate(4) = Kstate(4)-1;
       Kstate(3) = Kstate(3)+1;
      elseif (r < rate(27)) 
       Kstate(3) = Kstate(3)-1;
       Kstate(2) = Kstate(2)+1;
      else
       Kstate(2) = Kstate(2)-1;
       Kstate(1) = Kstate(1)+1;
      end % End if statement

    end % end if tswitch<dt
 
end % end while tswitch<dt
  
NaStateOut = Nastate;
KStateOut = Kstate;
end % end Markov chain Gillespie update function