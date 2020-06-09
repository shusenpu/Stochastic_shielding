function [EV_Na2,EV_K2,CoeffK,CoeffNa,CoeffNah,ErrorK,ErrorNa,GL] = HHeigens_new(V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Shusen Pu and Peter Thomas
% June, 2020
% find eigen values and eigenvectors for Na and K drift matrix
% model: 14D HH 

%%% Inputs
% V voltage (mV)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Na Current
gNa = 120; % mS/cm^2

% K Current
gK = 36; % mS/cm^2

% Passive Leak
gL = 0.3; % mS / cm^2

nt=length(V);

EV_Na=nan(nt,1);
EV_K=nan(nt,1);

EV_Na2=nan(nt,1);
EV_K2=nan(nt,1);

EV_Na3=nan(nt,1);
EV_K3=nan(nt,1);

EVec_Na=nan(8,nt);
V_dN=nan(5,nt);
V_dm=nan(8,nt);
V_dh=nan(8,nt);
EVec_K=nan(5,nt);

CoeffK=nan(5,nt);
CoeffNa=nan(2,nt);
CoeffNah=nan(8,nt);

ErrorK=nan(nt,1);
ErrorNa=nan(nt,1);
GL=nan(nt,1);

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
        [-4*alphan(V), betan(V)             , 0                , 0                  , 0;
          4*alphan(V), -3*alphan(V)-betan(V),  2*betan(V)               , 0,                   0;
          0,        3*alphan(V),        -2*alphan(V)-2*betan(V), 3*betan(V),          0;
          0,        0,               2*alphan(V),          -alphan(V)-3*betan(V), 4*betan(V);
          0,        0,               0,                 alphan(V),          -4*betan(V)];
    
    

for i=1:nt
 disp(i)
  % Input Current
  V0 = V(i);
  % update the gates

  [V_Na,D_Na] = eig(ANa(V0));
  [V_K,D_K] = eig(AK(V0));
  
%   d_Na=diag(D_Na);
%   D2=sort(d_Na,'descend'); % make diagonal matrix out of sorted diagonal values of input D
  [D2_Na, ind1]=sort(diag(D_Na),'descend'); % store the indices of which columns the sorted eigenvalues come from
  P_Na=V_Na(:,ind1); % arrange the columns in this order

  [D2_K, ind2]=sort(diag(D_K),'descend'); % store the indices of which columns the sorted eigenvalues come from
  P_K=V_K(:,ind2); % arrange the columns in this order
  
  EV_Na(i)=D2_Na(2);
  EV_K(i)=D2_K(2);
  
   EV_Na2(i)=D2_Na(3);
  EV_K2(i)=D2_K(3);
  
     EV_Na3(i)=D2_Na(4);
  EV_K3(i)=D2_K(4);
  
  EVec_Na(:,i)=P_Na(:,2);
  EVec_K(:,i)=P_K(:,2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m = alpham(V0) / (alpham(V0) + betam(V0)); % m
  h = alphah(V0) / (alphah(V0) + betah(V0)); % h
  n = alphan(V0) / (alphan(V0) + betan(V0)); % n
  
  V_dN(:,i)=[-4*(1-n)^3;4*(1-n)^2*(1-4*n);12*n*(1-n)*(1-2*n);4*n^2*(3-4*n);4*n^3];
  
  %----check whether dn is parall to the second eigenvector----
  %find the coefficients
  CoeffK(:,i)=V_dN(:,i)./P_K(:,2);
  %find the error if they are linearly dependent
  ErrorK(i)=norm(P_K(:,2)*mean(CoeffK(:,i))-V_dN(:,i));
  
  V_dm(:,i)=[-3*(1-m)^2*(1-h);3*(1-h)*(3*m^2-4*m+1);3*(1-h)*(2*m-3*m^2);3*(1-h)*m^2;...
      -3*h*(1-m)^2;3*h*(3*m^2-4*m+1);3*h*(2*m-3*m^2);3*h*m^2];
  
  V_dh(:,i)=[-(1-m)^3;-3*(1-m)^2*m;-3*(1-m)*m^2;-m^3;...
      (1-m)^3;3*(1-m)^2*m;3*(1-m)*m^2;m^3];
  %combine dm and dh into a matrix
  V_temp=[V_dm(:,i),V_dh(:,i)];
  %find the lease square solution
  CoeffNa(:,i)=V_temp\P_Na(:,2);
  CoeffNah(:,i)=V_dh(:,i)./P_Na(:,2);
  
  ErrorNa(i)=norm(P_Na(:,2)*mean(CoeffNah(:,i))-V_dh(:,i));
  GL(i)=-(gL+gNa*m^3*h+gK*n^4);
%   keyboard
 
  
end

end % End Function Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF SOLVER                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define functions used above

% subunit kinetics (Hodgkin and Huxley parameters, modified from XPP code by PJT OCt, 2018)
function out = alpham(V)

if V==-40
      out=1;  
else
    out =  .1*(V+40)/(1-exp(-(V+40)/10));
  %0.1 * (25-V)/ (exp((25-V)/10)-1);
end

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

if V==-55
      out=0.1;  
else
  out = .01*(V+55)/(1-exp(-(V+55)/10));
  %0.01 * (10-V) / (exp((10-V)/10)-1);
end

end

function out = betan(V)
  out = .125*exp(-(V+65)/80);
  %0.125 * exp(-V/80);
end

