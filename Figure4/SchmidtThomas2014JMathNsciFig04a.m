% April 8, 2013
% Computing R_k error terms for the Hodgkin-Huxley potassium channel (5-
% state chain with non-symmetric transition rates) - This version uses
% right and left eigenvectors to compute R_k.

clear all


%% Basic parameters - input number of states and number of reactions below

nstates=5;  % number of states
nrxns=8;    % number of reactions (transitions)
Ntot=500; % number of random walkers

% Define range of voltage values
vol_list = -100:10:100;

%% Definitions, initializations

% Define measurement vector C (e.g. conductance)
C = zeros(1,nstates);
C(nstates) = 1;
% Want C to be a column vector
C = C';

% Initialize a matrix in which to store each R vector for different voltage
% values (rows of the matrix)
R_vol = zeros(length(vol_list),nrxns);
R_abs_vol = zeros(length(vol_list),nrxns);

% Initialize a matrix to keep the leading eigenvector as a function of
% voltage (rows of V_vol)
V_vol = zeros(length(vol_list),nstates);
lambda_vol = zeros(length(vol_list),nstates);

% Initialize parameters alpha and beta, functions of voltage
alpha = zeros(1,length(vol_list));
beta = zeros(1,length(vol_list));

%% Computation of error terms R_k over a range of voltage values

for jlist=1:length(vol_list)   % Loop over all i_vol values
    i_vol = vol_list(jlist);

    % Compute parameters alpha and beta which are functions of voltage
    alpha(jlist) = 0.01*(i_vol+55)/(1-exp(-0.1*(i_vol+55)));
    if(i_vol==-55)
        alpha(jlist)=(0.01*(-54+55)/(1-exp(-0.1*(i_vol+55))) + 0.01*(-56+55)/(1-exp(-0.1*(i_vol+55))))/2;
    end
    
    beta(jlist) = 0.125*exp(-(i_vol+65)/80);
    if(i_vol==-65)
        beta(jlist)=(0.125*exp(-(-64+65)/80) + 0.125*exp(-(-66+65)/80))/2;
    end
    
    % Compute the WEIGHTED adjacency matrix A for this 5-state chain, given a particular
    % voltage value i_vol
    A = zeros(nstates);
    A(1,2)=4*alpha(jlist); 
    A(2,3)=3*alpha(jlist); 
    A(3,4)=2*alpha(jlist); 
    A(4,5)=alpha(jlist);

    A(5,4)=4*beta(jlist);
    A(4,3)=3*beta(jlist);
    A(3,2)=2*beta(jlist);
    A(2,1)=beta(jlist);

    D=diag(sum(A,2));   % degree matrix, row sums of weighted A yield out-degrees
    L=(A-D)'; % use the weighted A to compute the "graph Laplacian" L, want the transpose
    
    %%% First, find the steady-state distribution
    [vv,dd]=eig(L);
    vv_idx=find(diag(dd)==max(diag(dd)));
    prob_ss=vv(:,vv_idx)/sum(vv(:,vv_idx)); % normalize the last vector -- corresponding to the null eigenvalue, we hope!
    N_ave_ss=Ntot*prob_ss; % number of random walkers, on average, at each node
    clear vv dd prob_ss 
    
    % Compute voltage dependent "noise" matrix B (sigma_k depends on 
    % voltage for each reaction k)
    
    B = zeros(nstates,nrxns);
    B(1,1)=-sqrt(4*alpha(jlist)*N_ave_ss(1));   B(1,2)=sqrt(beta(jlist)*N_ave_ss(2));  
    B(2,1)=sqrt(4*alpha(jlist)*N_ave_ss(1));    B(2,2)=-sqrt(beta(jlist)*N_ave_ss(2));
    
    B(2,3)=-sqrt(3*alpha(jlist)*N_ave_ss(2));   B(2,4)=sqrt(2*beta(jlist)*N_ave_ss(3));
    B(3,3)=sqrt(3*alpha(jlist)*N_ave_ss(2));    B(3,4)=-sqrt(2*beta(jlist)*N_ave_ss(3));
    
    B(3,5)=-sqrt(2*alpha(jlist)*N_ave_ss(3));   B(3,6)=sqrt(3*beta(jlist)*N_ave_ss(4));    
    B(4,5)=sqrt(2*alpha(jlist)*N_ave_ss(3));    B(4,6)=-sqrt(3*beta(jlist)*N_ave_ss(4));
    
    B(4,7)=-sqrt(alpha(jlist)*N_ave_ss(4));     B(4,8)=sqrt(4*beta(jlist)*N_ave_ss(5));    
    B(5,7)=sqrt(alpha(jlist)*N_ave_ss(4));      B(5,8)=-sqrt(4*beta(jlist)*N_ave_ss(5));    
              

    %% Computation of eigenvalues and eigenvectors

    n = nstates;    % Number of states in the system
    m = nrxns;      % Number of stochastic fluctuations (reactions) to sum over 

    % RIGHT eigenvectors
    [V,vd] = eig(L);

    lambda=0;
    % Eigenvalues 
    lambda = diag(vd);

    % Sort the eigenvalues: 0 is first, then list in descending order
    [lambda,index] = sort(lambda,'descend'); 
    lambda_vol(jlist,:) = lambda(:);

    % Sort the eigenvectors to match the ordering of the eigenvalues
    V = V(:,index);

    % Normalize the leading eigenvector (given by V(:,index(1))) and save
    % it in the jlist row of matrix V_vol
    Vnorm = zeros(nstates,1);   % Initialize
    Vnorm(:) = V(:,1)/sum(V(:,1));
    V_vol(jlist,:) = Vnorm(:);
       
    % LEFT eigenvectors
    [U,ud] = eig(L');

    lambdaU=0;
    % Eigenvalues 
    lambdaU = diag(ud);
    
    % Sort the eigenvalues: 0 is first, then list in descending order
    [lambdaU,indexU] = sort(lambdaU,'descend');

    % Sort the eigenvectors to match the ordering of the eigenvalues
    U = U(:,indexU);
    
    % Normalize U and call it W so that W'V = delta_ij
    for j=1:n
       W(:,j) = U(:,j)/(U(:,j)'*V(:,j)); 
    end
    
           
    
    %% Computation of errors R_k
    
    % Initialize vector for error terms R = [R_1 ... R_k ... R_m], for each
    % voltage value
    R = 0;
    R_abs = 0;
    
    % Sum over each stochastic fluctuation (given by the B_k's) 
    for k=1:m

        % Define projection matrix P_k so that B_k = B * P_k and B = \sum_k B_k
        P_k = zeros(m);
        P_k(k,k) = 1;

        % Compute B_k for this value of k
        B_k = B*P_k;

        % Compute R_k for each edge k
        S = 0;
        Sabs = 0;
        for j=2:n  % exclude the i=j=1 term since lambda(1)=0
            for i=1:n
                S = S + (1/(lambda(i)+lambda(j)))*C'*V(:,i)*W(:,i)'*B_k*B_k'*W(:,j)*V(:,j)'*C;

                % Also sum over absolute value of each eigenpair of Q - compare with S above
                Sabs = Sabs + abs(1/(lambda(i)+lambda(j))*C'*V(:,i)*W(:,i)'*B_k*B_k'*W(:,j)*V(:,j)'*C);
            end
        end
        R(k) = S;
        R_abs(k) = Sabs;
    end
    R;
    R_abs;
    
    R_vol(jlist,:) = R;
    R_abs_vol(jlist,:) = R_abs;    
       
    
end   % end of for loop for vol_list



%% Plot of each R_k as a function of voltage
figure
h=plot(-R_vol, '*-');
set(gca,'FontSize',14)
str1 = sprintf('HH K channel: R_k as a function of voltage');
title(str1);
xlabel('Voltage')
ylabel('Importance R_k');
xlim([0,length(vol_list)+1]);  %ylim([-10,300]);
set(gca,'xTick',1:2:length(vol_list))
set(gca,'xTickLabel',{'-100','-80','-60','-40','-20','0','20','40','60','80','100'})
legend([h(2),h(4),h(6),h(8)],'R_1,R_2','R_3,R_4','R_5,R_6','R_7,R_8')


