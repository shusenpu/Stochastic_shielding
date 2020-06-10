% January 21, 2012
% Computing R_k error terms for the Hodgkin-Huxley sodium channel (2 tiered 4-
% state chain with non-symmetric transition rates): 8 states, 20 reactions.
% NOTE: This version uses right and left eigenvectors to compute R_k.
clear all

%% Basic parameters - input number of states and number of reactions below

nstates=8;  % number of states
nrxns=20;    % number of reactions (transitions)
Ntot=500; % number of random walkers

% Define range of voltage values 
vol_list = -100:10:100;


%% Definitions, initializations

% Define measurement vector C (e.g. conductance)
C = zeros(1,nstates);
C(nstates) = 1;
% C(6) = 1;
% C(2) = 1;
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
alpha_m = zeros(1,length(vol_list));
beta_m = zeros(1,length(vol_list));
alpha_h = zeros(1,length(vol_list));
beta_h = zeros(1,length(vol_list));


%% Computation of error terms R_k over a range of voltage values

for jlist=1:length(vol_list)   % Loop over all i_vol values
    i_vol = vol_list(jlist);

    % Compute parameters alpha_m and beta_m which are functions of voltage,
    % corresponding to horizontal transitions between each 4-state chain
    alpha_m(jlist) = 0.1*(i_vol+40)/(1-exp(-(i_vol+40)/10));
    if(i_vol==-40)
        alpha_m(jlist) = (0.1*(-39+40)/(1-exp(-(-39+40)/10)) + 0.1*(-41+40)/(1-exp(-(-41+40)/10)))/2;
    end
    beta_m(jlist) = 4*exp(-(i_vol+65)/18);
    
    % Compute parameters alpha_h and beta_h which are functions of voltage,
    % corresponding to vertical transitions between each 4-state chain
    alpha_h(jlist) = 0.07*exp(-(i_vol+65)/20);
    beta_h(jlist) = 1/(1+exp(-(i_vol+35)/10));
    if(i_vol==-35)
        beta_h(jlist) = (1/(1+exp(-(-34+35)/10)) + 1/(1+exp(-(-36+35)/10)))/2;
    end
    
    % Compute the WEIGHTED adjacency matrix A for this 5-state chain, given a particular
    % voltage value i_vol
    A = zeros(nstates);
    
    A(2,1)=beta_m(jlist);   A(1,2)=3*alpha_m(jlist);
    A(3,2)=2*beta_m(jlist); A(2,3)=2*alpha_m(jlist);
    A(4,3)=3*beta_m(jlist); A(3,4)=alpha_m(jlist);
    A(6,5)=beta_m(jlist);   A(5,6)=3*alpha_m(jlist);
    A(7,6)=2*beta_m(jlist); A(6,7)=2*alpha_m(jlist);
    A(8,7)=3*beta_m(jlist); A(7,8)=alpha_m(jlist);
    A(1,5)=alpha_h(jlist);  A(5,1)=beta_h(jlist);
    A(2,6)=alpha_h(jlist);  A(6,2)=beta_h(jlist);
    A(3,7)=alpha_h(jlist);  A(7,3)=beta_h(jlist);
    A(4,8)=alpha_h(jlist);  A(8,4)=beta_h(jlist);

    D=diag(sum(A,2));   % degree matrix, row sums of weighted A yield out-degrees
    L=(A-D)'; % use the weighted A to compute the "graph Laplacian" L, want the transpose
    
    % Need to weight the noise matrix by the population. First, find the steady-state distribution
    [vv,dd]=eig(L);
    vv_idx=find(diag(dd)==max(diag(dd)));
    prob_ss=vv(:,vv_idx)/sum(vv(:,vv_idx)); % normalize the last vector -- corresponding to the null eigenvalue, we hope!
    N_ave_ss=Ntot*prob_ss; % number of random walkers, on average, at each node
    clear vv dd prob_ss 
    
    % Compute voltage dependent matrix B (sigma_k depends on voltage for
    % each reaction k) - add factors of N_ave_ss(i) to each i->j reaction
    B = zeros(nstates,nrxns);
    
    B(1,1)=-sqrt(3*alpha_m(jlist)*N_ave_ss(1));     B(1,2)=sqrt(beta_m(jlist)*N_ave_ss(2));  
    B(2,1)=sqrt(3*alpha_m(jlist)*N_ave_ss(1));      B(2,2)=-sqrt(beta_m(jlist)*N_ave_ss(2));
    
    B(2,3)=-sqrt(2*alpha_m(jlist)*N_ave_ss(2));     B(2,4)=sqrt(2*beta_m(jlist)*N_ave_ss(3));
    B(3,3)=sqrt(2*alpha_m(jlist)*N_ave_ss(2));      B(3,4)=-sqrt(2*beta_m(jlist)*N_ave_ss(3));
    
    B(3,5)=-sqrt(alpha_m(jlist)*N_ave_ss(3));       B(3,6)=sqrt(3*beta_m(jlist)*N_ave_ss(4));    
    B(4,5)=sqrt(alpha_m(jlist)*N_ave_ss(3));        B(4,6)=-sqrt(3*beta_m(jlist)*N_ave_ss(4));
    
    B(5,7)=-sqrt(3*alpha_m(jlist)*N_ave_ss(5));     B(5,8)=sqrt(beta_m(jlist)*N_ave_ss(6));    
    B(6,7)=sqrt(3*alpha_m(jlist)*N_ave_ss(5));      B(6,8)=-sqrt(beta_m(jlist)*N_ave_ss(6));    

    B(6,9)=-sqrt(2*alpha_m(jlist)*N_ave_ss(6));     B(6,10)=sqrt(2*beta_m(jlist)*N_ave_ss(7));    
    B(7,9)=sqrt(2*alpha_m(jlist)*N_ave_ss(6));      B(7,10)=-sqrt(2*beta_m(jlist)*N_ave_ss(7));
    
    B(7,11)=-sqrt(alpha_m(jlist)*N_ave_ss(7));      B(7,12)=sqrt(3*beta_m(jlist)*N_ave_ss(8));    
    B(8,11)=sqrt(alpha_m(jlist)*N_ave_ss(7));       B(8,12)=-sqrt(3*beta_m(jlist)*N_ave_ss(8));

    B(1,13)=-sqrt(alpha_h(jlist)*N_ave_ss(1));      B(5,14)=-sqrt(beta_h(jlist)*N_ave_ss(5));    
    B(5,13)=sqrt(alpha_h(jlist)*N_ave_ss(1));       B(1,14)=sqrt(beta_h(jlist)*N_ave_ss(5));
    
    B(2,15)=-sqrt(alpha_h(jlist)*N_ave_ss(2));      B(6,16)=-sqrt(beta_h(jlist)*N_ave_ss(6));    
    B(6,15)=sqrt(alpha_h(jlist)*N_ave_ss(2));       B(2,16)=sqrt(beta_h(jlist)*N_ave_ss(6));
    
    B(3,17)=-sqrt(alpha_h(jlist)*N_ave_ss(3));      B(7,18)=-sqrt(beta_h(jlist)*N_ave_ss(7));    
    B(7,17)=sqrt(alpha_h(jlist)*N_ave_ss(3));       B(3,18)=sqrt(beta_h(jlist)*N_ave_ss(7));
    
    B(4,19)=-sqrt(alpha_h(jlist)*N_ave_ss(4));      B(8,20)=-sqrt(beta_h(jlist)*N_ave_ss(8));    
    B(8,19)=sqrt(alpha_h(jlist)*N_ave_ss(4));       B(4,20)=sqrt(beta_h(jlist)*N_ave_ss(8));

    
    %% Computation of eigenvalues and eigenvectors

    n = nstates;    % Number of states in the system
    m = nrxns;      % Number of stochastic fluctuations (B_k's) to sum over 

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
        for j=2:n
            for i=1:n
                S = S + 1/(lambda(i)+lambda(j))*C'*V(:,i)*W(:,i)'*B_k*B_k'*W(:,j)*V(:,j)'*C;

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

%% Plot of each R_k as functions of voltage.
figure
h = plot(-R_vol, '*-');
set(gca,'FontSize',14)
str1 = sprintf('HH Na channel: R_k as a function of voltage');
title(str1);
xlabel('Voltage')
ylabel('Importance R_k');
xlim([0.5,length(vol_list)+0.5]);
set(gca,'xTick',1:2:length(vol_list))
set(gca,'xTickLabel',{'-100','-80','-60','-40','-20','0','20','40','60','80','100'})
legend([h(10),h(12),h(16),h(18),h(20)],'R_9,R_{10}','R_{11},R_{12}','R_{15},R_{16}','R_{17},R_{18}','R_{19},R_{20}')


