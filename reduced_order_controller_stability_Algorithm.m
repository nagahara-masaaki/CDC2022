%% MATLAB program for the numerical example in
%% Nagahara, Iwata, and Sebe, Algorithms 2022

%% Initialization
clear;
ops=sdpsettings('solver','sdpt3');

%% Syste matrices
% COMPleib benchmark
system_model = 'NN1';

[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib(system_model);

B=B2; C=C2;
D=zeros(size(C2,1),size(B2,2));

%%%%%%%
Ps = ss(A,B,C,D);
[n,m] = size(B);
rb = rank(B);
rc = rank(C);

%% Order of controller
nc = 0;

%% Reduced order controller by minimizing the nuclear norm
X = sdpvar(n,n);
Y = sdpvar(n,n);

Bp = null(B')'; % B "perp"
Cp = null(C)'; % C' "perp"

In = eye(n);

epsil1 = 1e-8;
epsil2 = 1e-8;

LMI1 = [X,In;In,Y] - epsil1*zeros(2*n,2*n);
LMI2 = -Bp*[In,A]*[zeros(n,n),X;X,epsil2*X]*[In;A']*Bp'  - epsil1*eye(n-rb);
LMI3 = Cp*[A', -In]*[-epsil2*Y, Y;Y, zeros(n,n)]*[A;-In]*Cp' - epsil1*eye(n-rc);

LMI = [LMI1>=0,LMI2>=0,LMI3>=0];

%% Nuclear norm minimization
diagnostic = optimize(LMI,norm([X,In;In,Y],'nuclear'),ops);

X_L1 = value(X);
Y_L1 = value(Y);

%% Initialization for projection
% initial guess for (X,Y)
X00 = zeros(n,n);
Y00 = zeros(n,n);
% maximum number of iterations
max_iter = 1e2;
% stopping criterion
EPS = 1e-6;

%% Alternating Projetion by Grigoriadis and Skelton '96 (GS96)
% parameters
X0 = X00; Y0 = Y00;
Px = X0; Py = Y0;
Qx = X0; Qy = Y0;

res = [];
res2= [];
for k = 1:max_iter
    fprintf('(GS96) k=%d\n',k);
    
    % projections for rank condition
    T = [X0,In;In,Y0];
    [U0,S0,V0] = svd(T);
    diag_S1 = s_sparse_operator(diag(S0),n+nc);
    S1 = diag(diag_S1);
    X1 = U0(1:n,:)*S1*V0(1:n,:)';
    Y1 = U0(n+1:2*n,:)*S1*V0(n+1:2*n,:)';
    
    X1_GS=X1; Y1_GS=Y1;
    
    res = [res,norm(X1*Y1-eye(n),'fro')];
    
    % feasibility check
    assign(X,X1,1); assign(Y,Y1,1);
    g1 = min(eig(double(LMI1)))+epsil1;
    g2 = min(eig(double(LMI2)))+epsil1;
    g3 = min(eig(double(LMI3)))+epsil1;
    
    if g1>0 & g2>0 & g3>0
        fprintf('normal exit from iteration (after projection rank)\n');
        break;
    end
    
    % projection for LMIs
    [X0,Y0] = Proj_LMI(X1,Y1,X,Y,LMI);
    
    X1_GS=X0; Y1_GS=Y0;
    
    % feasibility check
    assign(X,X0); assign(Y,Y0);
    g1 = min(eig(double(LMI1)))+epsil1;
    g2 = min(eig(double(LMI2)))+epsil1;
    g3 = min(eig(double(LMI3)))+epsil1;
    
    res2tmp=norm(X0*Y0-eye(n),'fro');
    res2=[res2,res2tmp];
    
    % rank check
    if res2tmp<EPS
        fprintf('normal exit from iteration (after projection LMI)\n');
        break;
    end
    
end

res_GS = res;
res2_GS = res2;


%% Alternating Projection (proposed)
X0 = X00; Y0 = Y00;
Px = X0; Py = Y0;
Qx = X0; Qy = Y0;

res = [];
res2= [];
for k = 1:max_iter
    fprintf('(proposed) k=%d\n',k)
    
    % projection for rank condition
    [X1,Y1] = Proj_rank(X0,Y0,n+nc);
    res = [res,norm(X1*Y1-eye(n),'fro')];
    
    X1_P=X1; Y1_P=Y1;
    
    % feasibility check
    assign(X,X1); assign(Y,Y1);
    g1 = min(eig(double(LMI1)))+epsil1;
    g2 = min(eig(double(LMI2)))+epsil1;
    g3 = min(eig(double(LMI3)))+epsil1;
    
    if g1>0 & g2>0 & g3>0
        fprintf('normal exit from projection rank\n');
        break;
    end
    
    % projection for LMIs
    [X0,Y0] = Proj_LMI(X1,Y1,X,Y,LMI);
    
    X1_P=X0; Y1_P=Y0;
    
    % feasibility check
    assign(X,X0); assign(Y,Y0);
    g1 = min(eig(double(LMI1)))+epsil1;
    g2 = min(eig(double(LMI2)))+epsil1;
    g3 = min(eig(double(LMI3)))+epsil1;
    
    res2tmp=norm(X0*Y0-eye(n),'fro');
    res2=[res2,res2tmp];
    
    % rank check
    if res2tmp<EPS
        fprintf('normal exit from iteration (after projection LMI)\n');
        break;
    end
    
end

res_P = res;
res2_P = res2;


%% Analysis
% ||XY-I||
figure;
semilogy(1:max_iter,res,'-');
hold on
semilogy(1:max_iter,res_GS,'-.');
xlim([1,max_iter])
xlabel('iteration number k')
title("||XY-I||")
grid on;

figure;
semilogy(1:max_iter,res2,'-');
hold on
semilogy(1:max_iter,res2_GS,'-.');
xlim([1,max_iter])
xlabel('iteration number k')
title("||XY-I||")
grid on;


% Controller
% nuclear norm
infeas = sum(sum(isnan(X_L1))) + sum(sum(isnan(Y_L1)));
if infeas > 0
    disp(sprintf('Residual (NNM): NaN'))
    disp('feedback system: unstable')
else
    [K1,Pcl1] = controller_realization(Ps,X_L1,Y_L1,nc);
    disp(sprintf('Residual (NNM): %f',norm(X_L1*Y_L1-eye(n),'fro')))
    if max(real(pole(Pcl1))) < 0
        disp('feedback system: stable')
    else
        disp('feedback system: unstable')
    end
end

% GS96
[K2,Pcl2] = controller_realization(Ps,X1_GS,Y1_GS,nc);
disp(sprintf('Residual (GS96): %f', norm(X1_GS*Y1_GS-eye(n),'fro')))
if max(real(pole(Pcl2))) < 0
    disp('feedback system: stable')
else
    disp('feedback system: unstable')
end

% Proposed method
[K3,Pcl3] = controller_realization(Ps,X1_P,Y1_P,nc);
disp(sprintf('Residual (proposed): %f', norm(X1_P*Y1_P-eye(n),'fro')))
if max(real(pole(Pcl3))) < 0
    disp('feedback system: stable')
else
    disp('feedback system: unstable')
end



return;
