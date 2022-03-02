%% MATLAB program for the numerical example in 
%% Nagahara and Sebe, CDC 2022

%% Initialization
clear;

%% Syste matrices
% Maciejowski, MPC, p.46
A = [-1.93, 0, 0, 0;
    0.394, -0.426, 0, 0;
    0, 0, -0.63, 0;
    0.82, -0.784, 0.413, -0.426];
B = [1.274, 1.274;
    0, 0;
    1.34, -0.65;
    0, 0];
C = [0,1,0,0;
    0,0,1,0;
    0,0,0,1];
[n,m] = size(B);
[p,n] = size(C);
D = zeros(p,m);

Ps = ss(A,B,C,D);
[n,m] = size(B);
r = rank(B);

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
LMI2 = -Bp*[In,A]*[zeros(n,n),X;X,epsil2*X]*[In;A']*Bp'  - epsil1*eye(n-r);
LMI3 = Cp*[A', -In]*[-epsil2*Y, Y;Y, zeros(n,n)]*[A;-In]*Cp' - epsil1*eye(n-r);

LMI = [LMI1>=0,LMI2>=0,LMI3>=0];

%% Nuclear norm minimization
diagnostic = optimize(LMI,norm([X,In;In,Y],'nuclear'));

X_L1 = value(X);
Y_L1 = value(Y);

%% Initialization for projection
% initial guess for (X,Y)
X00 = zeros(n,n);
Y00 = zeros(n,n);
% maximum number of iterations
max_iter = 15;
% stopping criterion
EPS = 1e-8;

%% Alternating Projetion by Grigoriadis and Skelton '96 (GS96)
% parameters
X0 = X00; Y0 = Y00;
Px = X0; Py = Y0;
Qx = X0; Qy = Y0;

res = [];
for k = 1:max_iter
    % projections for rank condition
    T = [X0,In;In,Y0];
    [U0,S0,V0] = svd(T);
    diag_S1 = s_sparse_operator(diag(S0),n+nc);
    S1 = diag(diag_S1);
    X1 = U0(1:n,:)*S1*V0(1:n,:)';
    Y1 = U0(n+1:2*n,:)*S1*V0(n+1:2*n,:)';
    
    Px = X0 + Px - X1;
    Py = Y0 + Py - Y1;

    res = [res,norm(X1*Y1-eye(n),'fro')];
    
    % feasibility check
    assign(X,X1); assign(Y,Y1);
    g1 = min(eig(double(LMI1)));
    g2 = min(eig(double(LMI2)));
    g3 = min(eig(double(LMI3)));
    if g1>0 & g2>0 & g3>0 & k>1 & abs(res(end)-res(end-1)) < EPS
        break;
    end
    
    % projection for LMIs
    [X0,Y0] = Proj_LMI(X1+Qx,Y1+Qy,X,Y,LMI);
    Qx = X1 + Qx - X0;
    Qy = Y1 + Qy - Y0;
    fprintf('(GS96) k=%d\n',k)
end

X1_GS = X1;
Y1_GS = Y1;
res_GS = res;


%% Iterative greedy LMI (proposed)
X0 = X00; Y0 = Y00;
Px = X0; Py = Y0;
Qx = X0; Qy = Y0;

res = [];
for k = 1:max_iter
    % projection for rank condition
    [X1,Y1] = Proj_rank(X0,Y0,n+nc);
    Px = X0 + Px - X1;
    Py = Y0 + Py - Y1;
    
    res = [res,norm(X1*Y1-eye(n),'fro')];
    
    % feasibility check
    assign(X,X1); assign(Y,Y1);
    g1 = min(eig(double(LMI1)));
    g2 = min(eig(double(LMI2)));
    g3 = min(eig(double(LMI3)));
    if g1>0 & g2>0 & g3>0 & k>1 & abs(res(end)-res(end-1)) < EPS
        break;
    end
    
    % projection for LMIs
    [X0,Y0] = Proj_LMI(X1+Qx,Y1+Qy,X,Y,LMI);
    Qx = X1 + Qx - X0;
    Qy = Y1 + Qy - Y0;
    fprintf('(proposed) k=%d\n',k)
end

X1_P = X1;
Y1_P = Y1;
res_P = res;

%% Analysis
% ||XY-I||
figure;
semilogy(1:max_iter,res,'o-');
hold on
semilogy(1:max_iter,res_GS,'x--');
title("||XY-I||")
xlim([1,max_iter])
xlabel('iteration number k')

% Controller
% nuclear norm
[K1,Pcl1] = controller_realization(Ps,X_L1,Y_L1,nc);
disp('Controller (nuclear norm minimization)')
K1
disp('Poles (nuclear norm minimization)')
pole(Pcl1)

% GS96
[K2,Pcl2] = controller_realization(Ps,X1_GS,Y1_GS,nc);
disp('Controller (GS96)')
K2
disp('Poles (GS96)')
pole(Pcl2)

% Proposed method
[K3,Pcl3] = controller_realization(Ps,X1_P,Y1_P,nc);
disp('Controller (proposed)')
K3
disp('Poles (proposed)')
pole(Pcl3)

return;