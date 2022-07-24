%% MATLAB program for the numerical example in 
%% Nagahara, Iwata, and Sebe, Algorithms 2022
%% H-inf optimal static controller by bisection search

%% Initialization
clear;
ops=sdpsettings;
Fsolver='sdpt3';
%Fsolver='sedumi';
Fsolver='mosek';

%% Syste matrices
% COMPleib benchmark
[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib('AC4');	% n=4
%[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib('NN1');	% n=3
%[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib('NN12');	% n=6
%[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib('HE6');	% n=20
BD=[B2;D12;zeros(size(B1',1),size(B2,2))];
CD=[C2,D21,zeros(size(C2',2),size(C1,1))];

D=zeros(size(C2,1),size(B2,2));


%%%%%%%
Ps = ss(A,B2,C2,D);
[n,m] = size(B2);
rb = rank(B2);
rc = rank(C2);
[m1,p1] = size(D11);

%% Order of controller
nc = 0;

%% Bisection search for gamma
gamma_b = 1000; %initial guess of upper bound of gamma
gamma_a = 0; %initial guess of lower bound of gamma
gamma = (gamma_a+gamma_b)/2;

EPS_res = 0.1;
res = inf;
while(res > EPS_res)
    
%% Reduced order controller by minimizing the nuclear norm
X = sdpvar(n,n);
Y = sdpvar(n,n);
%gamma = sdpvar(1,1);

BDp = null(BD')'; % BD "perp"
CDp = null(CD)'; % CD' "perp"

rb = rank(BD);
rc = rank(CD);

In = eye(n);
Ip = eye(p1);
Im = eye(m1);

epsil1 = 1e-8;
epsil2 = 1e-8;

LMI1 = [X,In;In,Y] - epsil1*zeros(2*n,2*n);
LMI2 = -BDp*[A*X+X*A', X*C1', B1; C1*X, -gamma*Im, D11;B1', D11', -gamma*Ip]*BDp'  - epsil1*eye(n+m1+p1-rb);
LMI3 = -CDp*[Y*A+A'*Y, Y*B1, C1';B1'*Y, -gamma*Ip, D11';C1, D11, -gamma*Im]*CDp' - epsil1*eye(n+m1+p1-rc);

LMI = [LMI1>=0,LMI2>=0,LMI3>=0];


%% Initialization for projection
% initial guess for (X,Y)
X00 = zeros(n,n);
Y00 = zeros(n,n);
% maximum number of iterations
max_iter = 15;
% stopping criterion
EPS = 1e-8;

%% Alternating projection (proposed)
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
    assign(X,X1,1); assign(Y,Y1,1);
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
    %fprintf('(proposed) k=%d\n',k)
end
fprintf('gamma=%f\n', gamma);
X1_P = X1;
Y1_P = Y1;

% Controller
[K3,Pcl3] = controller_realization(Ps,X1_P,Y1_P,nc);
stab_check = all(real(pole(Pcl3))<0);
if stab_check
    gamma_b = gamma;  
else
    gamma_a = gamma;
end
gamma2 = (gamma_a+gamma_b)/2;
res = abs(gamma-gamma2);
gamma = gamma2;
end

return;
