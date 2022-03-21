%% MATLAB program for the numerical example in 
%% Nagahara and Sebe, CDC 2022
%% H-inf optimal static controller by bisection search

%% Initialization
clear;
ops=sdpsettings;

%% Syste matrices
% COMPleib benchmark
system_model = 'AC4';
%system_model = 'NN1';
%system_model = 'NN12';
system_model = 'HE6';

[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib(system_model);	% n=4
BD=[B2;D12;zeros(size(B1',1),size(B2,2))];
CD=[C2,D21,zeros(size(C2',2),size(C1,1))];

D=zeros(size(C2,1),size(B2,2));


%%%%%%%
Ps = ss(A,B2,C2,D);
Pg = ss(A,[B1,B2],[C1;C2],[D11,D12;D21,D]);
[n,m] = size(B2);
rb = rank(B2);
rc = rank(C2);
[m1,p1] = size(D11);

%% Order of controller
nc = 0;

%% Bisection search for gamma
if strcmp(system_model, 'AC4')
    gamma_a =  0; %initial guess of lower bound of gamma (AC4)
    gamma_b =   4; %initial guess of upper bound of gamma (AC4)
elseif strcmp(system_model, 'NN1')
    gamma_a =  0; %initial guess of lower bound of gamma (NN1)
    gamma_b = 100; %initial guess of upper bound of gamma (NN1)
elseif strcmp(system_model, 'NN12')
    gamma_a =  0; %initial guess of lower bound of gamma (NN12)
    gamma_b = 100; %initial guess of upper bound of gamma (NN12)
else
    gamma_a =  0; %initial guess of lower bound of gamma
    gamma_b = 1000; %initial guess of upper bound of gamma
end

% 
% 
% 
%

% SDP variables
X = sdpvar(n,n);
Y = sdpvar(n,n);
%gamma = sdpvar(1,1);

% other constant matrices
BDp = null(BD')'; % BD "perp"
CDp = null(CD)'; % CD' "perp"

rb = rank(BD);
rc = rank(CD);

In = eye(n);
Ip = eye(p1);
Im = eye(m1);


% margins, and stopping criteria
epsil1 = 1e-6;
epsil2 = 1e-8;

EPS_res =1e-2;
res = inf;

%% Initial guess for bisection
X1_Psav=zeros(n);
Y1_Psav=zeros(n);

gamall=gamma_b;

norm_check = inf;
Ksav = [];

while(gamma_b-gamma_a > EPS_res)
  gamma = (gamma_a+gamma_b)/2;
    
  %% Reduced order controller
  %% Definition of LMIs
  LMI1 = [X,In;In,Y] - epsil1*zeros(2*n,2*n);
  LMI2 = -BDp*[A*X+X*A', X*C1', B1; C1*X, -gamma*Im, D11; B1',D11', -gamma*Ip]*BDp' - epsil1*eye(n+m1+p1-rb);
  LMI3 = -CDp*[Y*A+A'*Y, Y*B1,  C1';B1'*Y,-gamma*Ip, D11';C1, D11,  -gamma*Im]*CDp' - epsil1*eye(n+m1+p1-rc);

  LMI = [LMI1>=0,LMI2>=0,LMI3>=0];



  %% Initialization for projection
  % initial guess for (X,Y)
  %  X00 = zeros(n,n);
  %  Y00 = zeros(n,n);
  X00 = X1_Psav;
  Y00 = Y1_Psav;

  % maximum number of iterations
  max_iter = 15;
  max_iter = 1e3;
  max_iter = 1e2;
  % stopping criterion
  EPS = 1e-8;

  %% Iterative greedy LMI (proposed)
  X0 = X00; Y0 = Y00;
  %  Px = X0; Py = Y0;
  %  Qx = X0; Qy = Y0;
  %Px = zeros(n); Py = zeros(n);
  %Qx = zeros(n); Qy = zeros(n);

  F_infeasible=0;

  res1 = [];
  for k = 1:max_iter
    % projection for rank condition
    [X1,Y1] = Proj_rank(X0,Y0,n+nc);
    X1_P = X1;
    Y1_P = Y1;

    % feasibility check
    assign(X,X1); assign(Y,Y1);
    g1 = min(eig(double(LMI1)))+epsil1;
    g2 = min(eig(double(LMI2)))+epsil1;
    g3 = min(eig(double(LMI3)))+epsil1;

    %fprintf('after projection onto low rank\n');
    res1tmp=norm(X1*Y1-eye(n),'fro');
    res1 = [res1,res1tmp];
    
    %g1,g2,g3

    if g1>0 & g2>0 & g3>0
      [K3,Pcl3] = controller_realization(Ps,X1_P,Y1_P,nc);
      stab_check = all(real(pole(Pcl3))<0);
      norm_check = norm(lft(Pg,K3),inf,1e-6);
      if stab_check && norm_check<gamma
        gamma=norm_check;
        gammasav=gamma;
        Ksav=K3;
        X1_Psav=X1_P;
        Y1_Psav=Y1_P;
        fprintf('normal exit after projaction rank\n');
        break;
      end
    end
    
    % projection for LMIs
    % fprintf('projection LMI\n');
    [X0,Y0,LMIinfo] = Proj_LMI(X1,Y1,X,Y,LMI);

    if LMIinfo.problem ~= 0
      F_infeasible=1;
      break;
    end
    X1_P = X0;
    Y1_P = Y0;

    assign(X,X0); assign(Y,Y0);
    g1 = min(eig(double(LMI1)))+epsil1;
    g2 = min(eig(double(LMI2)))+epsil1;
    g3 = min(eig(double(LMI3)))+epsil1;


%    fprintf('after projection onto LMIs\n');
    res2tmp=norm(X0*Y0-eye(n),'fro');
%    g1,g2,g3

    [K3,Pcl3] = controller_realization(Ps,X1_P,Y1_P,nc);
    stab_check = all(real(pole(Pcl3))<0);
    norm_check = norm(lft(Pg,K3),inf,1e-6);
    if stab_check && norm_check<gamma
      gamma=norm_check;
      gammasav=gamma;
      Ksav=K3;
      X1_Psav=X1_P;
      Y1_Psav=Y1_P;
      fprintf('normal exit after projaction LMI\n');
      break;
    end

    %    pause

  end

  fprintf('gamma=%f\n', gamma);

  if norm_check>gamma
    % projection loop is terminated by loop counter
    % i.e., no feasible solutions are found
    F_infeasible=1;
  end

  if F_infeasible
    gamma_a=gamma;
  else
    gamma_b=gamma;
  end
  %gamma_a,gamma_b

  %  pause

  gamall=[gamall,gamma_b];
end


Ksav
if isempty(Ksav)
    fprintf('No solution was found');
else
    sysCL=lft(Pg,Ksav);
    max(real(pole(sysCL)))

    fprintf('guaranteed H-infinity norm: %9.6f\n',gammasav);
    fprintf('achieved H-infinity norm:   %9.6f\n',norm(sysCL,inf,1e-6));
end


