%% MATLAB program for the numerical example in
%% Nagahara, Iwai, and Sebe, Algorithms 2022

%% Initialization
clear;
ops=sdpsettings('solver','sdpt3');

%% Syste matrices
% COMPleib benchmark
system_model = 'NN1';

[A,B1,B2,C1,C2,D11,D12,D21]=COMPleib(system_model);


%% Controller design by hinfstruct
BD=[B2;D12;zeros(size(B1',1),size(B2,2))];
CD=[C2,D21,zeros(size(C2',2),size(C1,1))];

D=zeros(size(C2,1),size(B2,2));

%%%%%%%
P_hinfstruct = ss(A,[B1,B2],[C1;C2],[D11,D12;D21,D]);
n = size(A,1);
m2= size(B2,2);
p2= size(C2,1);

%%%%%%%
% Design by hinfstruct
k0=realp('k0',zeros(m2,p2));
tic;
[kp,gammap]=hinfstruct(P_hinfstruct,k0);
time_hinfstruct = toc;
K_hinfstruct=getValue(kp);
clear BD CD D n m2 p2
%%%%%%%%%%


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

%% CCL (Cone Complementary Linearinzagion)
tic

diagnostic = optimize(LMI,1,ops);
X00 = value(X);
Y00 = value(Y);

% maximum number of iterations
max_iter = 1e2;
% stopping criterion
EPS = 1e-6;

res = [];
res2= [];
X0 = X00;
Y0 = Y00;

for k = 1:max_iter
    fprintf('(CCL) k=%d\n',k);
    diagnostic = optimize(LMI,trace(Y0*X+X0*Y),ops);
    X1 = value(X); Y1 = value(Y);
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
    
    % rank check
    res2tmp=norm(X0*Y0-eye(n),'fro');    
    if res2tmp<EPS
        fprintf('normal exit from iteration (after projection LMI)\n');
        break;
    end
    X0 = X1; Y0 = Y1;
end
time_CCL = toc;
res_CCL = res;
X1_CCL = X1; Y1_CCL = Y1;


%% Nuclear norm minimization
tic;
diagnostic = optimize(LMI,norm([X,In;In,Y],'nuclear'),ops);
time_NNM = toc;

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
tic
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
time_GS = toc;

res_GS = res;
res2_GS = res2;


%% Alternating Projection (proposed)
X0 = X00; Y0 = Y00;
Px = X0; Py = Y0;
Qx = X0; Qy = Y0;

res = [];
res2= [];
tic
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
time_proposed = toc;

res_P = res;
res2_P = res2;


%% Analysis
% ||XY-I||
% figure;
% semilogy(1:max_iter,res,'-');
% hold on
% semilogy(1:max_iter,res_GS,'-.');
% xlim([1,max_iter])
% xlabel('iteration number k')
% title("||XY-I||")
% grid on;
% 
% figure;
% semilogy(1:max_iter,res2,'-');
% hold on
% semilogy(1:max_iter,res2_GS,'-.');
% xlim([1,max_iter])
% xlabel('iteration number k')
% title("||XY-I||")
% grid on;
% 
disp(sprintf('\n'))


% Controller
% hinfstruct
disp(sprintf('**HF(Hinfstruct)**'))
disp(sprintf('time(HF): %f', time_hinfstruct))
disp(sprintf('Residual(HF): not available'))
cp_hinf = pole(lft(P_hinfstruct,K_hinfstruct));
if max(real(cp_hinf))>=0
  disp('feedback system: UNstable')
else
  disp('feedback system: stable')
end
disp(sprintf('\n'))

% CCL
disp(sprintf('**CCL**'))
disp(sprintf('time(CCL): %f', time_CCL))
[K_CCL,Pcl_CCL] = controller_realization(Ps,X1_CCL,Y1_CCL,nc);
disp(sprintf('Residual (CCL): %f', norm(X1_CCL*Y1_CCL-eye(n),'fro')))
if max(real(pole(Pcl_CCL))) < 0
    disp('feedback system: stable')
else
    disp('feedback system: UNstable')
end
disp(sprintf('\n'))

% nuclear norm
disp(sprintf('**NNM**'))
disp(sprintf('time(NNM): %f', time_NNM))
infeas = sum(sum(isnan(X_L1))) + sum(sum(isnan(Y_L1)));
if infeas > 0
    disp(sprintf('Residual (NNM): NaN'))
    disp('feedback system: UNstable')
else
    [K1,Pcl1] = controller_realization(Ps,X_L1,Y_L1,nc);
    disp(sprintf('Residual (NNM): %f',norm(X_L1*Y_L1-eye(n),'fro')))
    if max(real(pole(Pcl1))) < 0
        disp('feedback system: stable')
    else
        disp('feedback system: UNstable')
    end
end
disp(sprintf('\n'))

% GS96
disp(sprintf('**GS96**'))
disp(sprintf('time(GS96): %f', time_GS))
[K2,Pcl2] = controller_realization(Ps,X1_GS,Y1_GS,nc);
disp(sprintf('Residual (GS96): %f', norm(X1_GS*Y1_GS-eye(n),'fro')))
if max(real(pole(Pcl2))) < 0
    disp('feedback system: stable')
else
    disp('feedback system: UNstable')
end
disp(sprintf('\n'))

% Proposed method
disp(sprintf('**Proposed**'))
disp(sprintf('time(proposed): %f', time_proposed))
[K3,Pcl3] = controller_realization(Ps,X1_P,Y1_P,nc);
disp(sprintf('Residual (proposed): %f', norm(X1_P*Y1_P-eye(n),'fro')))
if max(real(pole(Pcl3))) < 0
    disp('feedback system: stable')
else
    disp('feedback system: UNstable')
end
disp(sprintf('\n'))



return;
