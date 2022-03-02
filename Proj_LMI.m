function [Xopt,Yopt] = Proj_LMI_4(X0,Y0,X,Y,LMI1)

ops=sdpsettings;
ops=sdpsettings(ops,'verbose',0);

[mx,nx] = size(X0);
[my,ny] = size(Y0);

U = blkdiag(X,Y);
U0 = blkdiag(X0,Y0);
gamm = sdpvar(1,1);
S = sdpvar(mx+my,mx+my);

M3 = [S,U-U0; (U-U0)', eye(nx+ny)];

LMI = [LMI1,M3>=0];

optimize(LMI,trace(S),ops);

Xopt = value(X);
Yopt = value(Y);

end

