 clear all
%[M,q,r,zk,wk,l,u] = parameters();
N = 30;

M = eye(N)-rand(N,N);
q = rand(N,1);
zk = rand(N,1);
wk = zeros(N,1);
l = -1e15*randi([0 1],N,1);
u = Inf*ones(N,1);
r = -(M*zk-eye(N)*wk+q);



% load('bugcycling.mat')
% M = A.M;
% q = A.q;
% l = A.l;
% u = A.u;
% zk = A.zk;
% wk = A.wk;
% r = A.r;

tic
[zk_n,wk_n,A] = MLCP(M,q,r,zk,wk,l);
toc
residual = M*zk_n-wk_n+q;

% tic
% z = pathlcp(M,q,l,u,zk);
% toc
