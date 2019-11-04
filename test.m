
clear all
%A = parameters(A);
N = 10;

M = eye(N)+rand(N,N);
q = rand(N,1);
zk = rand(N,1);
wk = zeros(N,1);
l = -1e15*randi([0 1],N,1);
u = Inf*ones(N,1);
r = -(M*zk-eye(N)*wk+q);





tic
[zk_n,wk_n,A] = MLCP(M,q,r,zk,wk,l);
toc
residual = M*zk_n-wk_n+q;

tic
z = pathlcp(M,q,l,u,zk);
toc
%         