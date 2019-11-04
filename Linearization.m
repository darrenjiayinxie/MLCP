function [M,q,r,zk,wk] = Linearization(xk,l,F,dF)
n = size(l,1);
I = eye(n);

zk = projection_operator(xk,l);
wk = projection_operator(zk-xk,zeros(n,1));

M = dF(zk) + eps*I;
q = F(zk)-M(zk)*zk;
r= normal_map(zk);


end