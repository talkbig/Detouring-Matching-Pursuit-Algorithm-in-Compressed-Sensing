close all
clear
M=128;
N=256;
K=64;
tol=eps;

PHI=randn(M,N);
PHI=PHI./repmat(sum(PHI.^2),M,1);

rank=randperm(N);
rank=rank(1:K);
xtrue=zeros(N,1);
xtrue(rank)=randn(K,1);
y=PHI*xtrue;

PSI=PHI'*PHI;
zeta=[PHI,y]'*y;

[support,x,res_norm,sf]=GPA(PSI,zeta,N,K,tol);

stem([xtrue,x]);