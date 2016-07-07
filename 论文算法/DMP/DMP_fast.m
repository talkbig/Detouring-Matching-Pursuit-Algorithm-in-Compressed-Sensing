close all
clear
M=128;
N=256;
K=32;
tol=eps;

PHI=randn(M,N);
PHI=PHI./repmat(sum(PHI.^2),M,1);

rank=randperm(N);
rank=rank(1:K);
xtrue=zeros(N,1);
xtrue(rank)=randn(K,1);
%xtrue(rank)=ones(K,1);
%xtrue(rank)=(1:K)';
y=PHI*xtrue;

PSI=[PHI,y]'*[PHI,y];

[support,x,res_norm,sf]=DMP(PSI,N,K,tol);

stem([xtrue,x]);