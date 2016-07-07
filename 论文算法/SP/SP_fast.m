close all
clear
M=128;
N=256;
K=40;

PHI=randn(M,N);
PHI=PHI./repmat(sum(PHI.^2),M,1);

rank=randperm(N);
rank=rank(1:K);
xtrue=zeros(N,1);
xtrue(rank)=randn(K,1);
y=PHI*xtrue;

tic;
[support,x,res_norm]=SP(PHI,y,K);
time=toc

stem([xtrue,x]);