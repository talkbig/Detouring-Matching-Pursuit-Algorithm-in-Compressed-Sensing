function [support,x]=OMP(PHI,y,K)

[~,N]=size(PHI);
x_mmse=zeros(K,1);
support=zeros(K,1);
y_res=y;

for t=1:K
    proj=PHI'*y_res;
    [~,pos]=max(abs(proj));
    support(t)=pos;
    x_mmse(1:t)=PHI(:,support(1:t))\y;
    y_res=y-PHI(:,support(1:t))*x_mmse(1:t);
end

x=zeros(N,1);
x(support)=x_mmse;