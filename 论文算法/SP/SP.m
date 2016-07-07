function [support,x,res]=SP(PHI,y,K)

[~,N]=size(PHI);

y_r=y;
support=[];
support_old=ones(K,1);
res_old=sum(y_r.^2);
frame=ones(N,1);
tol=1;

while true
    [~,add_index]=sort(abs(y_r'*PHI),'descend');
    add_index=add_index(1:K);
    support=union(support,add_index);
    x_mmse=PHI(:,support)\y;
    [~,add_index]=sort(abs(x_mmse),'descend');
    support=support(add_index(1:K)); 
    x_mmse=PHI(:,support)\y;
    y_r=y-PHI(:,support)*x_mmse;
    res=sum(y_r.^2);
    if res<eps
        break;
    end
    rate=res/res_old;
    if rate>=1+tol || sum(frame(support))==0
        support=support_old;
        x_mmse=x_mmse_old;
        res=res_old;
        break;
    end
    if rate>=1
        tol=tol/2;
    end
    support_old=support;
    x_mmse_old=x_mmse;
    res_old=res;
    frame=ones(N,1);
    frame(support_old)=zeros(K,1);
end

x=zeros(N,1);
x(support)=x_mmse;