function [support,x,res,sf]=GPA(PSI,zeta,N,M,K,tol)

support=[];
support=GPA_add(PSI,zeta,support,N,K);
x_mmse=PSI(support,support)\zeta(support);
res=zeta(end)-zeta(support)'*x_mmse;

if res<tol
    sf=1;
    zoom=0;
else
    sf=0;
    zoom=min(K,M-K-2);
end
while zoom>0
    support_new=GPA_add(PSI,zeta,support,N,zoom);
    [support_new,x_mmse_new]=GPA_remove(PSI,zeta,support_new,zoom);
    res_new=zeta(end)-zeta(support_new)'*x_mmse_new;
    if res_new>=res 
        zoom=zoom-1;
    else
        support=support_new;
        res=res_new;
        x_mmse=x_mmse_new;
        if res<=tol
            sf=1;
            break;
        end
    end
end

x=zeros(N,1);
x(support)=x_mmse;


