function [support,x,res,sf]=DMP(PSI,N,K,tol)

[coff_matrix,support,complement]=DMP_add(PSI,N,K);
x_mmse=coff_matrix(:,end);
res=PSI(end,end)-PSI(end,support)*x_mmse;
sf=(res<tol);

if sf==0
    product_inv=inv(PSI(support,support));
    [support,x_mmse,res,sf]=DMP_back_forth(PSI,coff_matrix,product_inv,support,complement,K,tol);
end   
    
x=zeros(N,1);
x(support)=x_mmse;