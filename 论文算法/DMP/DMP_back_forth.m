function [support_update,x_mmse,res_update,sf]...
    =DMP_back_forth(PSI,coff_matrix,product_inv,support,complement,K,tol)

[~,N]=size(PSI);
N=N-1;
support_update=support;
x_mmse=coff_matrix(:,end);
res_update=PSI(end,end)-PSI(end,support)*coff_matrix(:,end);
sf=0;

back=0;
remove_point=zeros(K,1);
add_point=zeros(K,1);
while back<K-1 && back>=0
    back=back+1;
    [~,remove_point(back),remove_index]=DMP_min(coff_matrix,product_inv,support);
    [coff_matrix,product_inv,support,complement]...
        =DMP_back(PSI,coff_matrix,product_inv,support,complement,remove_point(back),remove_index);
    if sum(support~=support_update(1:K-back))==0
        continue;        
    end    
    coff_matrix_new=coff_matrix;
    product_inv_new=product_inv;
    support_new=support;
    complement_new=complement;
    update=1;
    frame=zeros(N,1);
    for t=back:-1:1
        [~,factor,add_point(t),add_index]=DMP_max(PSI,coff_matrix_new,support_new,complement_new);
        frame(remove_point(t))=1;
        if sum(frame(add_point(t:back))==0)==0
            update=0;
            break;
        end
        [coff_matrix_new,product_inv_new,support_new,complement_new]=...
            DMP_forth(PSI,coff_matrix_new,product_inv_new,support_new,complement_new,factor,add_point(t),add_index);
    end
    if update==1
        res=PSI(end,end)-PSI(end,support_new)*coff_matrix_new(:,end);
        if res<res_update
            support_update=support_new;
            res_update=res;
            x_mmse=coff_matrix_new(:,end);
            if res<tol
                sf=1;
                break;
            end
            coff_matrix=coff_matrix_new;
            product_inv=product_inv_new;
            support=support_new;
            complement=complement_new;
            back=0;
        end 
    end
end