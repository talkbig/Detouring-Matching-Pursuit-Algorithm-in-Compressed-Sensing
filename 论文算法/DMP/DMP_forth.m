function [coff_matrix_new,product_inv_new,support_new,complement_new]=...
    DMP_forth(PSI,coff_matrix,product_inv,support,complement,factor,point,index)

complement(index)=[];
complement_new=complement;
coff_col=coff_matrix(:,index);
coff_matrix(:,index)=[];
coff_row=(PSI(point,[complement;end])-PSI(point,support)*coff_matrix)/factor;
coff_matrix_new=[coff_matrix-coff_col*coff_row;coff_row];
support_new=[support;point];        
coff_factor=coff_col/factor;
product_inv_new=[product_inv+coff_factor*coff_col',-coff_factor;-coff_factor',1/factor];
