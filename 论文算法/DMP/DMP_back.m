function [coff_matrix_new,product_inv_new,support_new,complement_new]...
    =DMP_back(PSI,coff_matrix,product_inv,support,complement,point,index)

product_element=product_inv(index,index);
product_row=product_inv(index,:);
product_inv(index,:)=[];
product_row(index)=[];
product_col=product_inv(:,index);
product_inv(:,index)=[];
product_inv_new=product_inv-product_col*(product_row/product_element);
complement_new=[point;complement];
support(index)=[];
support_new=support;

coff_row=coff_matrix(index,:);
coff_matrix(index,:)=[];
coff_col=product_inv_new*PSI(support_new,point);
coff_matrix=coff_matrix+coff_col*coff_row;
coff_matrix_new=[coff_col,coff_matrix];


