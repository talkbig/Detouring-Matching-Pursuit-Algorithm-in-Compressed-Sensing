function [support_new,x_mmse] = GPA_remove(PSI,zeta,support,zoom)

product_inv=inv(PSI(support,support));

for t=1:zoom
    x_mmse=product_inv*zeta(support);
    [~,remove_index]=min(x_mmse.^2./diag(product_inv));
    product_element=product_inv(remove_index,remove_index);
    product_row=product_inv(remove_index,:);
    product_inv(remove_index,:)=[];
    product_row(remove_index)=[];
    product_col=product_inv(:,remove_index);
    product_inv(:,remove_index)=[];
    product_inv=product_inv-product_col*(product_row/product_element);
    support(remove_index)=[];
end
support_new=support;
x_mmse=product_inv*zeta(support);
    
    