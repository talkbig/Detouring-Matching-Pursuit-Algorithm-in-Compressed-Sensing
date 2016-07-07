function [delta,remove_point,remove_index]=DMP_min(coff_matrix,product_inv,support)

[delta,remove_index]=min(coff_matrix(:,end).^2./diag(product_inv));
remove_point=support(remove_index);