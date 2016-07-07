function [delta,add_factor,add_point,add_index]=DMP_max(PSI,coff_matrix,support,complement)

self_relative=diag(PSI);

factor=self_relative(complement)-sum(PSI(support,complement).*coff_matrix(:,1:end-1),1)';
[delta,add_index]=max((PSI(complement,end)-(PSI(end,support)*coff_matrix(:,1:end-1))').^2./factor);
add_point=complement(add_index);
add_factor=factor(add_index);