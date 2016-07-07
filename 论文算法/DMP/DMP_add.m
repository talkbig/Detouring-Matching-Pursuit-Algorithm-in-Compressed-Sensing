function [coff_matrix,support,complement]=DMP_add(PSI,N,zoom)

self_relative=diag(PSI);
[~,add_index]=max((PSI(1:end-1,end).^2./(self_relative(1:end-1))));
support=add_index;
complement=(1:N)';
complement(add_index)=[];
coff_matrix=PSI(add_index,add_index)\PSI(add_index,[complement;end]);
if zoom~=1
    for t=1:zoom-1
        factor=self_relative(complement)-sum(PSI(support,complement).*coff_matrix(:,1:end-1),1)';
        [~,add_index]=max((PSI(complement,end)-(PSI(end,support)*coff_matrix(:,1:end-1))').^2./factor);
        add_point=complement(add_index);
        complement(add_index)=[];
        coff_col=coff_matrix(:,add_index);
        coff_matrix(:,add_index)=[];
        coff_row=(PSI(add_point,[complement;end])-PSI(add_point,support)*coff_matrix)/factor(add_index);
        coff_matrix=[coff_matrix-coff_col*coff_row;coff_row];
        support=[support;add_point];        
    end
end
