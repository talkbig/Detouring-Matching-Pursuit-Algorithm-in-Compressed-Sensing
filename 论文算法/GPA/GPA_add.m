function support_new = GPA_add(PSI,zeta,support,N,zoom)

self_relative=diag(PSI);
if isempty(support)
    [~,add_index]=max(zeta(1:N).^2./self_relative);
    support=add_index;
    zoom=zoom-1;   
end
if zoom ~= 0
    complement=(1:N)';
    complement(support)=[];
    coff_matrix=PSI(support,support)\PSI(support,:);
    for t=1:zoom
        factor=self_relative(complement)-sum(PSI(support,complement).*coff_matrix(:,complement),1)';
        [~,add_index]=max((zeta(complement)-(zeta(support)'*coff_matrix(:,complement))').^2./factor);
        add_point=complement(add_index);
        coff_row=(PSI(add_point,:)-coff_matrix(:,add_point)'*PSI(support,:))/factor(add_index);
        coff_matrix=coff_matrix-coff_matrix(:,add_point)*coff_row;
        coff_matrix=[coff_matrix;coff_row]; 
        support=[support;add_point];
        complement(add_index)=[];
    end
end   
support_new=support;
