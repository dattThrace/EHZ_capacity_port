function [Aeq,beq,numOfSol] = prepareAeqAndBeq(Aeq,beq)
Aeq(abs(Aeq)<1e-10) = 0;
[~,R,p] = qr(Aeq','vector');
AeqRank = find(abs(diag(R)) > 1e-8,1,'last');
[~,R2,p2] = qr([Aeq,beq]','vector');
if (AeqRank == find(abs(diag(R2)) > 1e-8,1,'last'))
    %dispose of linearly dependant constraints.
    I = p(1:AeqRank);
    Aeq = Aeq(I,:);
    beq = beq(I);
    if (size(Aeq,1) == size(Aeq,2) && size(Aeq,1) == AeqRank)
        numOfSol = 1;
    else
        numOfSol = Inf;
    end
else
    %beq is not in the span of Aeq.
    numOfSol = 0;
end
end