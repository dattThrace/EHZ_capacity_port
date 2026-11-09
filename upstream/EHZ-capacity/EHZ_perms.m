function cap = EHZ_perms(K,varargin)
% EHZ -- calculate the Ekeland-Hofer-Zehnder of a polytope K in the 2n Euclidean space.
% The coordinate system is of the form q_1,q_2,...,q_n,p_1,p_2,...,p_n.
% K - A matrix with the vertices of the polytope in its rows. 
% to disable the wait bar run EHZ(K,'ShowWaitBar',false).
showWB = true;
if (nargin > 1)
    if (strcmp(varargin{1},'ShowWaitBar'))
        showWB = varargin{2};
    end
end
wb = waitbar(0,'Staring...','Visible',showWB);
n = size(K,2)/2;
% translate K to be centered.
C = barycenter(K);
K = K - repmat(C,size(K,1),1);
% get the normals and which vertices are in each faces.
[normals,inds,K] = calcNormalsAndFaces(K);
heights = calcHeights(K,inds,normals);

k = size(inds,1);
allOmegas = zeros(k);
for i = 1:k
    for j = i+1:k
        allOmegas(i,j) = omega(normals([i,j],:),n);
        allOmegas(j,i) = -allOmegas(i,j);
    end
end
minCap = Inf;
timeEstimate = factorial(k)./((2:k).*factorial(k-(2:k)));
timeEstimate = cumsum(timeEstimate)./sum(timeEstimate);
for j = 2:k
    allChoices = nchoosek(1:k,j);
    allCaps = inf(size(allChoices,1),1);
    parfor i = 1:size(allChoices,1)
        currentPerm = allChoices(i,:);
        Aeq = [normals(currentPerm,:)'; heights(currentPerm)];
        [Aeq,beq,numOfSol] = prepareAeqAndBeq(Aeq,[zeros(2*n,1);1]);
        if (isinf(numOfSol))
            allCaps(i) = calculateForPerm(allOmegas(currentPerm,currentPerm),Aeq,beq);
        elseif(numOfSol == 1)
            betaArr = Aeq\beq;
            if (all(betaArr > - 1e-10))
                allCaps(i) = calculateForPermSingleSol(allOmegas(currentPerm,currentPerm),betaArr);
            end
        end
%         if (mod(i,1000) == 0)
%             disp(i);
%         end
    end
    minCap = min(minCap,min(allCaps));
    waitbar(timeEstimate(j-1),wb,['Please wait (',num2str(j),'/',num2str(k),') ...']);
end

cap = -1/2/minCap;
close(wb);
end

function [minCap, minPerms] = calculateForPerm(allOmegas,AeqOrg,beq)
routes = [ones(factorial(size(allOmegas,1)-1),1),perms(2:size(allOmegas,1))];
allCaps = inf(factorial(size(allOmegas,1)-1),1);
for i = 1:size(routes,1)
    currentPerm = routes(i,:);
    H = allOmegas(currentPerm,currentPerm);
    H = H - triu(H) + tril(H)';
    Aeq = AeqOrg(:,currentPerm);
    matWithLagrangeCoef = zeros(size(H,1)+size(Aeq,1));
    matWithLagrangeCoef(1:size(H,1),1:size(H,1)) = H;
    matWithLagrangeCoef(1:size(H,1),size(H,1)+1:end) = Aeq';
    matWithLagrangeCoef(size(H,1)+1:end,1:size(H,1)) = Aeq;
    % If the matrix is not of full rank, this means that either there is
    % no solution or there are infinite solutions. If there are infinite
    % solutions, one can find a solution where \beta_i=0 for some i, and
    % we will find this solution using a smaller length permutation.
    if (abs(det(matWithLagrangeCoef)) < 1e-10)
        continue;
    end
    betaArr = matWithLagrangeCoef\[zeros(size(H,1),1);beq];
    betaArr = betaArr(1:size(H,1));
    newCap = betaArr' * H * betaArr/2;
    if (all(betaArr > - 1e-10))
        allCaps(i) = newCap;
    end
end
minCap = min(allCaps);
if (isinf(minCap))
    minPerms = [];
else
    minPerms = routes(abs(allCaps-minCap)<1e-8,:);
end
end

function minCap = calculateForPermSingleSol(allOmegas,betaArr)
routes = [ones(factorial(size(allOmegas,1)-1),1),perms(2:size(allOmegas,1))];
allCaps = inf(factorial(size(allOmegas,1)-1),1);
for i = 1:size(routes,1)
    currentPerm = routes(i,:);
    H = allOmegas(currentPerm,currentPerm);
    H = H - triu(H) + tril(H)';
    newCap = betaArr(currentPerm)' * H * betaArr(currentPerm)/2;
    allCaps(i) = newCap;
end
minCap = min(allCaps);
end

function res = omega(V,n)
res = V(1,1:n)*V(2,n+1:2*n)' - V(1,n+1:2*n)*V(2,1:n)';
end
