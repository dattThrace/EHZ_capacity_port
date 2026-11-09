function [normals,inds,K,verToFace] = calcNormalsAndFaces(K)
n = size(K,2);
 
% inds = convhulln(K,{'QJ'});
% K = K(unique(inds(:)),:);
% [inds, vol] = convhulln(K,{'QJ'});

inds = convhulln(K,{'QJ'});
[uinds, ~, ic] = unique(inds(:));
inds(:) = ic;
K = K(uinds,:);

inds = sort(inds,2);
inds = sortrows(inds,1:size(inds,2));
removeInds = [];
% remove "faces" with rank smaller than n-1.
for i=1:size(inds,1)
    if (rank(K(inds(i,2:end),:) - repmat(K(inds(i,1),:),n-1,1),1e-10) < n-1)
        removeInds = [removeInds,i];
    end
end
inds(removeInds,:) = [];

normals = zeros(size(inds));
for i = 1:size(inds,1)
    normals(i,:) = getNormal(n,K,inds(i,:));
end

[normals,IA, IC]= uniquetol(normals,1e-8,'ByRows',true);

if (nargout > 3)
    verToFace = false(size(K,1),size(normals,1));
    for i = 1:size(inds,1)
        for k=1:length(inds(i,:))
            verToFace(inds(i,k),IC(i)) = true;
        end
    end
end

inds = inds(IA,:);
wasRemoved = ~isempty(removeInds) | length(IA) ~= size(normals,1);

% if (nargout > 3)
%     verToFace2 = zeros(size(K,1),n);
%     for i=1:size(K,1)
%         k = 1;
%         for j=1:size(inds,1)
%             if (rank(K(inds(j,:),:) - repmat(K(i,:),n,1),1e-10) < n)
%                 verToFace2(i,k) = j;
%                 k = k + 1;
%             end
%         end
%     end
%     for i=1:size(verToFace,1)
%         if ~all(sort(verToFace2(i,verToFace2(i,:) > 0)) == sort(find(verToFace(i,:))))
%             disp('oy oy oy !!')
%         end
%     end
% end

if (wasRemoved)
    disp('DEBUG: not all faces are simplices');
end

end

function normal = getNormal(n,K,inds)
tmp = K(inds(2:end),:) - repmat(K(inds(1),:),n-1,1);
normal = zeros(1,n);
for i = 1:n
    normal(i) = (-1)^i * det(tmp(:,[1:i-1,i+1:end]));
end
normal = normal./norm(normal);

%normal is orthogonal to the face, we are left to find its orientation.
det1 = det(K(inds,:) - repmat(K(inds(1),:)+normal,n,1));
tmp1 = setdiff(1:size(K,1),inds);
i = 1;
%find a vertex of K which is not on the current face.
while( rank(K(inds,:) - repmat(K(tmp1(i),:),n,1),1e-10) < length(inds))
    i = i + 1;
end
det2 = det(K(inds,:) - repmat(K(tmp1(i),:),n,1));
if (sign(det1) == sign(det2))
    normal = -normal;
end

end
