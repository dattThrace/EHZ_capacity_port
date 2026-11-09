function C = barycenter(K)
%barycenter - Finds the barycenter of K
n = size(K,2);
%find triangulation of K
inds = delaunayn(K,{'QJ'});
TrgVols = zeros(size(inds,1),1);
% We find the barycenter of each simplex in the triangulation (the mean of the vertices) and take a 
% weighted mean (by the volume of the simplices) of those barycenters.
C = zeros(1,n);
for i = 1:size(inds,1)
    TrgVols(i) = abs(det(K(inds(i,2:end),:) - K(inds(i,1),:)));
    C = C + TrgVols(i)*sum(K(inds(i,:),:),1);
end
C = C/sum(TrgVols);
C = C/(n+1);
end