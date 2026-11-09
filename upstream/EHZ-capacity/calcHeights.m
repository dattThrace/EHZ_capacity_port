function heights = calcHeights(K,inds,normals)
heights = zeros(1,size(inds,1));
for i = 1:size(inds,1)
    heights(i) = dot(K(inds(i,1),:),normals(i,:));
end
end