function KXT = KtimesT(K,T)
KXT = [repmat(K,size(T,1),1), repelem(T,size(K,1),1)];
end