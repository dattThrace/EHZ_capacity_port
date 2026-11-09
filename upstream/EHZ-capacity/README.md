# EHZ-capacity

This is an implementation for calculating the EHZ capacity of  a polytope in Matlab using the formula appearing in https://link.springer.com/article/10.1007/s00039-019-00486-4.
A polytope K in R^{2n} with m vertices is represented as a matrix with m rows and 2n cols, such that each row is a vertex. The coordinate system is q_1,...q_n,p_1,...,p_n.

Currently, this implementation is not very efficient and it is very slow when the polytope has a large number of faces.
However, one can significantly improve the running time by eliminating many permutations based on the fact that the minimal permutation should correspond to a closed characteristic on the boundary of K (see Remark 3.11 in the aforementioned paper).
I will upload an updated version soon.

This software was developed in Matlab R2021b.

An example for calculating the capacity of the Lagrangian product of the standard cube and the cross polytope
>> Q = [-1 -1; -1 1; 1 -1; 1 1];
>> 
>> CP = [eye(2); -eye(2)];
>> 
>> QXCP = KtimesT(Q,CP);
>> 
>> c = EHZ_perms(QXCP)
