function A = stiffnessMatrixAssembler2DB1(p,t)

np = size(p,2);
nt = size(t,2);
A = sparse(np,np); % allocate stiffness matrix, why is it a sparse datatype and what is that 
for K = 1:nt
    loc2glb = t(1:3,K); % local-to-global map
    x = p(1,loc2glb); % node x-coordinates
    y = p(2,loc2glb); % node y-coordinates
    [area,b,c] = HatGradients(x,y); % calculates slopes on K
    AK = (b*b'+c*c')*area; % element stiffness matrix 
    A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK; % add locall stiffnesses to global A
end
end