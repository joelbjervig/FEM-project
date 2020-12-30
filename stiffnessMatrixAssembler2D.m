function A = stiffnessMatrixAssembler2D(p,t)
    np = size(p,2); % number of nodes
    nt = size(t,2); % number of elements

    A = zeros(np,np); % allocate stiffness matrix
    
    for i = 1:nt    % loop over elements
        
        nodes = t(1:3,i);   % local matrix to global matrix mapping
        
        x = p(1,nodes);     % node x-coordinates
        y = p(2,nodes);     % node y-coordinates
        
        [area,b,c] = HatGradients(x,y); % compute gradients in x&y direction
                                        % and the area of traingle K.
                                        
        localA = (b*b'+c*c').*area;     % local element stiffness matrix
        
        
        A(nodes,nodes) = A(nodes,nodes) + localA;   % add local matrix stiffnesses
                                                    % to global matrix A
    end
end