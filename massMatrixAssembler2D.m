function M = massMatrixAssembler2D(p,t)
    np = size(p,2); % number of nodes
    nt = size(t,2); % number of elements

    M = zeros(np,np); % allocate mass matrix

    for i = 1:nt % loop over elements
        
        nodes = t(1:3,i);   % local matrix to global matrix mapping
        
        x = p(1,nodes);     % node x-coordinates
        y = p(2,nodes);     % node y-coordinates

        area = polyarea(x,y);   % triangle area with matlabs built-in fucntion
        
        
        localM =    [2 1 1;     % local element mass matrix
                     1 2 1;
                     1 1 2]/12*area; 
        M(nodes,nodes) = M(nodes,nodes)+ localM;    % add local element mass matrices
                                                    % to the global matrix M
    end
end
