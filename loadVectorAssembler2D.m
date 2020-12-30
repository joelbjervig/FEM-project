function b = loadVectorAssembler2D(p,t,u)
    np = size(p,2); % number of nodes
    nt = size(t,2); % number of elements
    
    b = zeros(np,1);% allocate load vector
    
    for i = 1:nt    % loop over all elements
        nodes = t(1:3,i);   % local vector to global vector mapping
        
        x = p(1,nodes);     % node x-coordinates
        y = p(2,nodes);     % node y-coordinates
        
        area = polyarea(x,y);   % triangle area with matlabs built-in fucntion
        
        localb = S(u(nodes))/3*area;    % local element load vector
        
        b(nodes) = b(nodes) + localb;   % add local element load vector
                                        % to global load vector b
     end
end