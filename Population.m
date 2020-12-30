function population = Population(p,e,t,u)
    population = 0;         % initialize population count
    for i=1:size(t,2)
        nodes = t(1:3,i);   % local to global mapping
        x = p(1,nodes);     % nodes x-coordinates
        y = p(2,nodes);     % nodes y-coordinates
        
        area = polyarea(x,y);   % triangle area with matlabs built-in function
        
        % increasing population by discrete integration ocer triangle by trapezodial rule
        population = population + area*(u(nodes(1)) + u(nodes(2)) + u(nodes(3)))/3;
    end
end