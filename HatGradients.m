function [area,b,c] = HatGradients(x,y)
    % computing area of triangle K
    area=polyarea(x,y);
    
    %Computing the slope
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area; % slope in x direction
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area; % slope in y direction
    end