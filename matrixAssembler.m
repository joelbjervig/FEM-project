function [A,M,b] = matrixAssembler(p,t,u)

    % call all assemblers
    A = stiffnessMatrixAssembler2D(p,t);
    M = massMatrixAssembler2D(p,t);
    b = loadVectorAssembler2D(p,t,u);
end