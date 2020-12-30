function M=my_mass_matrix_assembler(x)
    N = length(x)-1;
    M = zeros(N+1,N+1);     % allocate mass matrix
    for i = 1:N % loop over subintervals
        h = x(i+1) - x(i); % interval length

        M(i,i) = M(i,i) + h/3; % add h/3 to M(i,i)
        M(i,i+1) = M(i,i+1) + h/6;
        M(i+1,i) = M(i+1,i) + h/6;
        M(i+1,i+1) = M(i+1,i+1) + h/3;
    end

end
