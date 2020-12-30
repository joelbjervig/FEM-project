function A=my_stiffness_matrix_assembler(x)

    N = length(x)-1;
    A = zeros(N+1, N+1);    % allocate stiffnes matrix
    for i = 1:N % loop over elements
    	h = x(i+1) - x(i); % element length
        n = [i i+1]; % nodes
        A(n,n) = A(n,n) + [1 -1; -1 1]/h; % assemble element stiffness
    end
        
    A(1,1) = 1;
    A(1,2) = 0;
    A(N+1,N+1)= 1;
    A(N+1,N)=0;

end
