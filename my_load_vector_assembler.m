function B=my_load_vector_assembler(x)

    N = length(x)-1;
    B = zeros(N+1, 1);      % allocate load vector
    for i = 1:N
            h = x(i+1) - x(i);
            n = [i i+1];
            B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
    end
    B(1) = 0;   % Apply BC strongly
    B(end) = 0; % Apply BC strongly
end