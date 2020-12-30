clc;
clear all;
close all;

global delta;
delta = 0.01;

a = -1; % left end point of interval
b = 1; % right
N = 12; %...so wachayall want?
h = 1/N;
x = a:h:b; % node coords

eta = ones(N,1);     % allocate element residuals
TOL = 1e-3;
MAX = 1e4;
alpha = 0.9;

% Will never reach 10000 nodes before reaching tolerance requirement of
% 1e-3



while N<MAX && sum(eta.^2)>TOL

% solution, assemble!
% stiffness matrix assembly
A = my_stiffness_matrix_assembler(x);

% load vector assembly
B = my_load_vector_assembler(x);

% mass matrix assembly
M = my_mass_matrix_assembler(x);

% solve system of equations   
xi = A\B;

% discrete laplace(u_h) approximation
Lxi = -M\A*xi;

%%%%%%%%%%%%%%%%%%%%
% compute residual %
%%%%%%%%%%%%%%%%%%%%
eta = zeros(N,1);     % allocate element residuals
for i = 1:length(x)-1 % loop over elements
    h = x(i+1) - x(i); % element length
    a1 = f(x(i))+Lxi(i); % temporary variables
    b1 = f(x(i+1))+Lxi(i+1);
    t = (a1^2+b1^2)*h/2; % integrate fˆ2. Trapezoidal rule
    eta(i) = sqrt(h*h*t); % element residual
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine select elements %
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:N
    if(eta(i)^2>alpha*max(eta.^2))      	% if residual is large
        x = [x (x(i+1)+x(i))/2];    % insert another node
    end
end
N = length(x);      % the size of the vector N increases
x = sort(x);        % however its inserted at the end, so we need to sort the vector

end


% solution, assemble!one last time!

% stiffness matrix assembly
A = my_stiffness_matrix_assembler(x);

% load vector assembly
B = my_load_vector_assembler(x);

% mass matrix assembly
M = my_mass_matrix_assembler(x);

% solve system of equations   
xi = A\B;

% discrete laplace(u_h) approximation
Lxi = -M\A*xi;

%%%%%%%%%%%%%%%%%%%%
% compute residual %
%%%%%%%%%%%%%%%%%%%%
eta = zeros(N,1);     % allocate element residuals
for i = 1:length(x)-1 % loop over elements
    h = x(i+1) - x(i); % element length
    a1 = f(x(i))+Lxi(i); % temporary variables
    b1 = f(x(i+1))+Lxi(i+1);
    t = (a1^2+b1^2)*h/2; % integrate fˆ2. Trapezoidal rule
    eta(i) = sqrt(h*h*t); % element residual
end

% evaluate the given function
for i=1:length(x)
    func(i,1) = f(x(i));
end

% residual
R = func+Lxi;

%plot everything
figure;
subplot(2,2,1)
plot(x,xi)
xlabel('x domain')
ylabel('Approximation u_h') 
title('FEM solution u_h');

subplot(2,2,2)
plot(x,R);
xlabel('x domain')
ylabel('Residual R')
title('Residual R = g+u_h'''', g(x) = delta^{-1} f(x)')

subplot(2,2,3)
plot(x,eta);
xlabel('x domain')
ylabel('Eta')
title('Error indicator eta, with g(x) = delta^{-1} f(x)')

subplot(2,2,4)
plot(x(2:end),[1./diff(x)])
xlabel('x domain')
ylabel('Meassure of node density')
title('mesh size distribution')



