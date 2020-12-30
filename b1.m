clear all;
close all;

% define geometry and meshsize vector to be iterate over
geometry = @circleg ;
hmax = [2 4 8 16 32].^-1;

% function handles
f = @(x,y)8*pi^2*sin(2*pi*x)*sin(2*pi*y);
u_exact = @(x,y)sin(2*pi*x).*sin(2*pi*y);

EnE = []; % create empty vector to store values of the energynorm in
for h = hmax

    [p,e,t]  = initmesh ( geometry , 'hmax' , h );	% initialize mesh
    
    A = stiffnessMatrixAssembler2DB1(p,t);          % assemble stiffnessmatrix
    b = loadVectorAssembler2DB1(p,t,f);             % assebmle loadvector
    
    Uexact  = zeros(length(b),1);   % construct a zerovector for the exact solution on u  
    I = eye(length(p));             % construct the identity matrix
    A(e(1 ,:) ,:) = I(e(1 ,:) ,:);  % replace the rows corresponding to the boundary nodes by corresponding rows of I
    b(e(1 ,:)) = u_exact(p(1,e(1,:)),p(2,e(1,:))); % put the boundary value into the RHS
    u_h = A\b;  % solve the system of equations

    % put the exact solution for all points in a vector
    for i = 1:length(b)
        Uexact(i) = u_exact(p(1,i), p(2,i));
    end
    
    %error
    err = u_h-Uexact;
    
    % Energynorm of error
    EnE = [EnE, sqrt(err'*A*err)];
    
    % plot for coarsest and finest meshes
    if(h == hmax(1)) % if meshsize is large: 1/2
        subplot(1,2,1)
        pdeplot(p,e,t, 'xydata', u_h, 'zdata', u_h) % then plot solution
        title('Solution of simplified pray population equation, meshsize h = 1/2','fontsize', 16 )
        hold on
    elseif(h == hmax(end)) % if meshsize is small: 1/32
        subplot(1,2,2)
        pdeplot(p,e,t, 'xydata', u_h, 'zdata', u_h) % then plot solution
        title('Solution of simplified pray population equation, meshsize h = 1/32','fontsize', 16)
    end

end

%plot results
figure;
plot(log(hmax),log(EnE), 'o')
hold on;
plot(log(hmax), log(hmax.^(1.409))+2.381)
hold off;
title('Energynorm of the error. Convergencerate p = 1.409', 'fontsize', 16)
xlabel('log( h_{max} )', 'fontsize',14);
ylabel('log( ||u-u_{h}||_e )', 'fontsize',14);



