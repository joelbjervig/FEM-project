close all; clear all; clc;
% define domain
geometry = @circleg ;
hmax = [5 20 40].^(-1); % REMOVE FIRST ELEMENT IN FINAL 

% define constants
delta = 0.01;       % physical meaning?
alpha = 4;          % physical meaning?
k = 0.01;           % temporal resolution
T = 2;              % timelimit

% define vectors
populations=zeros(length(hmax),T/k+1); %store population rate for each h
a=1;    %indexing for populationrates, increases by one for each h
time = [0:k:T];     % time vector


    
% spin it!
for h=hmax
    
    % Initialize mesh
	[p,e,t]  = initmesh (geometry , 'hmax' , h);
    
    % begin a figure for the solutions
    figure
    sgtitle(['Solutios at times t = 0, 0.5, 2, with meshsize h = 1/', num2str(h^(-1))],'fontsize',16)
    
    u_h = zeros(length(p),1); % create empty vector for all solutions
    
    % impose initial condition
    for i=1:length(p)
        u_h(i) = InitialCondition();
    end
    
    % Initial condition solution plot
    subplot(1,3,1)
    pdeplot(p,e,t, 'xydata', u_h, 'zdata', u_h)
    title('Solution at time t = 0','fontsize',14)
    
    % Matrices, Assemble!
    [A,M,b] = matrixAssembler(p,t,u_h); % u_h set from initial condition above
    
    % CRANK-NICHOLSON TIME DISCRETIZATION
    LHS = M/k+0.5*delta*A;     % matrix computation for the upcoming timestep u_{h}^{n+1}
    RHS = M/k-0.5*delta*A;     % matrix copmutation for the current timestep u_{h}^{n}
    
    population = [];% create empty vector for storage of population
    for i=time

        u_h = LHS\(RHS*u_h-M*S(u_h));  % solve for the upcoming u_{u,n+1},
                                % using the matrices and the previous solution u_{h,n}

        
        % population
        population = [population  Population(p,e,t,u_h)];
%         pdeplot(p,e,t, 'xydata', u_h, 'zdata', u_h)
%         pause(0.001);
        
        % solution plot halfway at t = 0.5
              if(i==0.5)
            subplot(1,3,2)
            pdeplot(p,e,t, 'xydata', u_h, 'zdata', u_h)
            title('solution at time t = 0.5','fontsize',14)
        end
    end
    
    %solution plot at ending, t=2
        subplot(1,3,3)
        pdeplot(p,e,t, 'xydata', u_h, 'zdata', u_h)
        title('solution at time t = 2.','fontsize',14)
        pause;
     
     populations(a,:)=population;    % store all populationvalues
                                    	% for all meshsizes
     a = a+1;    % increace indexing
end

% differentiate the population vector to obtain the rate
% populationrates = zeros(length(hmax),length(populations));
% for j=2:length(population)
%     populationrates(:,j) = (populations(:,j)-populations(:,j-1))/k;
% end

% Task is to only plot the population for the two first meshes
figure;
plot(time(2:end), populations(1,2:end),time(2:end), populations(2,2:end),time(2:end), populations(3,2:end))
legend({'h = 1/5','h = 1/20', 'h = 1/40'});
title('Population of species over time','fontsize',16)
xlabel('Time')
ylabel('Population')

% sgtitle('Population of species over time','fontsize',16)
% subplot(3,1,1)
% plot(time(2:end), populations(1,2:end))
% title(['Meshsize h = 1/',num2str(hmax(1)^(-1))],'fontsize',14)
% xlabel('time in seconds','fontsize',14)
% ylabel('Population','fontsize',14)
% subplot(3,1,2)
% plot(time(2:end), populations(2,2:end))
% title(['Meshsize h = 1/', num2str(hmax(2)^(-1))],'fontsize',14)
% xlabel('time in seconds','fontsize',14)
% ylabel('Population','fontsize',14)
% subplot(3,1,3)
% plot(time(2:end), populations(3,2:end))
% title(['Meshsize h = 1/', num2str(hmax(3)^(-1))],'fontsize',14)
% xlabel('time in seconds','fontsize',14)
% ylabel('population,','fontsize',14)

