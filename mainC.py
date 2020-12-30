#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 10:35:53 2020

@author: joelbjervig
"""

import matplotlib.pyplot as plt
from fenics import *
import random

# Create mesh and define function space
mesh = Mesh("circle.xml")
#x = SpatialCoordinate(mesh)

# Construct the finite element space polynomial of degree 1
V = VectorFunctionSpace(mesh, 'P', 1)

# Define parameters
T = 1000
dt = 0.5
alpha = 0.4
beta = 2
gamma = 0.8
delta1 = 1
delta2 = 1
theta = 0.5


# Class representing the intial conditions TASK 1
class InitialConditions1(UserExpression):
    def eval (self , values , x):
        values[0] = (4/15)-2*pow(10,-7)*(x[0]-0.1*x[1]-255)*(x[0]-0.1*x[1]-675)
        values[1] = (22/45)-3*pow(10,-5)*(x[0]-450)-1.2*pow(10,-4)*(x[1]-150)
        
    def value_shape(self):
        return (2 ,)
   
# Class representing the intial conditions TASK 2
class InitialConditions2(UserExpression):
    def eval (self , values , x):
        r = random.uniform(0,1)
        values[0] = 0.5*(1-r)
        values[1] = 0.25+0.5*r
        
    def value_shape(self):
        return (2 ,)
    
# test and trial functions
u, v = TrialFunction(V), TestFunction(V)

# Define initial condition chosse IC1 or IC2 depending on task
u_initial = Function(V)
u_initial = InitialConditions2(degree=2)

# interpolate initial condition
u0 = Function(V)
u0 = project(u_initial, V) # or project u initial onto V space as project(u_initial)


# Create bilinear and linear forms. first term from backward euler
a0 = inner(u[0],v[0])*dx - dt*theta*inner(u[0],v[0])*dx + dt*theta*delta1*inner(grad(u[0]),grad(v[0]))*dx
a1 = inner(u[1],v[1])*dx + dt*theta*gamma*inner(u[1],v[1])*dx + dt*theta*delta2*inner(grad(u[1]),grad(v[1]))*dx
# a0 = u[0]*v[0]*dx - dt*u[0]*v[0]*dx + theta*delta1*dt*inner(grad(u[0]),grad(v[0]))*dx
# a1 = u[1]*v[1]*dx - dt*gamma*u[1]*v[1]*dx + theta*delta2*dt*inner(grad(u[1]),grad(v[1]))*dx 

L0 = inner(u0[0],v[0])*dx + dt*(1-theta)*inner(u0[0],v[0])*dx - dt*(1-theta)*delta1*inner(grad(u0[0]),grad(v[0]))*dx - dt*inner(u0[0]**2 + u0[0]*u0[1]/(u0[0]+alpha),v[0])*dx
L1 = inner(u0[1],v[1])*dx - dt*(1-theta)*gamma*inner(u0[1],v[1])*dx - dt*(1-theta)*delta2*inner(grad(u0[1]),grad(v[1]))*dx + dt*beta*inner(u0[0]*u0[1]/(u0[0]+alpha),v[1])*dx


# L0 = u0[0]*v[0]*dx - theta*delta1*dt*inner(grad(u0[0]),grad(v[0]))*dx + ((-u0[0]**2)+(u0[0]*u0[1])/(u0[0]+alpha))*v[0]*dx
# L1 = u0[1]*v[1]*dx - theta*delta2*dt*inner(grad(u0[1]),grad(v[1]))*dx + ((u0[0]*u0[1])/(u0[0]+alpha))*v[1]*dx

# L0 = u0[0]*v[0]*dx - theta*delta1*dt*inner(grad(u0[0]),grad(v[0]))*dx + u0[0]*(1-u0[0])*v[0]*dx
# L1 = u0[1]*v[1]*dx - theta*delta2*dt*inner(grad(u0[1]),grad(v[1]))*dx + u0[0]*u0[1]*v[1]*dx

a = a0+a1
L = L0+L1
# set up boundary conditions
# Because this is a homogeneous Neumann problem, no boundary conditions
# are specified.
# bc = []

# assemble stiffness matrix
A = assemble(a)


# Set output file
out_file = File("results/poisson1.pvd","compressed")
# open txt file to write population in
population_file = open("results/preyPredPop.txt","w+")

# Set initial condition
u = Function(V)
u.assign(u0) #u = u0


u_initial = Function(V)
u_initial.assign(u0) 

t_save = 0;
num_samples = 1000

t = 0.0 # initial time
# times to plot at
t0 = 0
t1 = 50
t2 = 100
t3 = 150
t4 = 500
t5 = 1000

# remember to apply the initial condition corresponding to each task.
task1 = True  # True = run task 1
task12 = False # True = run task 1.2: second instance of task 1
task2 = False  # True = run task 2


# initial solution, basically plotting IC
out_file<<(u,t)
# categories
population_file.write("Time, Prey, Predator\n")

# Timestepping
while t <= T:

    # u0 = u
    u0.assign(u)
    
    # assemble mattrix and vector
    b = assemble(L)
    
    # solve linear system with LU factorization
    solve(A,u.vector(),b,"lu", "default")
    
    t_save += dt
    # only save a certain amount of samples, might get big otherwise
    if t_save > T/num_samples or t >= T-dt:
        print("saving solution")
        
        # save file
        out_file << (u,t)
        
        t_save = 0
        
    # PRINTING
    if ( (task1==True and (t == t0 or t == t1 or t == t2 or t == t3)) or ((task12 == True) and (t==t5)) or ((task2 == True) and (t == t0 or t == t1 or t == t2 or t == t4 or t == t5)) ):
        # plot solution prey
        plot(u[0])
        plt.title('Population density of prey at time t = %i' % t)
        plt.xlabel('X-direction')
        plt.ylabel('Y-direction')
        plt.show() 
        # plot solution predators    
        plot(u[1])
        plt.title('Population density of predators at time t = %i' % t)
        plt.xlabel('X-direction')
        plt.ylabel('Y-direction')
        plt.show()
       
    # move to next timestep and adjust boundary condition.
    t = t + dt
    
    # compute the functional
    population_u = assemble(u[0] * dx)
    population_v = assemble(u[1] * dx)


    # print population for prey and predator
    print("",t ,", ",population_u,", ",population_v)
    tempstring = "" + str(t) + ", " + str(population_u) + ", " + str(population_v) + "\n"
    population_file.write(tempstring)



