#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:56:45 2020

@author: jenstandstad
"""

import numpy as np
import matplotlib.pyplot as plt
# Modified but very slightly stolen from 
# Solving Ordinary Differential Equations in Python, by Joakim Sundnes 
Increment = 1
MinT = 0
MaxT = 8
N = int( (MaxT - MinT)/Increment) 
InitialCondition = 0.1

"""
Version 1 of the ODESolver class hierarchy. This is the simplest version
of the class, which only works for scalar ODEs.
"""

class ODESolver:
    def __init__(self, f):
        # Wrap user's f in a new function that always
        # converts list/tuple to array (or let array be array)
        self.f = lambda u, t: np.asarray(f(u, t), float)

    def set_initial_condition(self, U0):
        if isinstance(U0, (float,int)):  # scalar ODE
            self.neq = 1                 # no of equations
            U0 = float(U0)
        else:                            # system of ODEs
            U0 = np.asarray(U0)
            self.neq = U0.size           # no of equations
        self.U0 = U0

    def solve(self, time_points):
        self.t = np.asarray(time_points)
        N = len(self.t)
        if self.neq == 1:  # scalar ODEs
            self.u = np.zeros(N)
        else:              # systems of ODEs
            self.u = np.zeros((N,self.neq))

        # Assume that self.t[0] corresponds to self.U0
        self.u[0] = self.U0

        # Time loop
        for n in range(N-1):
            self.n = n
            self.u[n+1] = self.advance()
        return self.u, self.t


class ForwardEuler_v2(ODESolver):
    def advance(self):
        """Advance the solution one time step."""
        # Create local variables to get rid of "self." in
        # the numerical formula
        u, f, n, t = self.u, self.f, self.n, self.t
        #dt is not necessarily constant:
        dt = t[n+1]-t[n]
        unew = u[n] + dt*f(u[n], t[n])
        return unew

class ForwardEuler(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        dt = t[n+1] - t[n]
        unew = u[n] + dt*f(u[n], t[n])
        return unew

class ExplicitMidpoint(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        dt2 = dt/2.0
        k1 = f(u[n], t[n])
        k2 = f(u[n] + dt2*k1, t[n] + dt2)
        unew = u[n] + dt*k2
        return unew

class RungeKutta4(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        dt2 = dt/2.0
        k1 = f(u[n], t[n])
        k2 = f(u[n] + dt2*k1, t[n] + dt2)
        k3 = f(u[n] + dt2*k2, t[n] + dt2)
        k4 = f(u[n] + dt*k3, t[n] + dt)
        unew = u[n] + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        return unew


registered_solver_classes = [
    ForwardEuler, ExplicitMidpoint, RungeKutta4]

def test_exact_numerical_solution():
    """
    Test the different methods for a problem
    where the analytical solution is known and linear.
    All the methods should be exact to machine precision
    for this choice.
    """
    a = 0.2; b = 3

    def f(u, t):
        return a + (u - u_exact(t))**5

    def u_exact(t):
        """Exact u(t) corresponding to f above."""
        return a*t + b

    U0 = u_exact(0)
    T = 8
    N = 10
    tol = 1E-15
    t_points = np.linspace(0, T, N)
    for solver_class in registered_solver_classes:
        solver = solver_class(f)
        solver.set_initial_condition(U0)
        u, t = solver.solve(t_points)
        u_e = u_exact(t)
        max_error = (u_e - u).max()
        msg = f'{solver.__class__.__name__} failed with max_error={max_error}'
        assert max_error < tol, msg


if __name__ == '__main__':
    test_exact_numerical_solution()



class DifferentialEquation:
    def __init__(self, f):
        self.F = f
        self.CurrentT = 0
    def __call__(self, u, t):
        self.CurrentT = t
        return self.F(u)

    

deq1 = DifferentialEquation(lambda u : u * (1.5/N))


# Benchmark analytic solution
def analytic_solution(t):
    return  0.1 * np.exp(0.2 *t)



def runAnPlot():
    fwEuler1 = ForwardEuler_v2(deq)
    fwEuler1.set_initial_condition(InitialCondition)
    
    fwEuler2 = ForwardEuler_v2(deq)
    fwEuler2.set_initial_condition(InitialCondition)
    
    fwEuler3 = ForwardEuler_v2(deq)
    fwEuler3.set_initial_condition(InitialCondition)
    
    tseries = np.linspace(MinT, MaxT, N)
    print(tseries)
    u_1, t_1 = fwEuler1.solve(np.linspace(MinT,MaxT, N))
    u_2, t_2 = fwEuler2.solve(np.linspace(MinT,MaxT, N*10))
    u_3, t_3 = fwEuler3.solve(np.linspace(MinT,MaxT, N*100))
    
    
    
    #plt.plot(tseries, analytic_solution(tseries), 'r', label='Analytic') 
    plt.plot(t_1, u_1,  label=f"Standard Euler(N = {N})")
    plt.plot(t_2, u_2,  label=f"Standard Euler(N = {N*10}")
    plt.plot(t_3, u_3,  label=f"Standard Euler(N = {N*100})")
    
    plt.title('Eulers method')
    
    plt.xlabel('x')
    plt.ylabel('y(x)') 
    plt.legend()
    plt.savefig('simple_ode_class_ODESolver.jpg')
    plt.show()
    
runAnPlot()

