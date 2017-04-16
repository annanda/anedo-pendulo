# -*- coding: utf-8 -*-

import numpy as np


class EDOPendulo(object):
    def __init__(self, g, L, t_inicial, t_final, n_intervalos):
        self.g=g
        self.L=L
        self.t_inicial=t_inicial
        self.t_final=t_final
        self.n_intervalos=n_intervalos
        self.atualizarH()
        
    def atualizarH(self):
        self.h = (self.t_final - self.t_inicial)/float(self.n_intervalos)        
        
    def eulerExplicitoNL(self, theta0, v0):
        theta = [theta0]
        v = [v0]
        t = [self.t_inicial]
        
        for i in xrange(0, self.n_intervalos):
            theta_next = theta[i]+self.h*v[i]
            v_next = v[i] - self.h*((self.g/self.L)*np.sin(theta_next))
            t_next = t[i] + self.h
            
            theta.append(theta_next)
            v.append(v_next)
            t.append(t_next)
        
        return np.array(t), np.array(theta), np.array(v)
    
    def eulerImplicitoL(self, theta0, v0):
        theta = [theta0]
        v = [v0]
        t = [self.t_inicial]
        
        for i in xrange(0, self.n_intervalos):
            
            v_next = (v[i] - (self.h*(self.g/self.L)*theta[i])) / (1. + (self.h**2)*self.g/self.L)
            theta_next = theta[i] + self.h*v_next
            t_next = t[i] + self.h
            
            theta.append(theta_next)
            v.append(v_next)
            t.append(t_next)        
        
        return np.array(t), np.array(theta), np.array(v)
    
    def eulerImplicitoNL(self, theta0, v0):
        theta = [theta0]
        v = [v0]
        t = [self.t_inicial]
        
        for i in xrange(0, self.n_intervalos):
            
            v_next = self.newton(theta[i], v[i])
            theta_next = theta[i] + self.h*v_next
            t_next = t[i] + self.h
            
            theta.append(theta_next)
            v.append(v_next)
            t.append(t_next)        
        
        return np.array(t), np.array(theta), np.array(v)
    
    def f(self, v_next, theta, v):
        return v_next - (v - (self.h*(self.g/self.L))*np.sin(theta + self.h*v_next))
    
    def df(self, v_next, theta, v):
        return self.g*(self.h**2)*np.cos(self.h*v_next + theta)/self.L + 1.
    
    def newton(self, theta, v):
        x = 10. #initial guess
        
        i=0
        while (i<100): #100 iterations
            x_next = x-self.f(x, theta, v)/self.df(x, theta, v)            
            x = x_next
            i+=1        
        
        return x


theta = np.pi/24
v = 0.
g = 10
L = 1
t_inicial = 0.
t_final = 10.
n_intervalos = 100

edo = EDOPendulo(g, L, t_inicial, t_final, n_intervalos)
t_explNL, p_explNL, v_explNL = edo.eulerImplicitoL(theta, v)
t_implL, p_implL, v_implL = edo.eulerImplicitoL(theta, v)
t_implNL, p_implNL, v_implNL = edo.eulerImplicitoNL(theta, v)

error = np.absolute(p_implNL - p_implL)

print "Erro em theta: "
print (np.amax(error))

error = np.absolute(v_implNL - v_implL)
print "Erro em v: "
print (np.amax(error))