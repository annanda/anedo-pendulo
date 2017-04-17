# -*- coding: utf-8 -*-

import numpy as np
import sys


class EDOPendulo(object):
    def __init__(self, t_inicial, t_final, n_intervalos, v_inicial, theta_inicial=np.pi/24,  g=10, L=1):
        self.g = g
        self.theta_inicial = theta_inicial
        self.L = L
        self.t_inicial = t_inicial
        self.t_final = t_final
        self.v_inicial = v_inicial
        self.n_intervalos = n_intervalos
        self.atualizar_h()
        
    def atualizar_h(self):
        self.h = (self.t_final - self.t_inicial)/float(self.n_intervalos)        
        
    def euler_explicito_nao_linear(self):
        theta = [self.theta_inicial]
        v = [self.v_inicial]
        t = [self.t_inicial]
        
        for i in xrange(0, self.n_intervalos):
            theta_next = theta[i]+self.h * v[i]
            v_next = v[i] - self.h * ((self.g/self.L) * np.sin(theta_next))
            t_next = t[i] + self.h
            
            theta.append(theta_next)
            v.append(v_next)
            t.append(t_next)
        
        return np.array(t), np.array(theta), np.array(v)
    
    def euler_implicito_linear(self):
        theta = [self.theta_inicial]
        v = [self.v_inicial]
        t = [self.t_inicial]
        
        for i in xrange(0, self.n_intervalos):
            
            v_next = (v[i] - (self.h * (self.g/self.L) * theta[i])) / (1. + (self.h**2) * self.g/self.L)
            theta_next = theta[i] + self.h * v_next
            t_next = t[i] + self.h
            
            theta.append(theta_next)
            v.append(v_next)
            t.append(t_next)        
        
        return np.array(t), np.array(theta), np.array(v)
    
    def euler_implicito_nao_linear(self):
        theta = [self.theta_inicial]
        v = [self.v_inicial]
        t = [self.t_inicial]
        
        for i in xrange(0, self.n_intervalos):
            
            v_next = self.newton(theta[i], v[i])
            theta_next = theta[i] + self.h * v_next
            t_next = t[i] + self.h
            
            theta.append(theta_next)
            v.append(v_next)
            t.append(t_next)        
        
        return np.array(t), np.array(theta), np.array(v)
    
    def f(self, v_next, theta, v):
        return v_next - (v - (self.h * (self.g/self.L)) * np.sin(theta + self.h * v_next))
    
    def df(self, v_next, theta, v):
        return self.g * (self.h**2) * np.cos(self.h * v_next + theta)/self.L + 1.
    
    def newton(self, theta, v):
        x = 10. #initial guess
        
        i=0
        while (i<100): #100 iterations
            x_next = x-self.f(x, theta, v) / self.df(x, theta, v)
            x = x_next
            i+=1        
        
        return x

if __name__ == '__main__':
    theta = np.pi/24
    v = 0.
    g = 10
    L = 1
    t_inicial = 0.
    t_final = 10.
    n_intervalos = 100

    edo = EDOPendulo(t_inicial, t_final, n_intervalos, v, theta, g, L)
    t_explNL, p_explNL, v_explNL = edo.euler_implicito_linear()
    t_implL, p_implL, v_implL = edo.euler_implicito_linear()
    t_implNL, p_implNL, v_implNL = edo.euler_implicito_nao_linear()

    error = np.absolute(p_implNL - p_implL)

    print "Erro em theta: "
    print (np.amax(error))

    error = np.absolute(v_implNL - v_implL)
    print "Erro em v: "
    print (np.amax(error))