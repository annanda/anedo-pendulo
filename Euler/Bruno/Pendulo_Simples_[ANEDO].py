'''BRUNO SCARAMUZZA DOS REIS'''

''' Resolução do problema do pêndulo simples linear e não-linear
através do método de passo único - Euler Explícito, Semi-explícito e
Implícito - para análise de convergência'''

''' O desenvolvimento se dará da seguinte forma:
    - Primeiro analisaremos o problema linearizado (sen(p) = p);
        - Resolução por Euler Explícito, Euler Semi-Explícito e Euler
            Implícito, com dois intervalos diferentes (h1 e h2);
        - Para cada método (e intervalo) será plotado os gráficos p x t e
            p x v;

    - Segundo analisaremos o problema não-linear;
        - Resolução por Euler Explícito, Euler Semi-Explícito e Euler
            Implícito, com dois intervalos diferentes (h1 e h2);
        - Para cada método (e intervalo) será plotado os gráficos p x t e
            p x v;
'''

import matplotlib.pyplot as plt
import numpy as np

''' CONSTANTES '''

p_init = np.pi/24;
v_init = 0;
g = -10;
L = 1;
N = 500;
t_init = 0;
t_fim = 10;
h1 = (t_fim - t_init)/N;
h2 = (t_fim - t_init)/(2*N);
tempos1 = np.linspace(t_init, t_fim, N);
tempos2 = np.linspace(t_init, t_fim, (2*N));
                     
# Todos os métodos, ao serem utilizados, foram abreviados para agilizar o processo de análise
# Linear Explícito -> lee
# Linear Semi-Explícito com P(k+1) usado em V(k+1) -> lesp
# Linear Semi-Explícito com V(k+1) usado em P(k+1) -> lesv
# Linear Implícito -> lei
# Não-Linear Explícito -> nlee
# Não-Linear Semi-Explícito com P(k+1) usado em V(k+1) -> nlesp
# Não-Linear Semi-Explícito com V(k+1) usado em P(k+1) -> nlesv
# Não-Linear Implícito -> nlei

# No caso, coloquei tempos1 e tempos2, h1 e h2 para analisar a convergência ao aumentar o
# número de nós
        
''' PROBLEMA LINEARIZADO '''

def linear_euler_exp(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        p_next = p_list[i-1] + h*v_list[i-1]
        v_next = v_list[i-1] + h*(g/L)*p_list[i-1]
        
        p_list.append(p_next)
        v_list.append(v_next)
    
    return [p_list,v_list]

lee1 = linear_euler_exp(p_init,v_init,g,L,N,h1)
lee2 = linear_euler_exp(p_init,v_init,g,L,2*N,h2)

def linear_euler_semip(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        p_next = p_list[i-1] + h*v_list[i-1]
        p_list.append(p_next)
        v_next = v_list[i-1] + h*(g/L)*p_list[i-1]
        
        v_list.append(v_next)
    
    return [p_list,v_list]

lesp1 = linear_euler_semip(p_init,v_init,g,L,N,h1)
lesp2 = linear_euler_semip(p_init,v_init,g,L,2*N,h2)

# Neste método semi-explícito temos o uso de p(k+1) para obter v(k+1)

def linear_euler_semiv(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        v_next = v_list[i-1] + h*(g/L)*p_list[i-1]
        v_list.append(v_next)
        p_next = p_list[i-1] + h*v_list[i-1]
        p_list.append(p_next)
    
    return [p_list,v_list]

lesv1 = linear_euler_semiv(p_init,v_init,g,L,N,h1)
lesv2 = linear_euler_semiv(p_init,v_init,g,L,2*N,h2)

# Neste método usa-se v(k+1) para obter p(k+1)

def linear_euler_imp(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        p_next = (p_list[i-1] + h*v_list[i-1])/(1 - ((h**2)*(g/L)))
        p_list.append(p_next)
        v_next = v_list[i-1] + h*(g/L)*p_list[i-1]
        
        v_list.append(v_next)
    
    return [p_list,v_list]

lei1 = linear_euler_imp(p_init,v_init,g,L,N,h1)
lei2 = linear_euler_imp(p_init,v_init,g,L,2*N,h2)



''' PROBLEMA NÃO LINEAR '''

def nonlinear_euler_exp(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        p_next = p_list[i-1] + h*v_list[i-1]
        v_next = v_list[i-1] + h*(g/L)*np.sin(p_list[i-1])
        
        p_list.append(p_next)
        v_list.append(v_next)
    
    return [p_list,v_list]

nlee1 = nonlinear_euler_exp(p_init,v_init,g,L,N,h1)
nlee2 = nonlinear_euler_exp(p_init,v_init,g,L,2*N,h2)

def nonlinear_euler_semip(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        p_next = p_list[i-1] + h*v_list[i-1]
        p_list.append(p_next)
        v_next = v_list[i-1] + h*(g/L)*np.sin(p_list[i-1])
        
        v_list.append(v_next)
    
    return [p_list,v_list]

nlesp1 = nonlinear_euler_semip(p_init,v_init,g,L,N,h1)
nlesp2 = nonlinear_euler_semip(p_init,v_init,g,L,2*N,h2)

def nonlinear_euler_semiv(p,v,g,L,N,h):
    p_list = [p]
    v_list = [v]
    
    for i in range(1,N):
        v_next = v_list[i-1] + h*(g/L)*np.sin(p_list[i-1])
        v_list.append(v_next)        
        p_next = p_list[i-1] + h*v_list[i-1]
        
        p_list.append(p_next)
    
    return [p_list,v_list]

nlesv1 = nonlinear_euler_semiv(p_init,v_init,g,L,N,h1)
nlesv2 = nonlinear_euler_semiv(p_init,v_init,g,L,2*N,h2)

''' MÉTODO DE NEWTON - SERÁ USADO NA RESOLUÇÃO IMPLÍCITA DO NÃO LINEAR '''
class NonLinear_Imp(object):
    
    def __init__(self,p,v,g,L,N,h):
        self.p = p
        self.v = v
        self.g = g
        self.l = L
        self.n = N
        self.h = h
        self.p1 = 0
        self.p2 = 0
        self.p_list_newton = []
        self.p_list = []
        self.v_list = []
    
    def newton(self):
        self.p_list_newton = [self.p_list[-1]] #Chute inicial é o útlimo valor de p
        self.p2 = self.p_list_newton[-1] - (self.p_list_newton[-1] - (self.h**2)*(self.g/self.l) * (np.sin(self.p_list_newton[-1])) - (self.p_list[-1] + self.h*self.v_list[-1]))/(1 - ((self.h**2)*(self.g/self.l) * (np.cos(self.p_list_newton[-1]))))
        self.p_list_newton.append(self.p2)
        i = 0
        
        while np.absolute(self.p_list_newton[-2] - self.p_list_newton[-1]) > 10**(-5):
            i = self.p_list_newton[-1] - (self.p_list_newton[-1] - (self.h**2)*(self.g/self.l) * (np.sin(self.p_list_newton[-1])) - (self.p_list[-1] + self.h*self.v_list[-1]))/(1 - ((self.h**2)*(self.g/self.l) * (np.cos(self.p_list_newton[-1]))))
            self.p_list_newton.append(i)
            
        return self.p_list_newton[-1]
    
    def nonlinear_euler_imp(self):
        self.p_list = [self.p]
        self.v_list = [self.v]
        
        for i in range(1,self.n):
            p_next = self.newton()
            self.p_list.append(p_next)
            
            v_next = self.v_list[i-1] + self.h * (self.g/self.l) * (np.sin(self.p_list[i-1]))
            self.v_list.append(v_next)
        
        return [self.p_list,self.v_list]

NonLinear1 = NonLinear_Imp(p_init,v_init,g,L,N,h1)
nlei1 = NonLinear1.nonlinear_euler_imp()
NonLinear2 = NonLinear_Imp(p_init,v_init,g,L,2*N,h2)
nlei2 = NonLinear2.nonlinear_euler_imp()
