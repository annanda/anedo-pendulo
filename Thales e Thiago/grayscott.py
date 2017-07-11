#!/usr/bin/env python3

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation

EULER = 1
RK2 = 2
RK4 = 4

method = RK2

N = 200
Nrows = N
Ncols = N
Da = 1.
Db = .5

# Mitose
# f = .0367
# k = .0649

# Corais
# f = .0545
# k = .0620

# Ondas
f = .014
k = .045

# Padrão maneiro com bolinhas
# Db = .3
# f = .035
# k = .065

# Padrão maneiro em U
# f = np.transpose(np.array([list(np.linspace(0.01, 0.10, Nrows)) for _ in range(Ncols)]))
# k = np.array([list(np.linspace(0.045, 0.07, Ncols)) for _ in range(Nrows)])

# Outro padrão maneiro
# Db = .3
# f = .078
# k = .061

lap = np.array([[.05, .20, .05],
                [.20, -1., .20],
                [.05, .20, .05]])

def laplacian(M):
	return signal.convolve2d(M, lap, mode='same', boundary='fill')

def fA(A, B, ABB):
	lapA = laplacian(A)
	return Da*lapA - ABB + f*(1-A)

def fB(A, B, ABB):
	lapB = laplacian(B)
	return Db*lapB + ABB - (k+f)*B

def grayscott(A, B, dt):

	ABB = A*B*B
	A1 = fA(A, B, ABB)
	B1 = fB(A, B, ABB)

	# Euler explicito
	if method is EULER:
		return A + A1*dt, B + B1*dt

	# Runge-Kutta segunda ordem
	if method is RK2:
		Aaux = A + A1*dt
		Baux = B + B1*dt
		ABB = Aaux*Baux*Baux
		A2 = fA(Aaux, Baux, ABB)
		B2 = fB(Aaux, Baux, ABB)

		return A + (A1 + A2)*dt/2, B + (B1 + B2)*dt/2

	Aaux = A + A1*dt/2
	Baux = B + B1*dt/2
	ABB = Aaux*Baux*Baux
	A2 = fA(Aaux, Baux, ABB)
	B2 = fB(Aaux, Baux, ABB)

	Aaux = A + A2*dt/2
	Baux = B + B2*dt/2
	ABB = Aaux*Baux*Baux
	A3 = fA(Aaux, Baux, ABB)
	B3 = fB(Aaux, Baux, ABB)
	
	Aaux = A + A3*dt
	Baux = B + B3*dt
	ABB = Aaux*Baux*Baux
	A4 = fA(Aaux, Baux, ABB)
	B4 = fB(Aaux, Baux, ABB)

	# Runge-Kutta quarta ordem
	A_ = A + (A1 + 2*A2 + 2*A3 + A4)*dt/6
	B_ = B + (B1 + 2*B2 + 2*B3 + B4)*dt/6

	return A_, B_

A, B = None, None

def run(*args):
	global A, B
	dt = 1.
	for _ in range(10):
		A, B = grayscott(A, B, dt)
	im.set_array(B)

	return im,

if __name__ == '__main__':

	A = np.ones((Nrows, Ncols))
	B = np.zeros((Nrows, Ncols))

	rngRows = range(3*Nrows//9, 6*Nrows//9)
	rngCols = range(3*Ncols//9, 6*Ncols//9)
	for i in rngRows:
		for j in rngCols:
			B[i, j] = .95

	fig = plt.figure()
	im = plt.imshow(B, vmax=0.5)

	ani = animation.FuncAnimation(fig, run, blit=True, interval=1, repeat=False)
	plt.show()
