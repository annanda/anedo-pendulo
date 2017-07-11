#!/usr/bin/env python3

from string import Template
import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit, pycuda.compiler
import matplotlib.pyplot as plt
import matplotlib.animation as animation

M = 500    # Mat. height
N = 500    # Mat. width
D_A = 1.0  # A diff. coef.
D_B = 0.5  # B diff. coef.
F = .037   # Feed rate
K = .060   # Kill rate
ITERS = 50 # Iters. per frame

template = """
    #define D_A $D_A
    #define D_B $D_B
    #define F $F
    #define K $K

    __global__ void laplacian(
        float *lap_a,
        float *a,
        int m,
        int n)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x + 1;
        int j = blockDim.y * blockIdx.y + threadIdx.y + 1;

        if (i >= m-1 || j >= n-1) return;

        float *a_i;
        float  val = 0.0f;

        a_i  = a + n*i - n;
        val += 0.05f * a_i[j-1] + 0.20f * a_i[j] + 0.05f * a_i[j+1];
        a_i += n;
        val += 0.20f * a_i[j-1] - 1.00f * a_i[j] + 0.20f * a_i[j+1];
        a_i += n;
        val += 0.05f * a_i[j-1] + 0.20f * a_i[j] + 0.05f * a_i[j+1];

        lap_a[n*i + j] = val;
    }

    __global__ void grayscott(
        float *a,
        float *b,
        float *lap_a,
        float *lap_b,
        float dt,
        int m,
        int n)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x + 1;
        int j = blockDim.y * blockIdx.y + threadIdx.y + 1;

        if (i >= m-1 || j >= n-1) return;

        // float F = 0.010 + (0.10 - 0.010) * i / m;
        // float K = 0.045 + (0.07 - 0.045) * j / n;

        float a_ij = a[n*i + j];
        float b_ij = b[n*i + j];
        float lap_a_ij = lap_a[n*i + j];
        float lap_b_ij = lap_b[n*i + j];

        float abb_ij = a_ij * b_ij * b_ij;

        a_ij += (D_A * lap_a_ij - abb_ij + F * (1 - a_ij)) * dt;
        b_ij += (D_B * lap_b_ij + abb_ij - (K + F) * b_ij) * dt;

        a[n*i + j] = min(a_ij, 1.0f);
        b[n*i + j] = min(b_ij, 1.0f);
    }
    """
code = Template(template).substitute(D_A=D_A, D_B=D_B, F=F, K=K)
mod = pycuda.compiler.SourceModule(code)

_laplacian = mod.get_function("laplacian");
_grayscott = mod.get_function("grayscott")

def laplacian(lap_a, a, m, n):
    BX = 1
    BY = 64

    gx = (m-2 + BX - 1) // BX
    gy = (n-2 + BY - 1) // BY
    _laplacian(lap_a, a, np.int32(m), np.int32(n),
        block=(BX, BY, 1), grid=(gx, gy, 1))

def grayscott(a, b, lap_a, lap_b, dt, m, n):
    laplacian(lap_a, a, m, n)
    laplacian(lap_b, b, m, n)

    BX = 1
    BY = 64

    gx = (m-2 + BX - 1) // BX
    gy = (n-2 + BY - 1) // BY
    _grayscott(a, b, lap_a, lap_b, np.float32(dt), np.int32(m), np.int32(n),
        block=(BX, BY, 1), grid=(gx, gy, 1))

'''
tmp = None
a = None
b = None
lap_a = None
lap_b = None
'''

def run(*args):
    for _ in range(ITERS):
        grayscott(a, b, lap_a, lap_b, 1.0, M+2, N+2)
    cuda.memcpy_dtoh(tmp, b)
    im.set_array(tmp)
    return im,

if __name__ == '__main__':
    tmp = np.empty((M+2, N+2), np.float32)
    a = cuda.mem_alloc(tmp.nbytes)
    b = cuda.mem_alloc(tmp.nbytes)
    lap_a = cuda.mem_alloc(tmp.nbytes)
    lap_b = cuda.mem_alloc(tmp.nbytes)

    tmp[...] = 0.0
    for i in range(1, M-1):
        for j in range(1, N-1):
            tmp[i, j] = 1.0
    cuda.memcpy_htod(a, tmp)

    tmp[...] = 0.0
    for i in range((M+2)//2 - 10, (M+2)//2 + 10):
        for j in range((N+2)//2 - 10, (N+2)//2 + 10):
            tmp[i, j] = 1.0
    cuda.memcpy_htod(b, tmp)

    fig = plt.figure()
    im = plt.imshow(tmp, vmin=0.0, vmax=0.5)
    ani = animation.FuncAnimation(fig, run, blit=True, interval=1, repeat=False)
    plt.show()