# notas:
# - (!) acho que isto so era interessante se pegasse nos best-fits de um resultado com dados reais qqrs, meter o erro esperado da LISA e mostrar o plot final que isto da
# - hardcoded para 4 processos e para floating point precision de 64 bit, pq, vamos ser sinceros, mais ngn vai usar isto sem seu eu no meu pc xD
# - na grid ele considera as entradas da list como o ponto do meio, nao faz diference se tiver mtos pontos
# - nao me consigo livrar do erro da obserrvacao, pq se tivermos uma observacao com um erro gigante, eh mais provavel encontra-la longe do que se tivesse um erro pequeno

# to-do:
# - (!) arranjar forma de sacar a posterior distribution (numerica) do MCMC
# - (!) assegurar que a probabilidade total esta normalizada a 1
# - (!) implementar o error em funcao do redshift esperado para um dado observatorio pq senao parece que redshifts distantes sao improvaveis, mas nao, so estao dispersos
# - meter uma loading bar (e estimated time)

# imports
from multiprocessing.shared_memory import SharedMemory
from scipy.integrate import dblquad, quad
from multiprocessing import Pool
import matplotlib.pyplot as plt
from datetime import timedelta
from numpy import linspace
from math import pi, exp
import numpy as np
import time

# gaussian
# (cuz nobody seems to have coded it yet?!)
def gauss(x, mu, sigma):
    return 1/((2*pi)**0.5 * sigma) * exp(-0.5*((x-mu)/sigma)**2)

# model
def model(z, h, Om):
    return (1+z) * (2.9979/h) * quad(lambda Z: 1/((Om*(1+Z)**3 + 1-Om)**0.5), 0, z)[0]

# likelihood for a single event
def likelihood(dLobs, h, Om, zobs, error):
    dLtheo = model(zobs, h, Om)
    return gauss(dLobs, dLtheo, error)

# posterior distribution
def posterior(h, Om):
    ph = gauss(h, 0.7, 0.0058)    # como eh que raio vou sacar a dist posterior do MCMC e meter aqui?
    pOm = gauss(Om, 0.284, 0.017) # .
    return ph*pOm

# loop that will be distributed between processes
def loop(zstart, zend, dLstart, dLend, N, buffer_name):
    # get shared memory and consider it as a numpy array
    shm = SharedMemory(name=buffer_name)
    pgrid = np.ndarray((N, N), dtype="f8", buffer=shm.buf)

    for i in range(zstart, zend):
        for j in range(dLstart, dLend):
            zobs = zgrid[i]
            dLobs = dLgrid[j]
            error = 0.5  # erro da observacao hipotetica

            # note: "swaping" x and y because of how matplotlib shows the results
            pgrid[j][i] = dblquad(lambda Om,h: likelihood(dLobs, h, Om, zobs, error)*posterior(h, Om), hmin, hmax, Ommin, Ommax)[0]

    # close connection to the buffer for this process
    shm.close()

    pass

# main
def main():
    # set read-only variables as global to avoid passing them directly to the children
    global zgrid
    global dLgrid
    global hmin
    global hmax
    global Ommin
    global Ommax

    # bounds for the parameters
    # (de futuro tem q vir diretamente da posteriori do MCMC)
    hmin = 0.7 - 3*0.0058
    hmax = 0.7 + 3*0.0058
    Ommin = 0.284 - 3*0.017
    Ommax = 0.284 + 3*0.017

    # initialize the probability grid
    # N must be multiple of 2
    N = 32
    zmin = 0
    zmax = 2
    dLmin = 0
    dLmax = 17.5
    zgrid = linspace(zmin, zmax, N).tolist()
    dLgrid = linspace(dLmin, dLmax, N).tolist()

    # create shared memory to share between child processes
    # hardcoded for 64 bits (8 bytes) floating numbers
    shm = SharedMemory(create=True, size=N*N*(64//8))
    pgrid = np.ndarray((N, N), dtype="f8", buffer=shm.buf)

    timestart = time.time()

    # compute the probability grid in parallel
    processes=4
    with Pool(processes=processes) as pool:
        args = [(0, N//2, 0, N//2, N, shm.name),
                (0, N//2, N//2, N, N, shm.name),
                (N//2, N, 0, N//2, N, shm.name),
                (N//2, N, N//2, N, N, shm.name)]
        pool.starmap(loop, args)

    timeend = time.time()
    timeelapsed = timedelta(seconds = round(timeend - timestart))
    print(f"Grid calculation time elapsed (h:m:s): {timeelapsed}")

    # plot probability grid
    fig, ax = plt.subplots()
    pcm = ax.pcolormesh(zgrid, dLgrid, pgrid)
    plt.colorbar(pcm, ax=ax)

    # plot best fit
    line = np.linspace(zmin, zmax, 100).tolist()
    plt.plot(line, [model(i, 0.7, 0.284) for i in line], label="Best fit", color="red")

    # configure and show plot
    plt.xlabel("z")
    plt.ylabel("$d_L(z)$")
    plt.xlim((zmin, zmax))
    plt.ylim((dLmin, dLmax))
    plt.grid(alpha=0.25)
    plt.legend()
    plt.show()

    # free and release shared memory block for all processes
    shm.close()
    shm.unlink()

# run
if __name__ == "__main__":
    main()
