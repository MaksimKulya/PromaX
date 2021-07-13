import numpy as np
import math

pi=math.pi

def THz(N,T,tau):
    
    # t = np.zeros(N)
    # E = np.zeros(N)
    # nu = np.zeros(N)
    # G = np.zeros(N)
      
    # for i in range(N):
    #     t[i] = (i - N / 2) * T / N
    #     E[i] = (t[i] / tau) * math.exp(-(t[i] ** 2) / (tau ** 2))
    #     nu[i] = (1 / T) * (i - N / 2)
      
    t = np.linspace(-T/2, T/2-T/N, N)  # create time grid
    nu = np.linspace(-N/2/T, N/2/T-1/T, N)  # create frequency grid
    E = np.exp(-(t / tau) ** 2) * t / tau  # initialize E(t)
    
    E = E / np.max(E)
    
    G = np.fft.fftshift(np.fft.fft(np.fft.fftshift(E)))
      
      
    # sumE = np.sum(np.square(E))
    # sumG = np.sum(np.square(np.abs(G)))/N
      
    return E,t,G,nu


def THz_noise(N,T,tau,noisepow):
    
    # t = np.zeros(N)
    # E = np.zeros(N)
    # nu = np.zeros(N)
    # G = np.zeros(N)
      
    # for i in range(N):
    #     t[i] = (i - N / 2) * T / N
    #     E[i] = (t[i] / tau) * math.exp(-(t[i] ** 2) / (tau ** 2))
    #     nu[i] = (1 / T) * (i - N / 2)
      
    t = np.linspace(-T/2, T/2-T/N, N)  # create time grid
    nu = np.linspace(-N/2/T, N/2/T-1/T, N)  # create frequency grid
    E = np.exp(-(t / tau) ** 2) * t / tau  # initialize E(t)
    
    E = E / np.max(E)
    
    noise = np.random.normal(0,noisepow,len(E))
    E = E + noise
    
    G = np.fft.fftshift(np.fft.fft(np.fft.fftshift(E)))
      
      
    # sumE = np.sum(np.square(E))
    # sumG = np.sum(np.square(np.abs(G)))/N
      
    return E,t,G,nu


def THz1_5(N,T,tau):
    
    # t = np.zeros(N)
    # E = np.zeros(N)
    # nu = np.zeros(N)
    # G = np.zeros(N)
    
    # for i in range(N):
    #     t[i] = (i - N / 2) * T / N
    #     E[i] = (t[i] / tau) * math.exp(-(t[i] ** 2) / (tau ** 2))
    #     nu[i] = (1 / T) * (i - N / 2)
    
    t = np.linspace(-T/2, T/2-T/N, N)  # create time grid
    nu = np.linspace(-N/2/T, N/2/T-1/T, N)  # create frequency grid
    E = np.exp(-(t / tau) ** 2) * np.cos(4.5*t/tau)  # initialize E(t)
    
    E = E / np.max(E)
    G = np.fft.fftshift(np.fft.fft(np.fft.fftshift(E)))
    
    
    # sumE = np.sum(np.square(E))
    # sumG = np.sum(np.square(np.abs(G)))/N
    
    return E,t,G,nu


def femto(N,T,tau,nuind):

    # t = np.zeros(N)
    # E = np.zeros(N)
    # nu = np.zeros(N)
    # G = np.zeros(N)
    
    # for i in range(N):
    #     t[i] = (i - N / 2) * T / N
    #     E[i] = (t[i] / tau) * math.exp(-(t[i] ** 2) / (tau ** 2))
    #     nu[i] = (1 / T) * (i - N / 2)
    
    t = np.linspace(-T/2, T/2-T/N, N)  # create time grid
    nu = np.linspace(-N/2/T, N/2/T-1/T, N)  # create frequency grid
    E = np.exp(-2*(t/tau)**2)*np.sin(2*pi*nu[nuind]*t)  # initialize E(t)
    
    E = E / np.max(E)
    G = np.fft.fftshift(np.fft.fft(np.fft.fftshift(E)))
    
    
    # sumE = np.sum(np.square(E))
    # sumG = np.sum(np.square(np.abs(G)))/N
    
    return E,t,G,nu





