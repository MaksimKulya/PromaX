from numpy import loadtxt



narr = loadtxt('D:\\SCIENCE\\2021\\RNF\\Mar\\one_layer_eff.txt', comments="#", delimiter=",", unpack=False)

narr2=narr[:,1]