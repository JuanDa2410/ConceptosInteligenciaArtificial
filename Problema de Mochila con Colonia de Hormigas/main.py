#%% @DAVID CASTAÑO
import csv
import funcion as func
import random
import copy
import matplotlib.pyplot as pl
import pandas as pd

#%% Obtener los datos
datafilename = "firstinstance.csv"
data = func.importfile(datafilename)

#%% Parámetros del sistema

knapsackcapacity = 970  # CAPACIDAD DE LA MOCHILA (HARD-CODED)
presentknapsackcapacity = knapsackcapacity

# initialprob = func.setinitialprobability(data, knapsackcapacity)
# neighbourhood = func.updateneighbour(neighbourhood, presentknapsackcapacity)
# Comprueba si los datos son incorrectos, es decir, si hay algún elemento con un peso mayor a la capacidad de la mochila.

iterations = 20
numants = 10
rho = 0.2

neighbourhood = {i[0]: [int(i[1]), int(i[2])] for i in data}
presentneigh = copy.deepcopy(neighbourhood)

tau = {i[0]: 10 for i in data}
mu = {i[0]: ((int(i[1]) * presentknapsackcapacity) / int(i[2])) for i in data}
alpha = 3
beta = 2

probmatrix = func.generatetransitionmatrix(presentneigh, tau, mu, alpha, beta)
presentprobmatrix = copy.deepcopy(probmatrix)

a = func.sampleitem(probmatrix)

#%% ACO (Optimización por colonias de hormigas)
func.aco(iterations, neighbourhood, knapsackcapacity, numants, probmatrix, tau, mu, alpha, beta)

#%% Gráficos de contorno para diferentes valores de alpha y beta.

alpharange = tuple((0.5, 10))
betarange = tuple((0.5, 10))
binsize = 3

func.getcountour(alpharange, betarange, binsize, iterations, neighbourhood, knapsackcapacity, numants, probmatrix, tau, mu)

#%% Diagrama de caja y resumen
ress = []
for _ in range(30):
    a = func.aco(iterations, neighbourhood, knapsackcapacity, numants, probmatrix, tau, mu, alpha, beta)
    ress.append(a)

#pl.boxplot(ress)
pdress = pd.DataFrame(ress)
summary = pdress.describe()
summary = summary.transpose()
summary.head()
#%%
