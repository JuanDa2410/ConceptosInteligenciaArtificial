# Función para leer el archivo de instancia.
import csv
import random
import copy
import numpy as np
import matplotlib.pyplot as plt

def importfile(filename):
    with open(filename) as filename:
        rows = csv.reader(filename)
        for i, row in enumerate(rows): 
            if i == 0:
                data = []
                data.append([i] + row[1:4])
            else: 
                #print(row)
                data.append([i] + row[1:4]) 
    return data

def setinitialprobability(data, knapsackcapacity):
    sumvalueperweight = 0
    for i in data:
        sumvalueperweight = sumvalueperweight + (int(i[1])/int(i[2]))
    res = []
    for i in data:
        res.append([i[0]]+ [(int(i[1])/(int(i[2])*sumvalueperweight))])
    return res

def updateneighbour(previousneighbourhood, presentknapsackcapacity, item):
    presentneightbourhood = copy.deepcopy(previousneighbourhood)
    removelist = []
    for i in previousneighbourhood.keys(): #previousneighbourhood es un diccionario de {índice: [valor, peso]} i en el índice aquí
        if presentknapsackcapacity < int(presentneightbourhood[i][1]):
            removelist.append(i)
    if item not in removelist: 
        removelist.append(item) 
    for i in removelist:
        del presentneightbourhood[i]
    return presentneightbourhood

def generatetransitionmatrix(presentneigh, tau, mu, alpha=3, beta=2):
    expresssum = 0
    for i in presentneigh.keys():
        expresssum = expresssum + ((tau[i]**alpha)*(mu[i]**beta))

    probmatrix = {
                    i: ((tau[i]**alpha)*(mu[i]**beta))/expresssum for i in presentneigh.keys()
                 }
    return probmatrix

def sampleitem(probmatrix):
    rand_val = random.random()
    total = 0
    for k, v in probmatrix.items():
        total += v
        if rand_val <= total:
            return k
    assert False, 'inaccesible'

def updatephero(resofkants, tau, globalprofit):
    for i in resofkants:
        z = resofkants[i][0]
        deltatau = (1/(1+((globalprofit - z)/globalprofit)))
        for j in resofkants[i][1].keys():
            tau[j] = tau[j] + deltatau
    return tau

def evaporate(tau, rho):
    for i in tau.keys():
        tau[i] = max(0.05, tau[i]*rho)
    return tau

def normalize(probmat):
    summ = 0
    for i in probmat:
        summ = summ + probmat[i]
    for i in probmat:
        probmat[i] = probmat[i]/summ
    return probmat

def aco(iterations, neighbourhood, knapsackcapacity, numants, probmatrix, tau, mu, alpha, beta):
    #random.seed(1)   SEED FUNCIONANDO
    #seed = seed + 1
    iternumber = 0
    globalprofit = -1
    globaltrack = []
    globalsample = {}

    while iternumber < iterations:
        presentneigh = neighbourhood
        presentknapsackcapacity = knapsackcapacity
        resofkants = {}
        #print("Número de iteración: ", iternumber)
        antprofit = -1
        antsample = {}

        for ant in range(numants):
            localprofit = 0
            solutionset = {}
            while (presentknapsackcapacity > 0) and (len(presentneigh.keys()) > 0):
                presentkeys = presentneigh.keys()
                #print("Longitud de los vecinos presentes= ", len(presentneigh))
                probmat = dict((k, probmatrix[k]) for k in presentkeys if k in probmatrix)
                probmat = normalize(probmat)

                item = sampleitem(probmat)
                #print("Elemento colocado en la mochila: ", item)
                solutionset.update({item: presentneigh[item]})
                presentknapsackcapacity = presentknapsackcapacity - presentneigh[item][1]
                localprofit = localprofit + presentneigh[item][0]
                presentneigh = updateneighbour(presentneigh, presentknapsackcapacity, item)
                #print("Longitud de los vecinos presentes después de colocar el elemento en la mochila: ", len(presentneigh))
                #print("El beneficio local es: ", localprofit)

            if localprofit > antprofit:
                antprofit = localprofit
                antsample = solutionset
                #print("El beneficio de la hormiga se actualiza a ", antprofit)

            resofkants.update({ant: [localprofit, solutionset]})

        if antprofit > globalprofit:
            #print("El mejor beneficio global se actualiza a ", globalprofit)
            globalprofit = antprofit
            globalsample = antsample
            globaltrack.append(globalprofit)

        tau = updatephero(resofkants, tau, globalprofit)
        probmatrix = generatetransitionmatrix(neighbourhood, tau, mu, alpha, beta)
        iternumber = iternumber + 1

    print("alpha =", alpha, "beta =", beta, "Mejor global = ", globalprofit)
    return globalprofit

# Función que "LLAMA A LA FUNCIÓN ACO" y traza los contornos.
def getcountour(alpharange, betarange, binsize, iterations, neighbourhood, knapsackcapacity, numants, probmatrix, tau, mu):
    # Tasa de cruce en el eje x:
    xlist = np.linspace(alpharange[0], alpharange[1], binsize)
    # Tasa de mutación en el eje y:
    ylist = np.linspace(betarange[0], betarange[1], binsize)
    X, Y = np.meshgrid(xlist, ylist)

    fullz = []
    for i in xlist:
        a1 = []
        for j in ylist:
            res = aco(iterations, neighbourhood, knapsackcapacity, numants, probmatrix, tau, mu, i, j)
            a1.append(res)
        fullz.append(a1)

    fig = plt.figure(figsize=(6, 5))
    left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
    ax = fig.add_axes([left, bottom, width, height])

    cp = ax.contour(X, Y, fullz)
    ax.clabel(cp, inline=True, fontsize=10)
    ax.set_title('Contour_13_50_1000.csv')
    ax.set_xlabel('alpha')
    ax.set_ylabel('beta')
    plt.show()
