import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

fine = 1./137.036

data = open("phys625-betadecay-prob1.spe", 'r')
lines = data.readlines()[1:]

energy = []
counts = []

for line in lines:
    p = line.split()
    energy.append(int(p[0]))
    counts.append(int(p[1]))
energy = np.array(energy)
counts = np.array(counts)


restMass = 510.998902
totalEnergy = []
momentum = []
eta = []

totalEnergy = np.array(energy+restMass)
momentum = np.array(np.sqrt(energy*2.*restMass))
eta = np.array((2*np.pi*fine*totalEnergy)/momentum)

def f(zp, b):                                               #for b 0 = + and 1 = -

    f = ((-1)**b * zp * 2 * np.pi * fine * totalEnergy/momentum)/(np.exp((-1)**b * zp * 2 * np.pi * fine * totalEnergy/momentum)-1)
    return f

y1 = np.array(np.sqrt(counts/(momentum**2 * f(21,0))))
y2 = np.array(np.sqrt(counts/(momentum**2 * f(19,1))))

yerr = []

def error(y):
    for value in y:
        if value > 0.:
            yerr.append(1/np.sqrt(value))
        else:
            yerr.append(value)
    return yerr

print min(error(counts))


fig, ax = plt.subplots()
# pt0 = ax.errorbar(energy,np.sqrt(counts/(((-1)**0 * 19 * eta)/(np.exp((-1)**0 * 19 * eta)-1)*momentum*momentum)), yerr=yerr, fmt='.', ecolor='g', capthick=.2)

# pt0 = ax.errorbar(energy,y2, yerr=error(counts), fmt='.', ecolor='g', capthick=.2)
#pt1 = ax.errorbar(energy, y1, yerr=yerr, fmt='.', ecolor='g', capthick=.2)


# print np.sqrt(counts/(f(19,0)*momentum*momentum))

plt.show()
