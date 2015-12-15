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
for line in lines:
    p = line.split()
    totalEnergy.append(float(p[0])+restMass)
    momentum.append(np.sqrt(float(p[0])*2*restMass))
    eta.append((2*np.pi*fine*(float(p[0])+restMass))/np.sqrt(float(p[0])*2*restMass))
totalEnergy = np.array(totalEnergy)
momentum = np.array(momentum)
eta = np.array(eta)

def square(list):
    return map(lambda x: x ** 2, list)


def f(zp, b):                                               #for b 0 = + and 1 = -

    f = ((-1)**b * zp * eta)/(np.exp((-1)**b * zp * eta)-1)
    return f
y1 = np.sqrt(counts/(momentum**2 * f(20,0)))
y2 = np.sqrt(counts/(momentum**2 * f(19,1)))

print y2

fig, ax = plt.subplots()
# pt0 = ax.errorbar(energy,np.sqrt(counts/(((-1)**0 * 19 * eta)/(np.exp((-1)**0 * 19 * eta)-1)*momentum*momentum)), yerr=yerr, fmt='.', ecolor='g', capthick=.2)

# pt0 = ax.errorbar(energy,y1)


# print np.sqrt(counts/(f(19,0)*momentum*momentum))

plt.show()
