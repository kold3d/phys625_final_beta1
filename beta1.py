import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chisquare

fine = 1./137.036

data = open("phys625-betadecay-prob1.spe", 'r')
lines = data.readlines()[1:]

energy = []
counts = []

for line in lines:
    p = line.split()
    energy.append(float(p[0]))
    counts.append(float(p[1]))
energy = np.array(energy)
counts = np.array(counts)


restMass = 510.998902
totalEnergy = []
momentum = []
eta = []

totalEnergy = np.array(energy+restMass)
momentum = np.array(np.sqrt(energy * (energy + 2. * restMass)))
eta = np.array((2. * np.pi * fine * totalEnergy) / momentum)

def f(zp, b):                                               #for b 0 = + and 1 = -
    f = ((-1.) ** b * zp * 2. * np.pi * fine * totalEnergy / momentum) / (
    np.exp((-1.) ** b * zp * 2. * np.pi * fine * totalEnergy / momentum) - 1.)
    return f


y1 = np.array(np.sqrt(counts / (momentum ** 2. * f(21., 0.))))
y2 = np.array(np.sqrt(counts / (momentum ** 2. * f(19., 1.))))

# def error(y):
#     yerr = []
#     for value in y:
#         if value > 0.:
#             yerr.append(1/np.sqrt(value))
#         else:
#             yerr.append(value)
#     return yerr
#
# print max(error(counts))

yerr1 = 1. / (2. * momentum * np.sqrt(f(21., 0.)))
yerr2 = 1. / (2. * momentum * np.sqrt(f(19., 1.)))

linearFunc = lambda x, *p: p[0] * x + p[1]

popt1, pcov1 = curve_fit(linearFunc, xdata=energy, ydata=y1, sigma=yerr1, p0=[-1, .12])
popt2, pcov2 = curve_fit(linearFunc, xdata=energy, ydata=y2, sigma=yerr2, p0=[-1, .01])

print popt2, pcov2

xxfit = np.linspace(min(energy), max(energy), len(energy))
yyfit1 = linearFunc(xxfit, *popt1)
yyfit2 = linearFunc(xxfit, *popt2)

residuals1 = (y1 - linearFunc(xxfit, *popt1)) / yerr1
residuals2 = (y1 - linearFunc(xxfit, *popt2)) / yerr2

### Chi square and p-value ###
chisq1, p1 = chisquare(f_obs=y1, f_exp=linearFunc(xxfit, *popt1), ddof=len(energy) - 2)
chisq2, p2 = chisquare(f_obs=y2, f_exp=linearFunc(xxfit, *popt2), ddof=len(energy) - 2)

print(
'chisq beta-: %f, p: %f, chisq beta+: %f, p: %f' % (chisq1 / (len(energy) - 2), p1, chisq2 / (len(energy) - 2), p2))

### Figures ###
gs = gridspec.GridSpec(2, 2, height_ratios=[3.5, 1])

ax1 = plt.subplot(gs[0])  # Kurie Plot1
ax2 = plt.subplot(gs[0])  # Fit1
ax3 = plt.subplot(gs[2])  # residuals1

ax4 = plt.subplot(gs[1])  # Kurie Plot2
ax5 = plt.subplot(gs[1])  # Fit2
ax6 = plt.subplot(gs[3])  # residuals2

pt1 = ax1.errorbar(xxfit, y1, yerr=yerr1, fmt='.', ecolor='g', capthick=.2)
pt2 = ax2.plot(xxfit, yyfit1, 'r-')
pt3 = ax3.plot(xxfit, residuals1, 'g-')

pt4 = ax4.errorbar(xxfit, y2, yerr=yerr2, fmt='.', ecolor='g', capthick=.2)
pt5 = ax5.plot(xxfit, yyfit2, 'r-')
pt6 = ax6.plot(xxfit, residuals2, 'g-')

plt.show()
