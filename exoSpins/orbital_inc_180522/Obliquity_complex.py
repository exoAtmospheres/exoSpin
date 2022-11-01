import numpy as np
import matplotlib.pyplot as plt
from pyqt_fit import kde,kde_methods
import pickle

G = 6.67e-20 #km**3 kg**-1 s**-2
rjup = 71492 #km
mjup = 1.898e27 #kg
h2s = 3600 #sec
dr = np.pi/180

# filename = 'result_nestle_ABPicb_B7.pic'
#
# with open(filename, 'rb') as f:
#     res  = pickle.load(f)
#     samp = res.samples
#     weig = res.weights
# radius = samp[:,4]*rjup #km
# vsini = samp[:,6] #km/s

file  = open("result_nestle_ABPicb_B7.txt", "r")
data = np.loadtxt(file, skiprows=1)
file.close()

radius = data[:,4]*rjup #km
vsini = data[:,6] #km/s

p = 2.12*h2s #s
v = 2*np.pi*radius/p # km/s

mp = 10*mjup #kg
pb = 2*np.pi*(radius**(3/2))/np.sqrt(G*mp) # break-up period (s)

print(f"Size initial sample = {len(v)}")
sample = (p > pb) # Removing solutions above the break-up velocity
vsini = vsini[sample]
v = v[sample]

vmax = v.max()
print(f"Maximum velocity with P>Pb = {vmax} km/s")

lv = kde.KDE1D(v)
lu = kde.KDE1D(vsini)

nv = 1000
vrange = np.linspace(0, vmax, nv)
i = np.linspace(0, np.pi, 100)

pcosi = np.zeros_like(i) # PDF of cos(ip)
for j, ci in enumerate(np.cos(i)):
    int_dv = lv(vrange)*lu(vrange*np.sqrt(1-ci*ci))
    pcosi[j] = np.trapz(int_dv, vrange)

pi = pcosi*np.sin(i) # PDF of ip

int_pi = np.trapz(pi, i) # normalization factor
pi /= int_pi

sample = (vsini < v)
nsample = np.sum(sample)
print(f"Size initial sample = {len(v)}")
print(f"Size new sample = {nsample} (vsini < v < vb)")
sini = vsini/v
ip = np.arcsin(sini[sample])
ip = np.concatenate((ip, np.pi-ip))
nsample *= 2

file = open("io.dat", "r")
io = np.loadtxt(file, skiprows=1)
file.close()

lo = kde.KDE1D(io, lower=0., upper=np.pi, method=kde_methods.cyclic) # PDF of io

### Here I compute proba(l) = proba(|io-ip|) = int pi(ip)*(pio(ip+l)+pio(ip-l))/2 (combination of probability functions)
l = np.linspace(0, np.pi, 100)
pl = np.zeros_like(l)
i = np.linspace(0, np.pi, 100)

def pio(x): # PDF of io, being sure it is defined everywhere (not only [0, pi])
    if np.isscalar(x):
        if (x > 0) and (x < np.pi):
            return(lo(x))
        else:
            return(0)
    else:
        res = np.zeros_like(x)
        for k, xx in enumerate(x):
            if (xx > 0) and (xx < np.pi):
                res[k] = lo(xx)
            else:
                res[k] = 0
        return(res)

for j, la in enumerate(l):
    int_dl = pi*(pio(i-la)+pio(i+la))/2
    pl[j] = np.trapz(int_dl, i)

int_pl = np.trapz(pl, l)
pl /= int_pl

# plot projected obliquity
plt.plot(l/dr, pl*dr, label="prior cos")

plt.xlim(xmin=0)
plt.ylim(ymin=0)

plt.xlabel(r"Projected Obliquity $|i_{\rm p} - i_{\rm o}|$ (deg)")
plt.ylabel("PDF")

plt.tight_layout()
plt.show()
