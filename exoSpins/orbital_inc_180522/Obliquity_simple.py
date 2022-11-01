import numpy as np
import matplotlib.pyplot as plt
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

sample = (vsini < v)
nsample = np.sum(sample)
print(f"Size new sample = {nsample} (vsini < v < vb)")
sini = vsini/v
ip = np.arcsin(sini[sample])
ip = np.concatenate((ip, np.pi-ip))
nsample *= 2

file = open("io.dat", "r")
io = np.loadtxt(file, skiprows=1)
file.close()

if len(io) > nsample:
    io = np.random.choice(io, nsample)
else:
    nsample = len(io)
    ip = np.random.choice(ip, nsample)
lambda_no_prior = abs(io-ip) # lambda =  projected obliquity |io-ip|

plt.hist(lambda_no_prior/dr, bins="auto", density=True)

plt.xlim(xmin=0)
plt.ylim(ymin=0)

plt.xlabel(r"Projected Obliquity $|i_{\rm p} - i_{\rm o}|$ (deg)")
plt.ylabel("PDF")

plt.tight_layout()
plt.show()
