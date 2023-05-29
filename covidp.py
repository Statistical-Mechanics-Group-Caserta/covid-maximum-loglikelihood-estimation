import numpy as np
import math

seed1 = -39814113
ZBQLPOI = lambda x: np.random.poisson(x)
nit = np.zeros(100001)
ntest = np.zeros(50001)
bra = np.zeros(1000001)
muaat = 0.0
time = 0.0
tt0 = 0.0
itime = np.zeros(1000000010)



taup = 0.15
tau = taup/100000
ap = 6.2 / taup
gnorm = 1.0 / (math.gamma(ap) * taup ** ap)
xt = 0.0
t_rand = np.zeros(10000001)
gp = 0.0
g = 0.0
stringa = ""
ifinal = 0
nevent = 0
jfinal = 0
ninf = 250
flag_test = False
i=0
# Loop 1
for j in range(1, int(taup / tau) * int(60 / taup)):
    xt = j * tau
    g = g + tau * (xt ** (ap - 1.0) * np.exp(-xt * 1.0 / taup))
    if g * gnorm >= gp + 1.0 / 6000:
        i = i + 1
        t_rand[i] = xt
        gp = g * gnorm
        print(xt, g * gnorm)

ifinal = i

# Read file
with open("rt.dat", "r") as f:
    stringa = f.readline()
    i = 0
    for line in f:
        data = line.split()
        j = int(data[0])
        bra[i] = float(data[1])
        xnimm = float(data[2])
        for k in range(nevent + 1, nevent + int(xnimm) + 1):
            itime[k] = i
        nit[i] = int(xnimm)
        i += 1
        print(j, bra[i])

jfinal = i - 1

# Generate initial infected people
for j in range(1, ninf + 1):
    nevent += 1
    time = 1
    itime[nevent] = time
    nit[int(time)] += 1

# Loop 2
for j in range(1, 100000000):
    tt0 = itime[j]
    muaat = bra[int(tt0)]
    nn = ZBQLPOI(muaat)
    if itime[j] < 0:
        break
    for i in range(1, nn + 1):
        rr = np.random.random()
        time = tt0 + t_rand[int(rr * ifinal)] + 0.5
        if time <= jfinal:
            nevent += 1
            itime[nevent] = time
            nit[int(time)] += 1

# Print results
with open("output.dat", "w") as f:
    for j in range(1, jfinal):
        f.write(str(j + 55) + " " + str(nit[j] * (1 + ntest[j] * flag_test)) + " " + str(nit[j]) + "\n")
