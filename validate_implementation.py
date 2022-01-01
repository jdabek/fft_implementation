import numpy as np
import matplotlib.pyplot as plt

fft = []
with open("fft.txt") as f:
    for line in f:
        l = line.strip()
        x = float(l.split("(")[1].split(",")[0])
        y = float(l.split(",")[1].split(")")[0])
        fft.append([x, y])

fft = np.array(fft)

plt.plot(fft)
plt.legend(['x', 'y'])
plt.savefig("fft.png")
