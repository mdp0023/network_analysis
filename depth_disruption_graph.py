# This python code creates the depth disruption graph from the pregnolato paper

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,150,15000)
y = (0.003864*x**2-1.1592*x+86.94)/86.94

fig,ax = plt.subplots(figsize=(6,2.5))

plt.plot(x,y, color='black')

ax.set_ylabel("Percent of Maximum\nSpeed Limit",  fontsize=12)
ax.set_xlabel("Depth of Water on Road (mm)", fontsize=12)

plt.grid()
fig.tight_layout()
plt.show()