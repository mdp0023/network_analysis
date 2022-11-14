import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
import numpy as np
from rasterio.enums import Resampling
from matplotlib.colors import ListedColormap

path = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Inundation/AN_Memorial_Day_Compound.tif'
spath = '/home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing/AOI_Inundation.tif'
fig, ax = plt.subplots(facecolor='grey', figsize=(12, 12))
ax.axis('off')


# Choose colormap
cmap = plt.cm.Blues
# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))
# Set alpha
my_cmap[:, -1] = np.linspace(0.5, 1, cmap.N)
my_cmap[0][3] = 0
# Create new colormap
my_cmap = ListedColormap(my_cmap)

# The grid of raster values can be accessed as a numpy array and plotted:
with rasterio.open(path) as src:
    oview = 1.5
   # NOTE this is using a 'decimated read' (http://rasterio.readthedocs.io/en/latest/topics/resampling.html)
    thumbnail = src.read(1, out_shape=(1, int(src.height // oview), int(src.width // oview)), resampling=Resampling.nearest)
    
    #thumbnail = thumbnail.astype('f4')
    thumbnail[thumbnail <= 0.01] = np.nan
    thumbnail[thumbnail >5] = 5

    


# plt.imshow(thumbnail, cmap='Blues')
# plt.colorbar()

show(thumbnail, cmap='Blues')


# for ji, w in raster.block_windows(1):
#     w = raster.read(1, window=w)
#     plt.imshow(w)

plt.show()
