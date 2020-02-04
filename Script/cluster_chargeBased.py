from cluster import HierarchicalClustering
from sklearn.cluster import DBSCAN
import numpy as np
from plotly import tools
import matplotlib.pyplot as plt

def _cluster_chargeBased_(position,acc):
    
   clustering = DBSCAN(eps=acc, min_samples=3).fit(position.reshape(-1, 1))
   s = np.linspace(-1,1)
   e = clustering.labels_
   center=clustering.components_
   #plt.plot(s, e)
   return e


def _cluster_dataBased(position,acc):
    return acc
