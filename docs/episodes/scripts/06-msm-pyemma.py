import numpy as np
import pyemma

x1 = np.random.normal(-1.0, 0.2, size=(1000, 1))
x2 = np.random.normal(1.0, 0.2, size=(1000, 1))
X = [np.vstack([x1, x2])]

tica = pyemma.coordinates.tica(X, lag=10)
Y = tica.get_output()

cl = pyemma.coordinates.cluster_kmeans(Y, k=20, max_iter=50)
msm = pyemma.msm.estimate_markov_model(cl.dtrajs, lag=10)
print('Timescales:', msm.timescales_[:5])
