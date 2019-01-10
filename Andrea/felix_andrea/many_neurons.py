
import numpy as np
import nengo
import matplotlib.pyplot as plt


model = nengo.Network(label='many neurons')

with model: 
    A = nengo.Ensemble(100, dimensions=1)

    sin = nengo.Node(lambda t: np.sin(8 * t))

    con2 = nengo.Connection(sin, A, synapse = 0.01)

    sin_probe = nengo.Probe(sin)
    A_probe = nengo.Probe(A, synapse=0.01)
    A_spikes = nengo.Probe(A.neurons)

with nengo.Simulator(model) as sim:
    sim.run(1)

from nengo.utils.matplotlib import rasterplot

plt.plot(sim.trange(), sim.data[sin_probe], 'r', label ='A in')
plt.plot(sim.trange(), sim.data[A_probe], label='A out')
plt.xlim(0,1)
plt.legend()
plt.show()
    