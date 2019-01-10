import numpy as np
import matplotlib.pyplot as plt

import nengo

from nengo.dists import Uniform

model = nengo.Network(label='Two Neurons')
with model:
    neurons = nengo.Ensemble(2, dimensions=1, intercepts=Uniform(-0.9, -0.9),
                             max_rates=Uniform(100,100), encoders=[[1], [-1]])
                             
                             
    sin = nengo.Node(lambda t: np.sin(8 * t))
    
    nengo.Connection(sin, neurons, synapse=0.01)
    
    sin_probe = nengo.Probe(sin)
    spikes = nengo.Probe(neurons.neurons)
    voltage = nengo.Probe(neurons.neurons, 'voltage')
    filtered = nengo.Probe(neurons)#, synapse=0.01)
    
with nengo.Simulator(model) as sim:
    sim.run(1)
    
    
plt.plot(sim.trange(), sim.data[filtered][:,0])
plt.plot(sim.trange(), sim.data[sin_probe])
plt.xlim(0, 1) 


from nengo.utils.matplotlib import rasterplot 
plt.figure(figsize=(10, 8))
plt.subplot(221)
rasterplot(sim.trange(), sim.data[spikes], colors=[(1, 0, 0), (0, 0, 0)])
plt.xlim(0, 1)
plt.yticks((0, 1), ("On", "Off"))

plt.subplot(222)
plt.plot(sim.trange(), sim.data[voltage][:, 0] + 1, 'r')
plt.plot(sim.trange(), sim.data[voltage][:, 1], 'b')
plt.yticks(())
plt.axis([0, 1, 0, 2])
plt.subplots_adjust(wspace=0.05);   




plt.show()