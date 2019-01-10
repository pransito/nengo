

import nengo
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

#%load_ext nengo.ipynb


from nengo.dists import Uniform 

model = nengo.Network(label = 'A Single Neuron')
with model:
    neuron = nengo.Ensemble(1, dimensions=1, 
                            intercepts=Uniform(-0.5, -0.5), 
                            max_rates=Uniform(100,100), 
                            encoders=[[1]])
with model:
    cos = nengo.Node(lambda t: np.cos(8 * t))

with model:
    con = nengo.Connection(cos, neuron)

with model:
    cos_probe = nengo.Probe(cos)
    con_probe = nengo.Probe(con)
    spikes = nengo.Probe(neuron.neurons)
    voltage = nengo.Probe(neuron.neurons, 'voltage')
    pre_filtered = nengo.Probe(neuron)
    filtered = nengo.Probe(neuron, synapse=0.01)

with nengo.Simulator(model) as sim:
    sim.run(1)


from nengo.utils.matplotlib import rasterplot 

plt.subplot(4,1,1)

plt.plot(sim.trange(), sim.data[pre_filtered])
plt.plot(sim.trange(), sim.data[cos_probe])
plt.xlim(0, 1)
plt.ylim(-2,10)
plt.ylabel('w/o filter')

plt.subplot(4,1,2)
plt.plot(sim.trange(), sim.data[filtered])
plt.plot(sim.trange(), sim.data[cos_probe])
plt.xlim(0, 1)
plt.ylim(-1,1)
plt.ylabel('filtered')

plt.subplot(4,1,3)
rasterplot(sim.trange(), sim.data[spikes])
#plt.plot(sim.trange(), sim.data[cos_probe])
#plt.ylim(0, 3000)
plt.xlim(0, 1)
plt.ylabel('spikes')

plt.subplot(4,1,4)
#plt.plot(sim.trange(), sim.data[con_probe])
plt.plot(sim.trange(), sim.data[voltage][:,0])
plt.xlim(0,1)
plt.ylim(-2,2)
plt.ylabel('voltage')


plt.show()    