import nengo
import random
import math

from nengo.utils.functions import piecewise
import matplotlib.pyplot as plt

model = nengo.Network(seed=200)
with model:
    nengo.neurons.LIF(tau_rc=0.01, tau_ref=0.001)
    
    # encoders for DA Ser
    DAencoders=[]              # create an empty list to store the encoders
    for i in range(1000):
        x=abs(random.uniform(0,1))  
        #y=random.uniform(-1,1)
        #z=random.uniform(-1,1)
        #DAencoders.append([x,y,z])  # add the encoder
        DAencoders.append([x])  # add the encoder

    Serencoders=[]              # create an empty list to store the encoders
    for i in range(1000):
        x=abs(random.uniform(0,1))*(-1)
        #y=random.uniform(-1,1)
        #z=random.uniform(-1,1)
        #Serencoders.append([x,y,z])  # add the encoder
        Serencoders.append([x])  # add the encoder
    
    # parameters
    eta   = 1   #@
    mu    = 0.6 #@DLPFC-ACC Ser signal (additional signal for cost of negative value)
    sigma = 2.0 #@DA and Ser PE rate output
    gamma = 0.6 #@Amy connection from Ser
    beta  = 0.3 #@Amy connection from DA
    alpha = 0.2 #@DA and Ser PE rate recurrent (learning rate)
    
    # Alex multiplication
    # we're using a 2D vector in OFC (representing [amy;everything else])
            # Define a function that computes the multiplication of two inputs
    def product(x):
        return x[0] * x[1]
        
    def no_amy(x):
        return x[0]
    
    # function computed by VS
    def E_func(x):
        return x[0] + x[1]
        
    # functions computed by ACC
    def R_func(x):
        Bf = 2*(x[0]>=0)-1
        return Bf/(eta+x[1])
        
    def ACC_id_func(x):
        return x[2]
        
    # functions computed by Amy
    def A_func(x):
        return x[0]+beta*x[1]+gamma*x[2]+x[3]
        
    # functions computed by DLPFC
    def C_func(x):
        return mu*x[1]

    # functions computed by DA
    def E_pl_func(x):
        return(x[0] - x[1])
        
    def DA_func(x):
        return sigma*x[0]
		
	def DA_func_c(x):
		return sigma*x[0]-(1-sigma)*x[1]
        
    def P_func_DA(x):
        return x[0] + alpha*x[1]
        
    # functions computed by Ser
    # use the same as in DA but the encoder will be "-1" to encode neg PE
    def E_mi_func(x):
        return(x[0] - x[1])

    def HT_func(x):
        return sigma*x[0]
		
	def HT_func_c(x):
		return sigma*x[0]-(1-sigma)*x[1]
        
    def P_func_Ser(x):
        return x[0] + alpha*x[1]
    
    stim_1 = nengo.Node([0])
    #stim_2 = nengo.Node([0])
    stim_2 = nengo.Node(piecewise({0: 0, 3: 2.4, 7: -2.4, 10: 1.0, 13: 0.0}))
    
    amygdala = nengo.Ensemble(n_neurons=1000, dimensions=4,radius=5)
    OFC = nengo.Ensemble(n_neurons=1000, dimensions=2,radius=6,max_rates=nengo.dists.Uniform(10, 80))
    ACC = nengo.Ensemble(n_neurons=1000, dimensions=3,max_rates=nengo.dists.Uniform(10, 80),radius=3)
    VS = nengo.Ensemble(n_neurons=1000, dimensions=2,radius=4)
    DLPFC = nengo.Ensemble(n_neurons=1000, dimensions=2, max_rates=nengo.dists.Uniform(10, 80))
    
    DA_Epl   = nengo.Ensemble(n_neurons=600, dimensions=2,radius=4)
    Ser_Eneg = nengo.Ensemble(n_neurons=600, dimensions=2,radius=4)
    
    DA = nengo.Ensemble(n_neurons=1000, dimensions=1, radius = 1,
    intercepts=nengo.dists.Uniform(0, 1), encoders = DAencoders)
    Ser = nengo.Ensemble(n_neurons=1000, dimensions=1,radius = 1,
    intercepts=nengo.dists.Uniform(0, 1),  encoders =Serencoders)
    
    P_Ser = nengo.Ensemble(n_neurons=2000, dimensions=2, radius = 10)
    P_DA  = nengo.Ensemble(n_neurons=1000, dimensions=2, radius = 10)
    
    output = nengo.Ensemble(n_neurons=600, dimensions=1,radius=4)
    
    ## Connections
    # stimulus
    nengo.Connection(stim_1, amygdala[0])
    nengo.Connection(stim_2, OFC[0])
    
    # Amygdala
    nengo.Connection(amygdala, OFC[1],function=A_func)
    
    # OFC
    nengo.Connection(OFC, ACC[0],function=no_amy)
    nengo.Connection(OFC, DA_Epl[0],function=no_amy)
    nengo.Connection(OFC, Ser_Eneg[0],function=no_amy,synapse=0.005)
    nengo.Connection(OFC, output[0],function=no_amy)
    
    # DA_Epl
    nengo.Connection(DA_Epl, DA[0],function= E_pl_func)
    nengo.Connection(DA_Epl,P_DA[1],function=E_pl_func)
    nengo.Connection(DA_Epl[1],P_DA[0])
    #nengo.Connection(DA_Epl,Ser[1],function= E_mi_func)
    
    # P_DA
    nengo.Connection(P_DA, DA_Epl[1],function=P_func_DA)

    # DA
    nengo.Connection(DA, VS[0],function=DA_func)
    nengo.Connection(DA, amygdala[1],function=DA_func)
    
    # Ser_Eneg
    nengo.Connection(Ser_Eneg, Ser[0],function= E_mi_func)
    nengo.Connection(Ser_Eneg,P_Ser[1],function=E_mi_func)
    nengo.Connection(Ser_Eneg[1],P_Ser[0])
    #nengo.Connection(DA_Epl,Ser[1],function= E_mi_func)
    
    # P_Ser
    nengo.Connection(P_Ser, Ser_Eneg[1],function=P_func_Ser)

    # Ser
    nengo.Connection(Ser, VS[1],function=HT_func)
    nengo.Connection(Ser, amygdala[2],function=HT_func)
    nengo.Connection(Ser, DLPFC[1],function=HT_func)
    
    nengo.Connection(VS, ACC[1], function = E_func)
    
    nengo.Connection(ACC, DLPFC[0],function=R_func)
    nengo.Connection(ACC, amygdala[3],function=ACC_id_func)
   
    nengo.Connection(DLPFC, ACC[2],function=C_func)
    # Add Probes
    # Anything that is probed will collect the data it produces over time, allowing us to analyze and visualize it later.
    stim_pb    = nengo.Probe(stim_2) # The raw spikes from the neuron
    spikes     = nengo.Probe(VS.neurons) # The raw spikes from the neuron
    voltage    = nengo.Probe(VS.neurons, 'voltage')  # Subthreshold soma voltage of the neuron
    filtered   = nengo.Probe(VS, synapse=0.01)       # Spikes filtered by a 10ms post-synaptic filter

# Run the Model
sim = nengo.Simulator(model) # Create the simulator
sim.run(15) # Run it for 1 seconds

# 6: Plot the Results
# Plot the decoded output of the ensemble
plt.plot(sim.trange(), sim.data[filtered])
plt.plot(sim.trange(), sim.data[stim_pb])
plt.xlim(0, 15)
plt.xlabel('time (s)')
plt.ylabel('electric potential')
plt.title('Postsynaptic Potentials at VS from VTA and DRN')

# # Plot the spiking output of the ensemble
# from nengo.utils.matplotlib import rasterplot
# plt.figure(figsize=(10, 8))
# plt.subplot(221)
# rasterplot(sim.trange(), sim.data[spikes])
# plt.ylabel("VS")
# plt.xlim(0, 15)

# # Plot the soma voltages of the neurons
# plt.subplot(222)
# plt.plot(sim.trange(), sim.data[voltage][:,0], 'r')
# plt.xlim(0, 15);
plt.show()