import nengo

model = nengo.Network()
with model:
    stim_1 = nengo.Node([0])
    stim_2 = nengo.Node([0])
    amygdala = nengo.Ensemble(n_neurons=1000, dimensions=1)
    nengo.Connection(stim_1, amygdala)
    
    OFC = nengo.Ensemble(n_neurons=1000, dimensions=2)
    nengo.Connection(amygdala, OFC[0])
    nengo.Connection(stim_2, OFC[1])
    
    prod = nengo.Ensemble(1000, dimensions=1, radius=20)
    
        # Define a function that computes the multiplication of two inputs
    def product(x):
        return x[0] * x[1]
    
    # Connect the combined ensemble to the output ensemble D
    nengo.Connection(OFC, prod, function=product)
