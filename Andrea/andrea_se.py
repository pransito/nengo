import nengo
import random
import math

model = nengo.Network()
with model:
    nengo.neurons.LIF(tau_rc=0.01, tau_ref=0.001)
    
    # encoders for DA Ser
    DAencoders=[]              # create an empty list to store the encoders
    for i in range(1200):
        x=random.uniform(0,1)  
        y=random.uniform(-1,1)
        z=random.uniform(-1,1)
        DAencoders.append([x,y,z])  # add the encoder

    Serencoders=[]              # create an empty list to store the encoders
    for i in range(1200):
        x=random.uniform(0,1)  
        y=random.uniform(-1,1)
        z=random.uniform(-1,1)
        Serencoders.append([x,y,z])  # add the encoder
    
    # parameters
    eta   = 1   #@
    mu    = 0.6 #@DLPFC-ACC Ser signal (additional signal for cost of negative value)
    sigma = 0.8 #@DA and Ser PE rate output
    gamma = 0.6 #@Amy connection from Ser
    beta  = 0.3 #@Amy connection from DA
    alpha = 0.3 #@DA and Ser PE rate recurrent (learning rate)
    
    # Alex multiplication
    # we're using a 2D vector in OFC (representing [amy;everything else])
            # Define a function that computes the multiplication of two inputs
    def product(x):
        return x[0] * x[1]
    
    # function computed by VS
    def E_func(x):
        return x[0] - x[1]
        
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
        E_pl = x[0]-x[2]
        if E_pl > 0:
            return E_pl
        else:
            return 0
    
    def DA_func(x):
        E_pl = E_pl_func(x)
        return sigma*E_pl
        
    def P_func_DA(x):
        E_pl = x[0]-x[2]
        if E_pl > 0:
            E_pl = E_pl
        else:
            E_pl = 0
        return x[2] + alpha*E_pl
        
    # functions computed by Ser
    def E_mi_func(x):
        E_mi = x[2]-x[0]
        if E_mi > 0:
            return E_mi
        else:
            return 0
    
    def HT_func(x):
        E_mi = E_mi_func(x)
        return sigma*E_mi
        
    def P_func_Ser(x):
        E_mi = x[2]-x[0]
        if E_mi > 0:
            E_mi = E_mi
        else:
            E_mi = 0
        return x[2] + alpha*E_mi
    
    stim_1 = nengo.Node([0])
    stim_2 = nengo.Node([0])
    
    
    amygdala = nengo.Ensemble(n_neurons=1000, dimensions=4,radius=3)
    OFC = nengo.Ensemble(n_neurons=1000, dimensions=2,radius=4,max_rates=nengo.dists.Uniform(10, 80))
    ACC = nengo.Ensemble(n_neurons=1000, dimensions=3,max_rates=nengo.dists.Uniform(10, 80))
    VS = nengo.Ensemble(n_neurons=1000, dimensions=2)
    DLPFC = nengo.Ensemble(n_neurons=1000, dimensions=2, max_rates=nengo.dists.Uniform(10, 80))
    
    DA = nengo.Ensemble(n_neurons=1200, dimensions=3, radius = 5)
    Ser = nengo.Ensemble(n_neurons=1200, dimensions=3,radius = 5)
    
    #DA = nengo.Ensemble(n_neurons=1200, dimensions=3, radius = 3,
    #intercepts=nengo.dists.Uniform(0, 3),  encoders =DAencoders)
    #Ser = nengo.Ensemble(n_neurons=1200, dimensions=3,radius = 3,
    #intercepts=nengo.dists.Uniform(0, 3),  encoders =Serencoders)
    
    output = nengo.Ensemble(n_neurons=1200, dimensions=1,radius=4)
    
    nengo.Connection(stim_1, amygdala[0])
    nengo.Connection(stim_2, OFC[0])
    
    
    #nengo.Connection(amygdala, OFC[1],function=A_func)
    
    nengo.Connection(OFC, ACC[0],function=product)
    nengo.Connection(OFC[0], DA[0])
    nengo.Connection(OFC[0], Ser[0])
    nengo.Connection(OFC, output[0],function=product)
    
    nengo.Connection(DA, DA[2],function= P_func_DA)
    nengo.Connection(DA, VS[0],function=DA_func)
    nengo.Connection(DA, amygdala[1],function=DA_func)
    nengo.Connection(DA, Ser[1],function=E_pl_func)
    
    nengo.Connection(Ser, Ser[2],function=P_func_Ser,synapse=0.005)
    nengo.Connection(Ser, VS[1],function=HT_func,synapse=0.005)
    nengo.Connection(Ser, DLPFC[1],function=HT_func,synapse=0.005)
    nengo.Connection(Ser, DA[1],function=E_mi_func,synapse=0.005)
    
    nengo.Connection(Ser, Ser[2],function=P_func_Ser,synapse=0.005)
    nengo.Connection(Ser, amygdala[2],function=HT_func,synapse=0.005)
    
    nengo.Connection(VS, ACC[1], function = E_func)
    
    nengo.Connection(ACC, DLPFC[0],function=R_func)
    nengo.Connection(ACC, amygdala[3],function=ACC_id_func)
   
    
    nengo.Connection(DLPFC, ACC[2],function=C_func)
    
    
    
    
    
