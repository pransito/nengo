import nengo
from nengo.dists import Uniform
import numpy as np
from nengo.utils.functions import piecewise

_beta=2 # related to DA, reaction to appetitive stimuli
_gamma=3 #related to 5HT, reaction to aversive stimuli (gamma > beta)?
# draft paper says sigma = 0.75
_sigma=0.75 # < 1, determines the extent to which rectified information is 
        # combined together
_alpha=0.73 # between 0 and 1
_mu=1


N = 1000
n = 400
rad = 1
syn=0.003

with nengo.Network() as model:
    
    #model.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
    # AMYG
    amyg = nengo.Ensemble(N, 1, radius = 4, max_rates=(np.random.rand(N)*100+100).astype(int))
    stim_amyg = nengo.Node([1])
    #stim_amyg = nengo.Node(lambda t: 0.3 if 0.0<t<0.25
    #                        else 0 if 0.25<t<0.43
    #                        else 1)
    nengo.Connection(stim_amyg, amyg)
    
    # OFC
    ofc = nengo.Ensemble(N, 2, max_rates=(np.random.rand(N)*100+100).astype(int), 
                        radius=3)
    stim_ofc = nengo.Node(lambda t: 1.4 if 0.0<t<0.2 
                            else 0.9 if 0.2<t<0.35 
                            else -1 if 0.35<t<0.5
                            else 0.5 if 0.5<t<0.65
                            else -0.6 if 0.65<t<0.8
                            else 1)
    #stim_ofc = nengo.Node(lambda t: 1.4 if 0.0<t<0.4 
    #                        else 0.9 if 0.4<t<0.7 
    #                        else -1 if 0.7<t<1.0
    #                        else 0.5 if 1.0<t<1.26
    #                        else -0.7 if 1.26<t<1.6
    #                        else 1)
    #stim_ofc = nengo.Node(lambda t: 1.5 if 0.0<t<0.3
    #                        else -1.5 if 0.3<t<0.6
    #                       else 1.5)
    #stim_ofc = nengo.Node(lambda t: 0 if 0.0<t<0.15
    #                        else 2 if 0.15<t<0.3
    #                        else 0 if 0.3<t<0.45
    #                        else 1)
    nengo.Connection(stim_ofc, ofc[0])
   
    # 5-HT   
    ht5_Eminus = nengo.Ensemble(n, 1, encoders=np.ones((n,1)), 
                            intercepts=Uniform(low=0.0, high=1.0), 
                            max_rates=(np.random.rand(n)*100+100).astype(int),
                            radius=rad)
    ht5_ht5 = nengo.Ensemble(n, 1, max_rates=(np.random.rand(n)*100+100).astype(int),
                            radius=rad)
    ht5_P = nengo.Ensemble(n, 1, max_rates=(np.random.rand(n)*100+100).astype(int),
                            radius=4)
    
    # DA
    da_Eplus = nengo.Ensemble(n, 1, encoders=np.ones((n,1)), 
                            intercepts=Uniform(low=0.0, high=1.0), 
                            max_rates=(np.random.rand(n)*100+100).astype(int),
                            radius=rad)
    da_da = nengo.Ensemble(n, 1, max_rates=(np.random.rand(n)*100+100).astype(int),
                            radius=rad)
    da_P = nengo.Ensemble(n, 1, max_rates=(np.random.rand(n)*100+100).astype(int),
                            radius=4)
    
    # VS
    vs = nengo.Ensemble(N, 1, max_rates=(np.random.rand(N)*100+100).astype(int),
                        radius=rad)
    
    # ACC - behaviour
    acc = nengo.Ensemble(N, 1, max_rates=(np.random.rand(N)*100+100).astype(int))
                        
    # DLPFC - which does nothing, really
    #dlpfc = nengo.Ensemble(N, 1, max_rates=(np.random.rand(N)*100+100).astype(int))
    
    
    # Connections
    
    nengo.Connection(amyg, ofc[1])
    def simple_product(x):
        return x[0]*x[1]
    def neg_product(x):
        return -x[0]*x[1]
    
    nengo.Connection(ofc, acc, function=simple_product, synapse=syn)
    
    nengo.Connection(ofc, ht5_Eminus, function=neg_product, synapse=syn)
    nengo.Connection(ht5_P, ht5_Eminus, synapse=syn)
    
    nengo.Connection(ofc, da_Eplus, function=simple_product, synapse=syn)
    nengo.Connection(da_P, da_Eplus, transform=-1, synapse=syn)
    
    nengo.Connection(ht5_Eminus, ht5_ht5, transform=_sigma, synapse=syn)
    nengo.Connection(da_Eplus, ht5_ht5, transform=-(1-_sigma), synapse=syn)
    
    nengo.Connection(da_Eplus, da_da, transform=_sigma, synapse=syn)
    nengo.Connection(ht5_Eminus, da_da, transform=-(1-_sigma), synapse=syn)
    
    nengo.Connection(da_P, da_P, synapse=0.01)
    nengo.Connection(vs, da_P, transform=_alpha, synapse=syn)
    
    nengo.Connection(ht5_P, ht5_P, synapse=0.01)
    nengo.Connection(vs, ht5_P, transform=_alpha, synapse=syn)
    
    nengo.Connection(da_da, vs, synapse=syn)
    nengo.Connection(ht5_ht5, vs, transform=-1, synapse=syn)
    
    nengo.Connection(da_da, amyg, transform=_beta, synapse=syn)
    nengo.Connection(ht5_ht5, amyg, transform=_gamma, synapse=syn)
    
    #nengo.Connection(ht5_ht5, dlpfc, transform=_mu)
    #nengo.Connection(dlpfc, acc[1])
    #nengo.Connection(acc[1], amyg)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    