from Neuron_Group import *

class Stochastic_Rate(Neuron_Group):
    
    def __init__(self,qty,activation_magnitude=1,tau_activation=100):
        
        super(Stochastic_Rate,self).__init__(qty)
        
        self.activation_magnitude=activation_magnitude
        self.tau_activation=tau_activation
    
    
        self.y=numpy.ones(qty,numpy.float)*0.0

        
    def _reset_(self):
        super(Stochastic_Rate,self)._reset_()
        self.y=numpy.ones(self.quantity,numpy.float)*0.0
        
    def update(self,t):
        
        
        self.swap_spikes()

        # connections to this neuron
        g_sum=0
        for cg in self.connections_to:
            active=numpy.flatnonzero(cg.incell.old_spike)
            cg.g[active,:]+=self.activation_magnitude*cg.weights[active,:]*cg.sign
            g_sum+=cg.g.sum()
            cg.g[:]-=(1.0/self.tau_activation)*cg.g*self.dt

    
        self.y+=(1.0/self.tau_activation)*(g_sum-self.y)*self.dt

        
        self.spike[:]=numpy.random.rand(self.quantity)<(self.y/1000.0/self.dt)
        
    # replacing the update method
#    update=wrapAsMethod(update_methods.Stochastic_Rate_update)
        
