from Neuron_Group import *

class Constant_Poisson(Neuron_Group):
    
    def __init__(self,qty,rate=50):
        
        super(Constant_Poisson,self).__init__(qty)

        self.rate=rate
        
    def update(self,t):
        
        self.swap_spikes()

        self.spike[:]=numpy.random.rand(self.quantity)<(self.rate/1000.0/self.dt)
    
    # replacing the update method
#    update=wrapAsMethod(update_methods.Constant_Poisson_update)
        
    
