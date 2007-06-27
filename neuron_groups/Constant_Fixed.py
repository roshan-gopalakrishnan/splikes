from Neuron_Group import *

class Constant_Fixed(Neuron_Group):
    
    def __init__(self,qty,rate=50,min_latency=0.0,max_latency=1000/50.0):
        
        super(Constant_Fixed,self).__init__(qty)
        
        
        self.rate=rate
        self.min_latency=min_latency
        self.max_latency=max_latency
        
        self.time_to_next_spike=numpy.zeros(self.quantity,numpy.float)
     
    def _reset_(self):
        super(Constant_Fixed,self)._reset_()
        self.time_to_next_spike=numpy.random.randint(self.min_latency,self.max_latency+1,self.quantity).astype(numpy.float)
        
    def update(self,t):
        
        self.swap_spikes()
        
        ttn=self.time_to_next_spike
        
        
        self.spike[:]=(t>=ttn)
        ttn[self.spike]+=1000.0*self.dt/self.rate


#    update=wrapAsMethod(update_methods.Constant_Fixed_update)

