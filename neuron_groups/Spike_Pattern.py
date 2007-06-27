from Neuron_Group import *

class Spike_Pattern(Neuron_Group):

    
    def __init__(self,qty,pattern=[],rate=1):
        
        super(Spike_Pattern,self).__init__(qty)
    
        self.pattern=pattern
        if len(self.pattern)!=self.quantity:
            raise ValueError
        
        self.rate=rate
       
        self._reset_()
        
    def _reset_(self):
        
        # reset for time t=0
        
        super(Spike_Pattern,self)._reset_()
        self.time_to_next_spike=numpy.zeros(self.quantity,numpy.float)
        self.pattern_count=numpy.zeros(self.quantity)
        self.pattern_length=numpy.array([len(x) for x in self.pattern])

        self.time_to_next_pattern=0
        self.time_for_last_pattern=0
        
                
    def update(self,t):

        self.swap_spikes()
        
        if t>=self.time_to_next_pattern:
            self.time_for_last_pattern=t
            self.time_to_next_pattern=self.time_to_next_pattern+1000*self.dt/self.rate
            
            self.pattern_count[:]=0
            
            for i,pattern in enumerate(self.pattern):
                self.time_to_next_spike[i]=pattern[0]+t
                
        ttn=self.time_to_next_spike
        
        
        for i in range(self.quantity):
            if t>=ttn[i]:
                self.spike[i]=1
                
                self.pattern_count[i]=self.pattern_count[i]+1
                if self.pattern_count[i]==self.pattern_length[i]: # the end
                    ttn[i]=numpy.inf
                else:
                    ttn[i]=self.pattern[i][self.pattern_count[i]]+self.time_for_last_pattern
            else:
                spike[i]=0

    # replacing the update method
#    update=wrapAsMethod(update_methods.Spike_Pattern_update)

