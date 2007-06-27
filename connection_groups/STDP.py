from Connection_Group import *

class STDP(Connection_Group):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),
                      sign=1,initial_weight_range=[0.0,1.0],
                      valid_weight_range=[0.0,1e500],
                      initial_weights=None):
        
        super(STDP,self).__init__(incell,outcell,
                                          sign,initial_weight_range,
                                          valid_weight_range,
                                          initial_weights)
                                          
        
        self.a_plus=0.005
        self.a_minus=0.005*1.05
        self.tau_plus=20
        self.tau_minus=20
        self.g_max=1

        self.nearest_neighbor=False
        
        self._reset_()
        
        
    def _reset_(self):
        super(STDP,self)._reset_()
        self.P=numpy.zeros(self.incell.quantity,numpy.float)
        self.M=numpy.zeros(self.outcell.quantity,numpy.float)
        
        
    def update(self,t):

        out_active=numpy.nonzero(self.outcell.spike)
        in_active=numpy.nonzero(self.incell.spike)

        
        
        
        self.M[:]-=self.M*self.dt/self.tau_minus;
        if self.nearest_neighbor:
            self.M[out_active]=-self.a_minus*self.sign
        else:
            self.M[out_active]-=self.a_minus*self.sign

        
        self.P[:]-=self.P*self.dt/self.tau_plus;
        if self.nearest_neighbor:
            self.P[in_active]=self.a_plus*self.sign
        else:
            self.P[in_active]+=self.a_plus*self.sign
            

        self.weights[:]+=(self.M*self.incell.spike[:,newaxis]*self.g_max+
                    self.P[:,newaxis]*self.outcell.spike*self.g_max)*self.sign
        
    # replacing the update method
#    update=wrapAsMethod(update_methods.STDP_update)

    
