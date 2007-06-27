from Neuron_Group import *

class Integrate_and_Fire(Neuron_Group):
    
    def __init__(self,qty):
        
        super(Integrate_and_Fire,self).__init__(qty)
        
        self.tau_m=20
        self.tau_ex=5
        self.tau_in=5
        self.V_rest=-65
        self.reset_adaptation=0
        self.adaptation_step=2
        self.tau_adaptation=100
        self.V_thresh=-55
        self.t_refract=0
        self.V_rev_exc=0
        self.V_rev_inh=-65
        self.g_exc_max=3e-2
        self.g_inh_max=1e-1
    

        self._reset_()
        

    def _reset_(self):
        super(Integrate_and_Fire,self)._reset_()
        
        self.V=numpy.ones(self.quantity,numpy.float)*self.V_rest
        self.V_reset=numpy.ones(self.quantity,numpy.float)*self.V_rest
        self.time_to_stop_refract=numpy.zeros(self.quantity,numpy.float)
        
    def update(self,t):
        
        
        self.swap_spikes()

        nr=(self.time_to_stop_refract<=t)  # not refract

        # connections to this neuron
        gg_ex=0
        gg_in=0
        for cg in self.connections_to:
            active=numpy.flatnonzero(cg.incell.old_spike)
            if cg.sign>0:
                gmax=self.g_exc_max
            else:
                gmax=self.g_inh_max
                
            cg.g[active,:]+=gmax*cg.weights[active,:]
            
            if cg.sign>0:
                gg_ex+=cg.g.sum()
            else:
                gg_in+=cg.g.sum()
        
        # only update voltages of non-refractory neurons
        self.V[nr]=(self.V[nr]-(self.V[nr]-self.V_reset)*self.dt/self.tau_m+
                      self.dt*gg_ex*(self.V_rev_exc-self.V[nr])/self.tau_m +
                      self.dt*gg_in*(self.V_rev_inh-self.V[nr])/self.tau_m)
          
        
        self.spike[:]=(self.V>self.V_thresh)
        if self.reset_adaptation:
            self.V_reset[self.spike]-=self.adaptation_step
            
        self.V[self.spike]=self.V_reset[self.spike]
        self.time_to_stop_refract[self.spike]=t+self.t_refract
        

        # decay conductances
        for cg in self.connections_to:
            if cg.sign>0:
                cg.g[:]-=cg.g*self.dt/self.tau_ex
            else:
                cg.g[:]-=cg.g*self.dt/self.tau_in
      
        if self.reset_adaptation:
            self.V_reset[self.spike]-=(self.V_reset[self.spike]-self.V_rest)*self.dt/self.tau_adaptation

    # replacing the update method
#    update=wrapAsMethod(update_methods.Integrate_and_Fire_update)
    
