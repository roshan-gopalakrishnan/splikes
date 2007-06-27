from Connection_Group import *

class Gerstner06(Connection_Group):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),name='hippo full',
                      sign=1,initial_weight_range=[0.0,1.0],
                      valid_weight_range=[0.0,1e500],
                      initial_weights=None):
        
        super(Gerstner06,self).__init__(incell,outcell,
                                        sign,initial_weight_range,
                                        valid_weight_range,
                                        initial_weights)

                                          
        self.nearest_neighbor=False
        
        
        if name=='hippo full':
            # hippo full 
            self.A2_plus=6.1e-3
            self.A3_plus=6.7e-3
            self.A2_minus=1.6e-3
            self.A3_minus=1.4e-4
            self.tau_x=946.
            self.tau_y=27.
        elif name=='hippo min':
            self.A2_plus=5.3e-3
            self.A3_plus=8e-3
            self.A2_minus=3.5e-3
            self.A3_minus=0.0
            self.tau_x=946.  # not needed?
            self.tau_y=40.
        elif name=='hippo full nearest':
            self.nearest_neighbor=True
            self.A2_plus=4.6e-3
            self.A3_plus=9.1e-3
            self.A2_minus=3e-3
            self.A3_minus=7.5e-9
            self.tau_x=575.
            self.tau_y=47.
        elif name=='hippo min nearest':
            self.nearest_neighbor=True
            self.A2_plus=4.6e-3
            self.A3_plus=9.1e-3
            self.A2_minus=3e-3
            self.A3_minus=0.0
            self.tau_x=575.  # not needed?
            self.tau_y=48.
        elif name=='vc full':
            # vc full 
            self.A2_plus=5e-10
            self.A3_plus=6.2e-3
            self.A2_minus=7e-3
            self.A3_minus=2.3e-4
            self.tau_x=101.
            self.tau_y=125.
        elif name=='vc min':
            self.A2_plus=0.
            self.A3_plus=6.5e-3
            self.A2_minus=7.1e-3
            self.A3_minus=0.0
            self.tau_x=101.  # not needed?
            self.tau_y=114.
        elif name=='vc full nearest':
            self.nearest_neighbor=True
            self.A2_plus=8.8e-11
            self.A3_plus=5.3e-2
            self.A2_minus=6.6e-3
            self.A3_minus=3.1e-3
            self.tau_x=714.
            self.tau_y=40.
        elif name=='vc min nearest':
            self.nearest_neighbor=True
            self.A2_plus=0.
            self.A3_plus=5e-2
            self.A2_minus=8e-3
            self.A3_minus=0
            self.tau_x=714.  # not needed?
            self.tau_y=40.
        else:
            print "Warning...using default for Gerstner06"
            # hippo full 
            self.A2_plus=6.1e-3
            self.A3_plus=6.7e-3
            self.A2_minus=1.6e-3
            self.A3_minus=1.4e-4
            self.tau_x=946
            self.tau_y=27
            

        self.tau_plus=16.8
        self.tau_minus=33.7
        self.incell=incell
        self.outcell=outcell
        
        self._reset_()
        
    def _reset_(self):
        super(Gerstner06,self)._reset_()

        self.R1=numpy.zeros(self.incell.quantity,numpy.float)
        self.R2=numpy.zeros(self.incell.quantity,numpy.float)
        self.O1=numpy.zeros(self.outcell.quantity,numpy.float)
        self.O2=numpy.zeros(self.outcell.quantity,numpy.float)
        
    def update(self, t):
    
        
        out_active=numpy.nonzero(self.outcell.spike)
        in_active=numpy.nonzero(self.incell.spike)

        # update r1 and o1
        self.O1[:]-=self.O1*self.dt/self.tau_minus
        
        if self.nearest_neighbor:
            self.O1[out_active]=1.0
        else:
            self.O1[out_active]+=1.0

        
        self.R1[:]-=self.R1*self.dt/self.tau_plus
        
        if self.nearest_neighbor:
            self.R1[in_active]=1.0
        else:
            self.R1[in_active]+=1.0

            
        # update the weights

        self.weights[:]+=( 
           (-self.O1*self.incell.spike[:,newaxis]*(
           self.A2_minus+self.A3_minus*self.R2[:,newaxis]))+ (
           
             self.R1[:,newaxis]*self.outcell.spike * (
           self.A2_plus+self.A3_plus*self.O2[:,newaxis])))*self.sign

                
                
        # update r2 and o2 afterwards
        self.O2[:]-=self.O2*self.dt/self.tau_y
        
        if self.nearest_neighbor:
            self.O2[out_active]=1.0
        else:
            self.O2[out_active]+=1.0

        
        self.R2[:]-=self.R2*self.dt/self.tau_x
        
        if self.nearest_neighbor:
            self.R2[in_active]=1.0
        else:
            self.R2[in_active]+=1.0
                

    # replacing the update method
#    update=wrapAsMethod(update_methods.Gerstner06_update)
