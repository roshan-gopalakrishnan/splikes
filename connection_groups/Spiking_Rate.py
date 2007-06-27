from Connection_Group import *

class Spiking_Rate(Connection_Group):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),
                      sign=1,initial_weight_range=[0.0,1.0],
                      valid_weight_range=[0.0,1e500],
                      initial_weights=None):
        
        super(Spiking_Rate,self).__init__(incell,outcell,
                                          sign,initial_weight_range,
                                          valid_weight_range,
                                          initial_weights)
        sec=1000
        self.tau_thresh=20*sec
        self.activation_magnitude=1
        self.tau_activation=100
        self.eta=1e-5
        self.thresh_o=25
        self.lambda_decay=0
        self.yo=0 # spontaneous level
        self.xo=0 # spontaneous level
        self.learning_rule=1  # 1 - bcm, 2 - lawcooper, 3 - hebb
            
        self.use_beta=False
        self.tau_beta=0
        
        self.x=numpy.zeros(incell.quantity,numpy.float)
        self.xtmp=numpy.zeros(incell.quantity,numpy.float)
        self.y=numpy.zeros(outcell.quantity,numpy.float)
        self.ytmp=numpy.zeros(outcell.quantity,numpy.float)
        
        self.beta=numpy.zeros(outcell.quantity,numpy.float)
        
        self.th=numpy.ones(outcell.quantity,numpy.float)*0.1
        
    def _reset_(self):
        super(Spiking_Rate,self)._reset_()
        
        self.x=numpy.zeros(self.incell.quantity,numpy.float)
        self.xtmp=numpy.zeros(self.incell.quantity,numpy.float)
        self.y=numpy.zeros(self.outcell.quantity,numpy.float)
        self.ytmp=numpy.zeros(self.outcell.quantity,numpy.float)
        
        self.th=numpy.ones(self.outcell.quantity,numpy.float)*0.1
        self.use_beta=self.tau_beta>0
        
        
    def update(self,t):

        out_active=numpy.nonzero(self.outcell.spike)
        in_active=numpy.nonzero(self.incell.spike)

        
        # convenience.  impact speed?
        th=self.th
        x=self.x
        xtmp=self.xtmp
        y=self.y
        ytmp=self.ytmp
        
        xtmp[in_active]+=self.activation_magnitude*self.sign
        
        x[:]+=(1.0/self.tau_activation)*(xtmp-x)*self.dt
        xtmp[:]-=(1.0/self.tau_activation)*xtmp*self.dt
        
        ytmp[out_active]+=self.activation_magnitude

        y[:]+=(1.0/self.tau_activation)*(ytmp-y)*self.dt
        ytmp[:]-=(1.0/self.tau_activation)*ytmp*self.dt


        th[:]+=(1.0/self.tau_thresh)*((y-self.yo)**2/self.thresh_o-th)
        
        # matches plasticity, not soctagon
        if self.learning_rule==1: # BCM
            self.weights+=self.eta*((x[:,newaxis]-self.xo)*(y-self.yo)*((y-self.yo)-th)-self.lambda_decay*self.weights)
        elif self.learning_rule==2: # law and cooper
            self.weights+=self.eta*((x[:,newaxis]-self.xo)*(y-self.yo)*((y-self.yo)-th)/th-self.lambda_decay*self.weights)
        elif self.learning_rule==3: # hebb w/oja
            self.weights+=self.eta*((x[:,newaxis]-self.xo)*(y-self.yo)-self.lambda_decay*weights-y**2*self.weights)
                
        

    # replacing the update method
#    update=wrapAsMethod(update_methods.Spiking_Rate_update)

