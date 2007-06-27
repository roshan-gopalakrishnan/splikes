from Neuron_Group import *

class Izh(Neuron_Group):
    
    def __init__(self,qty,name='class 1'):
        """
        Types of firing patterns:
        
             tonic spiking 0
             phasic spiking 1
             tonic bursting 2
             phasic bursting 3
             mixed mode 4
             spike frequency adaptation 5
             Class 1  6
             Class 2  7
             spike latency 8
             subthreshold oscillations 9
             resonator 10
             integrator 11
             rebound spike 12
             rebound burst 13
             threshold variability 14
             bistability 15
             DAP 16
             accomodation 17
             inhibition-induced spiking 18
             inhibition-induced bursting 19
        
        
        """
        
        super(Izh,self).__init__(qty)

        pars=([0.02  ,    0.2  ,   -65 ,     6 ,      14 ],    # tonic spiking 0
              [0.02  ,    0.25 ,   -65 ,     6 ,      0.5 ],   # phasic spiking 1
              [0.02  ,    0.2  ,   -50 ,     2 ,      15 ],    # tonic bursting 2
              [0.02  ,    0.25 ,   -55 ,    0.05,     0.6 ],   # phasic bursting 3
              [0.02  ,    0.2  ,   -55 ,    4  ,      10 ],    # mixed mode 4
              [0.01  ,    0.2  ,   -65 ,    8  ,      30 ],    # spike frequency adaptation 5
              [0.02  ,    -0.1 ,   -55 ,    6  ,      0  ],    # Class 1  6
              [0.2   ,    0.26 ,   -65 ,    0  ,      0  ],    # Class 2  7
              [0.02  ,    0.2  ,   -65 ,    6  ,      7  ],    # spike latency 8
              [0.05  ,    0.26 ,   -60 ,    0  ,      0  ],    # subthreshold oscillations 9
              [0.1   ,    0.26 ,   -60 ,    -1 ,      0  ],    # resonator 10
              [0.02  ,    -0.1 ,   -55 ,    6  ,      0  ],    # integrator 11
              [0.03  ,    0.25 ,   -60 ,    4  ,      0],      # rebound spike 12
              [0.03  ,    0.25 ,   -52 ,    0  ,      0],      # rebound burst 13
              [0.03  ,    0.25 ,   -60 ,    4  ,      0  ],    # threshold variability 14
              [1     ,    1.5  ,   -60 ,    0  ,    -65  ],    # bistability 15
              [  1   ,    0.2  ,   -60 ,    -21,      0  ],    # DAP 16
              [0.02  ,    1    ,   -55 ,    4  ,      0  ],    # accomodation 17
             [-0.02  ,    -1   ,   -60 ,    8  ,      80 ],    # inhibition-induced spiking 18
             [-0.026 ,    -1   ,   -45 ,    0  ,      80])       # inhibition-induced bursting 19
                        
                                     
        d={}
        d['tonic spiking']=0
        d['phasic spiking']=1
        d['tonic bursting']=2
        d['phasic bursting']=3
        d['mixed mode']=4
        d['spike frequency adaptation']=5
        d['class 1']=6
        d['class 2']=7
        d['spike latency']=8
        d['subthreshold oscillations']=9
        d['resonator']=10
        d['integrator']=11
        d['rebound spike']=12
        d['rebound burst']=13
        d['threshold variability']=14
        d['bistability']=15
        d['dap']=16  # depolarizing after potentials
        d['accomodation']=17
        d['inhibition-induced spiking']=18
        d['inhibition-induced bursting']=19
             

        try:
            idx=d[name]
        except KeyError:
            print "Warning...using default class 1"
            idx=d['class 1']
            
            

        self.a=pars[idx][0]
        self.b=pars[idx][1]
        self.c=pars[idx][2]
        self.d=pars[idx][3]
        self.I=pars[idx][4]

        self.tau_epsp1=50
        self.tau_epsp2=5
        self.epsp_scale=1.0
    
        self.V_reset=self.c
        self.V_peak=30

        self.initial_u=-20.0
        self.initial_V=-70.0
        
        self._reset_()
        
        import warnings

    def _reset_(self):
        super(Izh,self)._reset_()
        
        self.V=numpy.ones(self.quantity,numpy.float)*self.initial_V
        self.u=numpy.ones(self.quantity,numpy.float)*self.initial_u
        self.epsp1=numpy.zeros(self.quantity,numpy.float)
        self.epsp2=numpy.zeros(self.quantity,numpy.float)
        self.ipsp1=numpy.zeros(self.quantity,numpy.float)
        self.ipsp2=numpy.zeros(self.quantity,numpy.float)
        
        
    def update(self,t):
        
        self.V[self.spike]=self.V_reset
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
                gg_ex=0.0

        self.V=(self.V +(0.04*self.V*self.V+5.0*self.V+140.0-self.u)*self.dt+
              self.dt*gg_ex*(self.V_rev_exc-self.V)/self.tau_m +
              self.dt*gg_in*(self.V_rev_inh-self.V)/self.tau_m)

                          
            
        self.u=self.u+self.a*(self.b*self.V-self.u)*dt
            
        self.spike[:]=(self.V>self.V_peak)
        self.V[self.spike]=self.V_peak

        # decay conductances
        for cg in self.connections_to:
            if cg.sign>0:
                cg.g[:]-=cg.g*self.dt/self.tau_ex
            else:
                cg.g[:]-=cg.g*self.dt/self.tau_in
      
    # replacing the update method
#    update=wrapAsMethod(update_methods.Izh_update)
