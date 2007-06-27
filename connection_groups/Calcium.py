from Connection_Group import *
from numpy import exp

class Calcium(Connection_Group):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),
                      sign=1,initial_weight_range=[0.0,1.0],
                      valid_weight_range=[0.0,1e500],
                      initial_weights=None):
        
        super(Calcium,self).__init__(incell,outcell,
                                          sign,initial_weight_range,
                                          valid_weight_range,
                                          initial_weights)

        self.tau_ep1=50
        self.tau_ep2=5
        self.backspike_delay=0
        self.tau_backspike_slow=30
        self.tau_backspike_fast=3
        self.tau_ca=20
        self.v_reversal=130
        self.i_nmda_f=0.75
        self.i_nmda_s=0.25
        self.tau_nmda_s=200
        self.tau_nmda_f=50
          # 1 = sigmoid omega, 2 = quadratic, 3 = harel's eta function
          # 4 = biocyb eta and omega
          # 
        self.learning_rule=1  
          #self.k_plus=3e-5
          #self.k_minus=3e-7
        self.k_plus=9e-6
        self.k_minus=9e-8
        self.g_t=-4.5e-3  # what is this?
        self.Vp=2
        self.Vo=-65
        self.lambda_decay=0
        self.omega_offset=0.0
        self.mg1=-0.062
        self.mg2=3.57
        self.i_nmda_mu=0.8
        self.eta_gamma0=5e-5 # 1/5000
        self.eta_up=0.5*12*self.eta_gamma0
        self.eta_down=0.2*12*self.eta_gamma0
        self.omega_max=0.5
        self.theta_o=4.8/12.0
        self.backspike_amplitude=60
        self.peak_backspike_slow=0.25
        self.peak_backspike_fast=0.75
          
          # yeung's
        self.alpha1=0.25
        self.alpha2=0.4
        self.beta1=60
        self.beta2=20
          
          #harel's
        self.alpha1=0.4
        self.alpha2=0.65
        self.beta1=30
        self.beta2=30
          
          
        self.rest_voltage=-65
        self.v_total_use_V=1
        self.v_total_use_BP=1
        self.eta_p1=1
        self.eta_p2=.6
        self.eta_p3=3
        self.eta_p4=1e-5

        
        self._reset_()
        
        
    def _reset_(self):
        super(Calcium,self)._reset_()
        
        self.v_backspike_slow=numpy.zeros(self.outcell.quantity,numpy.float)
        self.v_backspike_fast=numpy.zeros(self.outcell.quantity,numpy.float)
        self.v_total=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.B=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.h=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.g_nmda=-1.0/400.0*numpy.ones((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.I_nmda=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.I_nmda_fast=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.I_nmda_slow=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.Ca=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.eta=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.omega=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        
        
        try:
            self.user_reset(self)
        except AttributeError:
            pass

    def update(self,t):

        
##  /*  from yeung's code 11/03 
##  
##      if(is_meta)
##                meta_dyn(g, v, t, p, ic); 
##    nmda_dyn(i_nmda, exc_spikes, g, v, t, p, ic);
##    ca_dyn(ca, i_nmda, t, p, ic);
##    learn_ca(w, ca, t, p, ic); 
##  
##    
##    yeung does: input, then learning, then output.
##    I do: input and output, then learning
##    
##    she also does her metaplasticity before the rest of the learning. 
##    
##    bkspike_dyn(&v_bk, post_spike_train, t, p, ic);
##    v = v_rest + v_bk; 
##
##    and she just does rest + v_bk
##  */

        out_active=numpy.nonzero(self.outcell.spike)
        in_active=numpy.nonzero(self.incell.spike)

        
        # set active backspikes to peak value
        self.v_backspike_slow[out_active]=self.peak_backspike_slow
        self.v_backspike_fast[out_active]=self.peak_backspike_fast

        self.v_total[:]*=0.0  # erase
        try:
            self.v_total[:]+=self.outcell.V
        except NameError:
            self.v_total[:]+=self.rest_voltage
            
        self.v_total[:]+=self.backspike_amplitude*(self.v_backspike_fast+self.v_backspike_slow)
        

        # meta plasticity
        
        if self.k_minus !=0 or self.k_plus!=0:  # do the metaplasticity

            ##      /* from yeung
            ##       tau_meta = 1.0/(k_minus*pow(v+65,2)+ k_plus);
            ##       g[i] = g[i] + t_scale/tau_meta * (g_t*k_plus*tau_meta - g[i]);
            ##      */

            self.g_nmda[:]+=self.dt*(self.k_plus*self.g_t-
             (self.k_plus+self.k_minus*(self.v_total-self.Vo)**self.Vp)*self.g_nmda)
            

             
             
        # NMDA dynamics
        

        self.I_nmda_slow[in_active,:]+=self.i_nmda_mu*(self.i_nmda_s-self.I_nmda_slow[in_active,:])
        self.I_nmda_fast[in_active,:]+=self.i_nmda_mu*(self.i_nmda_f-self.I_nmda_fast[in_active,:])
        
        self.B[:] = 1.0/(1.0 + (exp(self.mg1 * self.v_total) / self.mg2))
	    
        # this h is negative 
        self.h[:]= self.B*(self.v_total-self.v_reversal)

        # this g_nmda is negative 
        self.I_nmda[:]= self.g_nmda*(self.I_nmda_fast + self.I_nmda_slow) * self.h

        self.Ca[:]+=(self.I_nmda - self.Ca/self.tau_ca)*self.dt
        
        # p1=0.1; p2=p1*1e-4; p3=3; p4=1; % new paper % p2 ~=p1/1e-4 */
        # p1=0.5; p2=1e-4; p3=3; p4=0.05; % old paper */

        #  tau=p1[i]/(p2[i]+pow(Ca[i],p3[i]))+p4[i];   tau in seconds 
	# eta=1.0/(tau*1000.0);   % eta in 1/ms */

	    
        #  eta=p1[i]*Ca[i]/1000; 
	    
        if self.learning_rule==1: # difference of sigmoids, omega 
            self.omega[:]=self.omega_offset+(sig(self.Ca-self.alpha2,self.beta2)-
                                      0.5*sig(self.Ca-self.alpha1,self.beta1))
            self.eta[:]=self.eta_gamma0*self.Ca
        elif self.learning_rule==2: # quadratic
            self.omega[:]=self.omega_offset+ self.Ca*(self.Ca-self.theta_o)
            self.omega[:]=numpy.where(self.omega>self.omega_max,self.omega_max,self.omega)
            self.eta[:]= numpy.where(self.omega<0.0,self.eta_down,self.eta_up)
        elif self.learning_rule==3 or self.learning_rule==4: 
            # harel's hill function or  # biocyb omega and eta

            ##       p1: global magnitude
            ##       p4: basal calcium 
            ##       p3: exponent 
            ##       p2: half max 
            ##       want saturation, otherwise too much LTP 
            ##      
            ##       if p4 is 0, and Ca = p2 then we have half max 
            
              #omega=P.omega_offset+sig(Ca(c)-P.alpha2,P.beta2)-
                  #0.25*sig(Ca(c)-P.alpha1,P.beta1);

            # shouldn't it be with .5? 
            self.omega[:]=self.omega_offset+(sig(self.Ca-self.alpha2,self.beta2)-
                                      0.5*sig(self.Ca-self.alpha1,self.beta1))
            self.eta[:]=self.eta_gamma0*self.Ca

            if self.learning_rule==3: # harels 
                self.eta[:]=self.p1*(self.Ca+self.p4)**self.p3/((self.Ca+self.p4)**self.p3+self.p2**self.p3)
            else:                     # biocyb 
                self.eta[:]=(self.p2+self.Ca**P.p3)/(self.p1+self.p4*(self.p2+self.Ca**self.p3))
                
        else:
            raise ValueError, "Unknown learning rule"
        
        self.weights+=self.eta*(self.omega-self.lambda_decay*self.weights)
        
        
        # decay the backspike, NMDA currents, etc.
        
        self.I_nmda_slow[:]-=self.I_nmda_slow/self.tau_nmda_s
        self.I_nmda_fast[:]-=self.I_nmda_fast/self.tau_nmda_f
        
        self.v_backspike_slow[:]-=self.v_backspike_slow/self.tau_backspike_slow
        self.v_backspike_fast[:]-=self.v_backspike_fast/self.tau_backspike_fast

    # replacing the update method
    #update=wrapAsMethod(update_methods.Calcium_update)
