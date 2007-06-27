import numpy
import zpickle
from copy import deepcopy

from numpy import newaxis
import simutils

import time

from neuron_groups import Silent_Neuron  # default value

def sig(x,beta):
    
  return((numpy.tanh(beta*x/2.0)+1.0)/2.0)

def wrapAsMethod(to_be_wrapped):

    def wrapped(*args, **kw):
        return to_be_wrapped(*args, **kw)
    
    return wrapped


class Connection_Group(object):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),
                      sign=1,initial_weight_range=[0.0,1.0],
                      valid_weight_range=[0.0,1e500],
                      initial_weights=None):
        
        self.quantity=incell.quantity*outcell.quantity
        
        qty=self.quantity
        
        self.incell=incell
        self.outcell=outcell
        self.dt=1.0
        self.type=0
        
        incell.connections_from.append(self)
        outcell.connections_to.append(self)

        self.initial_weight_range=initial_weight_range
        self.initial_weights=initial_weights
        
        self.valid_weight_range=valid_weight_range
        
        self.saved_vars={'t':[]}
        self.save_var_times=[]
        self.save_count=0
        
        self.sign=sign
        self.g=numpy.zeros((incell.quantity,outcell.quantity),numpy.float)
        self.weights=numpy.random.uniform(initial_weight_range[0],
                    initial_weight_range[1],(incell.quantity,outcell.quantity))

        self.weight_saturation=numpy.array(valid_weight_range,numpy.float)
    
        
    def _reset_(self):
        if self.initial_weights is None:
            self.weights=numpy.random.uniform(self.initial_weight_range[0],
                                              self.initial_weight_range[1],
                                              (self.incell.quantity,
                                              self.outcell.quantity))
        
        else:
            self.weights=numpy.array(self.initial_weights).reshape(self.incell.quantity,self.outcell.quantity)

        self.g=numpy.zeros((self.incell.quantity,
                            self.outcell.quantity),numpy.float)
        

        
    def saturate(self):
        
        self.weights[self.weights<self.valid_weight_range[0]]=self.valid_weight_range[0]
        self.weights[self.weights>self.valid_weight_range[1]]=self.valid_weight_range[1]
        
        
    saturate=wrapAsMethod(simutils.saturate)
        
        
    def update(self,t):
        pass
        
    def savevar(self,var):
        for v in var:
            self.saved_vars[v]=[]
    
    def save(self,t):
        self.t=t
        for v in self.saved_vars:
            self.saved_vars[v].append(deepcopy(self.__getattribute__(v)))

            
class Constant_Connection(Connection_Group):
    
    pass

        
        
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
            
        self.x=numpy.zeros(incell.quantity,numpy.float)
        self.xtmp=numpy.zeros(incell.quantity,numpy.float)
        self.y=numpy.zeros(outcell.quantity,numpy.float)
        self.ytmp=numpy.zeros(outcell.quantity,numpy.float)
        
        self.th=numpy.ones(outcell.quantity,numpy.float)*0.1
        
    def _reset_(self):
        super(Spiking_Rate,self)._reset_()
        
        self.x=numpy.zeros(self.incell.quantity,numpy.float)
        self.xtmp=numpy.zeros(self.incell.quantity,numpy.float)
        self.y=numpy.zeros(self.outcell.quantity,numpy.float)
        self.ytmp=numpy.zeros(self.outcell.quantity,numpy.float)
        
        self.th=numpy.ones(self.outcell.quantity,numpy.float)*0.1
        
        
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
    update=wrapAsMethod(simutils.Spiking_Rate_update)

    
    
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
        
        
    def _reset(self):
        
        self.v_backspike_slow=numpy.zeros(self.outcell.quantity,numpy.float)
        self.v_backspike_fast=numpy.zeros(self.outcell.quantity,numpy.float)
        self.v_total=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.B=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.h=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.g_nmda=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.I_nmda=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.I_nmda_fast=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.I_nmda_slow=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.Ca=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.eta=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        self.omega=numpy.zeros((self.incell.quantity,self.outcell.quantity),numpy.float)
        
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
        self.v_backspike_slow[:]=where(self.outcell.spike,self.peak_backspike_slow,self.v_backspike_slow)
        self.v_backspike_fast[:]=where(self.outcell.spike,self.peak_backspike_fast,self.v_backspike_fast)

        self.v_total[:]*=0.0  # erase
        try:
            self.v_total[:]+=self.outcell.V
        except NameError:
            self.v_total[:]+=self.rest_voltage
            
        self.v_total[:]+=self.backspike_amplitude*(v_backspike_fast+v_backspike_slow)
        

        # meta plasticity
        
        if self.k_minus !=0 or self.k_plus!=0:  # do the metaplasticity

            ##      /* from yeung
            ##       tau_meta = 1.0/(k_minus*pow(v+65,2)+ k_plus);
            ##       g[i] = g[i] + t_scale/tau_meta * (g_t*k_plus*tau_meta - g[i]);
            ##      */

            g_nmda[:]+=self.dt*(self.k_plus*self.g_t-
             (self.k_plus+self.k_minus*(v_total-self.Vo)**self.Vp)*self.g_nmda)
            

             
             
        # NMDA dynamics
        

        self.I_nmda_slow[:]+=where(self.incell.spike,
                self.i_nmda_mu*(self.i_nmda_s-self.I_nmda_slow),0.0)
        self.I_nmda_fast[:]+=where(self.incell.spike,
                self.i_nmda_mu*(self.i_nmda_f-self.I_nmda_fast),0.0)
        
        
        self.B[:] = 1.0/(1.0 + (exp(self.mg1 * self.v_total) / self.mg2))
	    
        # this h is negative 
        self.h[:]= self.B*(self.v_total-self.v_reversal)

        # this g_nmda is negative 
        self.I_nmda[:]= self.g_nmda*(self.I_nmda_fast + self.I_nmda_slow) * self.h

        self.Ca[:]+=(self.I_nmda - self.Ca/self.tau_ca)*dt
        
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
            self.omega[:]=where(self.omega>self.omega_max,self.omega_max,self.omega)
            self.eta[:]= where(self.omega<0.0,self.eta_down,self.eta_up)
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
    update=wrapAsMethod(simutils.STDP_update)

    
class Gerstner06(Connection_Group):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),
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
        
        self.R1=numpy.zeros(incell.quantity,numpy.float)
        self.R2=numpy.zeros(incell.quantity,numpy.float)
        self.O1=numpy.zeros(outcell.quantity,numpy.float)
        self.O2=numpy.zeros(outcell.quantity,numpy.float)
        
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
                
