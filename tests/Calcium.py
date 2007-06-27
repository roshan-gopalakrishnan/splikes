import sys
import pdb

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

from tools import *

from numpy import *
import numpy

def calcium_reset(self):
    
    self.g_nmda=-.0015*numpy.ones((self.incell.quantity,self.outcell.quantity),numpy.float)
    

def pairing(meta=False):
    sec=1000
    min=60*sec
    hr=60*min

    n1=Constant_Fixed(1,rate=1,min_latency=50,max_latency=50)
    n2=Silent_Neuron(1)
    n2.rest_voltage=20
    
    c=Calcium(n1,n2,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1])
    c.savevar(['weights','v_total','g_nmda',
               'I_nmda_slow','I_nmda_fast','I_nmda','Ca'])

               
               
    c.eta_gamma0=0.002
    c.beta1=60
    c.beta2=20
    c.alpha1=.25
    c.alpha2=.4
    
    if meta:
        c.k_plus=9e-6
        c.k_minus=9e-8
    else:
        c.k_plus=0
        c.k_minus=0
    
    c.user_reset=calcium_reset
    
    
    total_time=50*sec
    sim_params=default_sim_params(total_time,[total_time-1])

    
    V_mat=linspace(-80,-10,36)
    
    w=[]
    for V in V_mat:
        n2.rest_voltage=V
        
        run_sim(sim_params,[n1,n2],[c])
        
        w.append(c.saved_vars['weights'][-1][0][0])
        sys.stdout.write('.')
        sys.stdout.flush()
        

    clf()
    plot(V_mat,w,'.-')
    
        
def stdp(mu=0.7,meta=False):
    
    sec=1000
    min=60*sec
    hr=60*min

    n1=Constant_Fixed(1,rate=1,min_latency=150,max_latency=150)
    n2=Constant_Fixed(1,rate=1,min_latency=50,max_latency=50)
    
    n1.save_spike_range=[0,10000]
    n2.save_spike_range=[0,10000]
    
    c=Calcium(n1,n2,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1])
    c.savevar(['weights','v_total','g_nmda',
               'I_nmda_slow','I_nmda_fast','I_nmda','Ca'])

    # pairing has a larger eta gamma0
    c.eta_gamma0=0.0002
    c.beta1=60
    c.beta2=20
    c.alpha1=.25
    c.alpha2=.4
    c.i_nmda_mu=mu
    
    if meta:
        c.k_plus=9e-6
        c.k_minus=9e-8
    else:
        c.k_plus=0
        c.k_minus=0
    
    
    total_time=100*sec
    
#    total_time=10000  # for debugging
    sim_params=default_sim_params(total_time,100)

    delta_t=linspace(-150,150,30)
#    delta_t=[10]  # for debugging
    
    
    w=[]
    clf()
    for dt in delta_t:
    
        n2.min_latency=dt+150
        n2.max_latency=dt+150
        
        run_sim(sim_params,[n1,n2],[c])
        
        w.append(c.saved_vars['weights'][-1][0][0])
        sys.stdout.write('.')
        sys.stdout.flush()
        
        plot([delta_t[0],delta_t[-1]],[0.5,0.5],'k--',linewidth=2)
        gca().hold(True)
        plot([0,0],[0,1],'k--',linewidth=2)
        plot(delta_t[0:len(w)],w,'.-',linewidth=2)
        gca().hold(False)
        show()

    clf()
    if False:
        plot_spikes([n1,n2])
        
    elif False:
        
        t=sim_params['times_to_save']
        Ca=c.saved_vars['Ca']
        
        plot(t,Ca,'.-')
    
    elif True:
        
        plot([delta_t[0],delta_t[-1]],[0.5,0.5],'k--',linewidth=2)
        gca().hold(True)
        plot([0,0],[0,1],'k--',linewidth=2)
        plot(delta_t[0:len(w)],w,'.-',linewidth=2)
        gca().hold(False)
        show()
        
    
def make_square(N=4,sz=100,rates=[10,40],display=False):
    
    idx=r_[0:sz]
    
    sq_sz=sz/N
    l=[]
    for i in range(N):
        
    
        r=ones(idx.shape)*rates[0]
        r[((sq_sz*i)<=idx) & (idx<(sq_sz*(i+1)))]=rates[1]
        
        l.append(r)
        
    a=array(l,numpy.float)
        
    if display:
        
        for r in a:
            plot(r)
    
    
    return a

def make_gaussian(N=4,sz=100,max_rate=30,sigma=10,display=False):
    
    
    
    centers=r_[0:N]*sz/N

    # reverse, so the shift works easier down below
    centers=sz-centers
    
    mid=sz/2
    idx=r_[0:sz]-mid
    
    g=exp(-idx**2/(2.0*sigma**2))*max_rate
    
    
    g=concatenate((g[mid:],g[0:mid]))
    l=[]
    for c in centers:
        
        r=concatenate((g[c:],g[0:c]))
        l.append(r)
        
    a=array(l,numpy.float)
        
    if display:
        
        for r in a:
            plot(r)
            
    
    
    
    return a

    
def plot_weights(n_list,c_list):
    
    
    c=c_list[0]
    w=c.saved_vars['weights'][-1]
    
    plot(w,'.-')
    gca().hold(False)
    draw()
    

    
def square(num_patterns=4,N=100):
    
    sec=1000
    min=60*sec
    hr=60*min

    # rows - number of rates
    # columns - number of inputs
    
    
    rates=make_square(num_patterns,N)
    
    total_time=2500000
#    total_time=1000
    
    n1=Variable_Poisson(N,rate=rates,time_between_rates=500,sequential=False)
    n1.save_spike_range=[0,1000]
    
    n2=Constant_Poisson(20,rate=10)
    n2.save_spike_range=[0,1000]
    
    
    n3=Integrate_and_Fire(1)
    n3.save_spike_range=[0,1000]
    n3.savevar(['V'])

    
    c2=Constant_Connection(n2,n3,sign=-1,initial_weight_range=[0.5,0.5])
    
    
    c=Calcium(n1,n3,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1e10])
    c.savevar(['weights','v_total','g_nmda',
               'I_nmda_slow','I_nmda_fast','I_nmda','Ca'])

    c.eta_gamma0=5e-4
    c.beta1=60
    c.beta2=20
    c.alpha1=.25
    c.alpha2=.4
    
    c.k_plus=9e-5
    c.k_minus=9e-7
    
    sim_params=default_sim_params(total_time,1000,
                        display=plot_weights,display_step=10000)

    
    run_sim(sim_params,[n1,n2,n3],[c,c2])
        
def gaussian(num_patterns=20,N=100):
    sec=1000
    min=60*sec
    hr=60*min

    # rows - number of rates
    # columns - number of inputs
    
    N=100
    rates=make_gaussian(num_patterns,N)
    
    total_time=2500000
#    total_time=1000
    
    n1=Variable_Poisson(N,rate=rates,time_between_rates=500,sequential=False)
    n1.save_spike_range=[0,1000]
    
    n2=Constant_Poisson(20,rate=10)
    n2.save_spike_range=[0,1000]
    
    
    n3=Integrate_and_Fire(1)
    n3.save_spike_range=[0,1000]
    n3.savevar(['V'])

    
    c2=Constant_Connection(n2,n3,sign=-1,initial_weight_range=[0.5,0.5])
    
    
    c=Calcium(n1,n3,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1e10])
    c.savevar(['weights','v_total','g_nmda',
               'I_nmda_slow','I_nmda_fast','I_nmda','Ca'])

    c.eta_gamma0=5e-4
    c.beta1=60
    c.beta2=20
    c.alpha1=.25
    c.alpha2=.4
    
    c.k_plus=9e-5
    c.k_minus=9e-7
    
    sim_params=default_sim_params(total_time,1000,
                        display=plot_weights,display_step=10000)

    
    run_sim(sim_params,[n1,n2,n3],[c,c2])
        

    
    
if __name__=="__main__":
        
    sec=1000
    min=60*sec
    hr=60*min

    # rows - number of rates
    # columns - number of inputs
    
    N=100
    rates=make_gaussian(20,N)
    
    total_time=2500000
#    total_time=1000
    
    n1=Variable_Poisson(N,rate=rates,time_between_rates=500,sequential=False)
    n1.save_spike_range=[0,1000]
    
    n2=Constant_Poisson(20,rate=10)
    n2.save_spike_range=[0,1000]
    
    
    n3=Integrate_and_Fire(1)
    n3.save_spike_range=[0,1000]
    n3.savevar(['V'])

    
    c2=Constant_Connection(n2,n3,sign=-1,initial_weight_range=[0.5,0.5])
    
    
    c=Calcium(n1,n3,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1e10])
    c.savevar(['weights','v_total','g_nmda',
               'I_nmda_slow','I_nmda_fast','I_nmda','Ca'])

    c.eta_gamma0=5e-4
    c.beta1=60
    c.beta2=20
    c.alpha1=.25
    c.alpha2=.4
    
    c.k_plus=9e-5
    c.k_minus=9e-7
    
    sim_params=default_sim_params(total_time,1000,
                        display=plot_weights,display_step=10000)

    
    run_sim(sim_params,[n1,n2,n3],[c,c2])
        

    
