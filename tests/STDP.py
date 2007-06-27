import sys
import time
import numpy

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

def test_default():
    sec=1000
    min=60*sec
    hr=60*min

    
    n1=Constant_Fixed(1,rate=1,min_latency=150,max_latency=150)
    n2=Constant_Fixed(1,rate=1,min_latency=50,max_latency=50)
    
    c=STDP(n1,n2,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1])
    c.savevar(['weights'])

    total_time=100*sec
    sim_params=default_sim_params(total_time,2)

    delta_t=linspace(-150,150,50)
#    delta_t=[10]  # for debugging
    
    
    w=[]
    clf()
    for dt in delta_t:
    
        n2.min_latency=dt+150
        n2.max_latency=dt+150
        
        t1=time.time()
        run_sim(sim_params,[n1,n2],[c])
        
        print time.time()-t1
        
        w.append(c.saved_vars['weights'][-1][0][0])
        sys.stdout.write('.')
        sys.stdout.flush()
        
        plot([delta_t[0],delta_t[-1]],[0.5,0.5],'k--',linewidth=2)
        gca().hold(True)
        plot([0,0],[0,1],'k--',linewidth=2)
        plot(delta_t[0:len(w)],w,'.-',linewidth=2)
        gca().hold(False)
        show()

    
    if True:
        
        plot([delta_t[0],delta_t[-1]],[0.5,0.5],'k--',linewidth=2)
        gca().hold(True)
        plot([0,0],[0,1],'k--',linewidth=2)
        plot(delta_t[0:len(w)],w,'.-',linewidth=2)
        gca().hold(False)
        show()
    



if __name__=="__main__":
    
    
    sec=1000
    min=60*sec
    hr=60*min

    
    n1=Constant_Fixed(1,rate=1,min_latency=150,max_latency=150)
    n2=Constant_Fixed(1,rate=1,min_latency=50,max_latency=50)
    
    n1.save_spike_range=[0,10000]
    n2.save_spike_range=[0,10000]
    c=STDP(n1,n2,initial_weight_range=[0.5,0.5],
                    valid_weight_range=[0,1])
    c.savevar(['weights'])

    total_time=100*sec
    sim_params=default_sim_params(total_time,2)

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

    
    if True:
        
        plot([delta_t[0],delta_t[-1]],[0.5,0.5],'k--',linewidth=2)
        gca().hold(True)
        plot([0,0],[0,1],'k--',linewidth=2)
        plot(delta_t[0:len(w)],w,'.-',linewidth=2)
        gca().hold(False)
        show()
        
