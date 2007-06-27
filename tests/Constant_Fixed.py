import sys

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *
from tools import *

def test_no_latency():
        
    n1=Constant_Fixed(10,rate=10,min_latency=50,max_latency=50)

    total_time=1000
    sim_params=default_sim_params(total_time)

    n1.save_spike_range=[0,total_time]
    

    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])

    print_rates(n1,total_time)
  
    

def test_latency():
    n1=Constant_Fixed(10,rate=10,min_latency=10,max_latency=90)

    total_time=1000
    sim_params=default_sim_params(total_time)

    n1.save_spike_range=[0,total_time]
    

    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print_rates(n1,total_time)
    
    
def test_default():

    n1=Constant_Fixed(10)

    total_time=1000
    sim_params=default_sim_params(total_time)

    n1.save_spike_range=[0,total_time]
    

    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print_rates(n1,total_time)
    

if __name__=="__main__":
    
    
    n1=Constant_Fixed(10)

    total_time=1000
    sim_params=default_sim_params(total_time)

    n1.save_spike_range=[0,total_time]
    

    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print_rates(n1,total_time)
    
