import sys

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

from tools import *

def test_default():
    pattern=[ [100,200,300],[150,250] ]
   
    
    n1=Spike_Pattern(2,pattern)

    total_time=10000
    sim_params=default_sim_params(total_time)

    n1.save_spike_range=[0,total_time]
    

    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])

    

if __name__=="__main__":
        
    pattern=[ [100,200,300],[150,250] ]
   
    
    n1=Spike_Pattern(2,pattern)

    total_time=10000
    sim_params=default_sim_params(total_time)

    n1.save_spike_range=[0,total_time]
    

    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])

