import sys

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

from tools import *



def test_1D():
    rates=array([10.0,60.0])
    total_time=10000
    sim_params=default_sim_params(total_time)
    
    n1=Variable_Poisson(1,rate=rates,time_between_rates=2000)
    n1.save_spike_range=[0,total_time]
    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print_rates(n1,[0,2000])
    print "-----"
    print_rates(n1,[2000,4000])


def test_2D():
    pat1=[80.0,40.0]
    pat2=[10.0,20.0]
    rates=array([ pat1,pat2 ])
    
    total_time=10000
    sim_params=default_sim_params(total_time)
    
    n1=Variable_Poisson(2,rate=rates,time_between_rates=2000)
    n1.save_spike_range=[0,total_time]
    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print pat1
    print_rates(n1,[0,2000])
    print "-----"
    print pat2
    print_rates(n1,[2000,4000])
    
def test_5D():
    
    # rows - number of rates
    # columns - number of inputs
    
    N=5
    rates=eye(N)*30.0+20.0
    
    total_time=10000
    sim_params=default_sim_params(total_time)
    
    n1=Variable_Poisson(N,rate=rates,time_between_rates=2000)
    n1.save_spike_range=[0,total_time]
    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print rates[0]
    print_rates(n1,[0,2000])
    print "-----"
    print rates[1]
    print_rates(n1,[2000,4000])
    print "-----"
    print rates[2]
    print_rates(n1,[4000,6000])
    

    

    
if __name__=="__main__":
    
    # rows - number of rates
    # columns - number of inputs
    
    N=5
    rates=eye(N)*30.0+20.0
    
    total_time=10000
    sim_params=default_sim_params(total_time)
    
    n1=Variable_Poisson(N,rate=rates,time_between_rates=2000)
    n1.save_spike_range=[0,total_time]
    run_sim(sim_params,[n1],[])
    
    clf()
    plot_spikes([n1])
    print rates[0]
    print_rates(n1,[0,2000])
    print "-----"
    print rates[1]
    print_rates(n1,[2000,4000])
    print "-----"
    print rates[2]
    print_rates(n1,[4000,6000])
    
    
