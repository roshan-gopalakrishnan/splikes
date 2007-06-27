import sys

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *
from numpy import *

from tools import *

if __name__=="__main__":
    
    M=mat('[1 1 ; 1 -1]')
    
    total_time=100000
    sim_params=default_sim_params(total_time)
    
    n1=Noise_Poisson(M=M,offset=[5,10],
                    time_between_rates=200,type='laplace',save_rates=True)
    n1.save_spike_range=[0,total_time]
    run_sim(sim_params,[n1],[])
    
    figure(1)

    clf()
    
    plot_spikes([n1])
    print_rates(n1,[0,2000])
    print "-----"
    print_rates(n1,[2000,4000])

    figure(2)
    clf()

    rates=squeeze(array(n1.saved_rates))

    l=rates.shape[1]
    
    for i in range(l):
        subplot(l,1,i+1)
        hist(rates[:,i],50)
        
        
    figure(3)
    clf()
    
    plot(rates[:,0],rates[:,1],'o')
    
    
    show()
    
