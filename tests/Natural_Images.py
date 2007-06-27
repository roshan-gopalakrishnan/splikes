import sys

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

from tools import *

def test_default():
    n1=Natural_Images(169,time_between_rates=500)

    total_time=500
    sim_params=default_sim_params(total_time)
    n1.save_spike_range=[0,total_time]
    
    
    for i in range(20):
    
    
        run_sim(sim_params,[n1],[])
        plot(n1.X.copy(),'.-')
        
def test_default2():
    n1=Natural_Images(169,time_between_rates=500)

    total_time=500
    sim_params=default_sim_params(total_time)
    n1.save_spike_range=[0,total_time]
    
    
    for i in range(20):
    
    
        run_sim(sim_params,[n1],[])
        
        subplot(4,5,i+1)
        im=n1.plasticity.utils.weights2image(n1.params,[n1.X])
        
        pcolor(im,cmap=gray())
        gca().set_axis_bgcolor('k')
        axis('equal')


if __name__=="__main__":

    n1=Natural_Images(169,time_between_rates=500)

    total_time=500
    sim_params=default_sim_params(total_time)
    n1.save_spike_range=[0,total_time]
    
    
    for i in range(20):
    
    
        run_sim(sim_params,[n1],[])
        plot(n1.X.copy(),'.-')
        
    
    
