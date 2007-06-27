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
    n1=Constant_Poisson(10,rate=20)
    n2=Izh(1)
    c=Constant_Connection(n1,n2,initial_weight_range=[1.0,1.0])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['epsp1','epsp2','V'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Izh')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(n2.saved_vars['t'])
    e1=numpy.array(n2.saved_vars['epsp1'])
    e2=numpy.array(n2.saved_vars['epsp2'])
    V=numpy.array(n2.saved_vars['V'])

    
    plot(t,e1,t,e2)

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()
    


def test_bursting():
    n1=Constant_Poisson(5,rate=20)
    n2=Izh(1,'tonic bursting')
    c=Constant_Connection(n1,n2,initial_weight_range=[1.0,1.0])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['epsp1','epsp2','V'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Izh')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(n2.saved_vars['t'])
    e1=numpy.array(n2.saved_vars['epsp1'])
    e2=numpy.array(n2.saved_vars['epsp2'])
    V=numpy.array(n2.saved_vars['V'])

    
    plot(t,e1,t,e2)

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()
    



if __name__=="__main__":
    n1=Constant_Poisson(5,rate=20)
    n2=Izh(1,'tonic bursting')
    c=Constant_Connection(n1,n2,initial_weight_range=[1.0,1.0])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['epsp1','epsp2','V'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Izh')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(n2.saved_vars['t'])
    e1=numpy.array(n2.saved_vars['epsp1'])
    e2=numpy.array(n2.saved_vars['epsp2'])
    V=numpy.array(n2.saved_vars['V'])

    
    plot(t,e1,t,e2)

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()
