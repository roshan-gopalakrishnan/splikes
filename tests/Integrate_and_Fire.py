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
        
    n1=Constant_Poisson(15,rate=200)
    n2=Integrate_and_Fire(1)
    c=Constant_Connection(n1,n2,initial_weight_range=[1.0,1.0])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['V'])
    c.savevar(['g'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Integrate and Fire')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(c.saved_vars['t'])
    g=numpy.array(c.saved_vars['g'])
    V=array(n2.saved_vars['V'])
    
    for i in range(g.shape[1]):
        plot(t,g[:,i])
    ylabel('conductance, g')

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()

def test_high_weights():
    n1=Constant_Poisson(15,rate=200)
    n2=Integrate_and_Fire(1)
    c=Constant_Connection(n1,n2,initial_weight_range=[2.0,2.0])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['V'])
    c.savevar(['g'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Integrate and Fire')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(c.saved_vars['t'])
    g=numpy.array(c.saved_vars['g'])
    V=array(n2.saved_vars['V'])
    
    for i in range(g.shape[1]):
        plot(t,g[:,i])
    ylabel('conductance, g')

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()
    

def test_low_weights():
    n1=Constant_Poisson(15,rate=200)
    n2=Integrate_and_Fire(1)
    c=Constant_Connection(n1,n2,initial_weight_range=[.5,0.5])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['V'])
    c.savevar(['g'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Integrate and Fire')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(c.saved_vars['t'])
    g=numpy.array(c.saved_vars['g'])
    V=array(n2.saved_vars['V'])
    
    for i in range(g.shape[1]):
        plot(t,g[:,i])
    ylabel('conductance, g')

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()
    
    
    
      
if __name__=="__main__":


    n1=Constant_Poisson(15,rate=200)
    n2=Integrate_and_Fire(1)
    c=Constant_Connection(n1,n2,initial_weight_range=[.1,4])
    

    total_time=200
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['V'])
    c.savevar(['g'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
   
    print "Time = ",time.time()-t1
    

    clf()
    
    subplot(3,1,1)
    plot_spikes([n1,n2])
    title('Integrate and Fire')
    ylabel('neuron number')
    
    subplot(3,1,2)
    
    t=numpy.array(c.saved_vars['t'])
    g=numpy.array(c.saved_vars['g'])
    V=array(n2.saved_vars['V'])
    
    for i in range(g.shape[1]):
        plot(t,g[:,i])
    ylabel('conductance, g')

    subplot(3,1,3)
    plot(t,V)
    ylabel('voltage, V')
    xlabel('time (ms)')
        
    show()
    
