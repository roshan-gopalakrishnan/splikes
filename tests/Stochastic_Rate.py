import sys

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

from tools import *


def test_default():
    n1=Constant_Poisson(2,rate=2)
    n2=Stochastic_Rate(1)
    n2.tau_activation=1000
    n2.activation_magnitude=1
    

    c=Constant_Connection(n1,n2,initial_weight_range=[1,1])

    total_time=10000
    sim_params=default_sim_params(total_time)
    
    n1.save_spike_range=[0,total_time]
    n2.save_spike_range=[0,total_time]
    n2.savevar(['y'])
    c.savevar(['g'])

    run_sim(sim_params,[n1,n2],[c])

    clf()
    
    subplot(2,1,1)
    plot_spikes([n1,n2])
    title('Stochastic Rate')
    ylabel('neuron number')
    
    subplot(2,1,2)
    
    t=array(n2.saved_vars['t'])
    y=array(n2.saved_vars['y'])
    
    plot(t,y)
    xlabel('time (ms)')
    ylabel('instantaneous rate, y')
    show()



def test_linear_input_output():
    sim_params=default_sim_params(10001,array([10000.0]))

    n1=Constant_Poisson(1,rate=10)
    n2=Stochastic_Rate(1)
    n2.tau_activation=1000
    n2.activation_magnitude=1
    
    n2.savevar(['y'])

    c=Constant_Connection(n1,n2,initial_weight_range=[1.0,1.0])

    rate_mat=arange(2,50,5.0)

    y_mat=[]
    for rate in rate_mat:
    
        sys.stdout.write('.')
        sys.stdout.flush()
        n1.rate=rate
    
    
        run_sim(sim_params,[n1,n2],[c])
    
        y=array(n2.saved_vars['y'])
        y_mat.append(y)

        
        
    plot(rate_mat,y_mat,'o')
    plot(rate_mat,rate_mat,'m--',linewidth=3)
    xlabel('Input Rate (Hz)')
    ylabel('Postsynaptic Variable, y')
    
    
def test_linear_weights_output():
    
    sim_params=default_sim_params(10001,array([10000.0]))
    n1=Constant_Poisson(1,rate=20)
    n2=Stochastic_Rate(1)
    n2.tau_activation=1000
    n2.activation_magnitude=1
    
    n2.savevar(['y'])

    c=Constant_Connection(n1,n2,initial_weight_range=[1.0,1.0])

    weight_mat=arange(0.1,2,0.2)

    y_mat=[]
    for weight in weight_mat:
    
        sys.stdout.write('.')
        sys.stdout.flush()
        c.initial_weight_range=[weight,weight]
    
    
        run_sim(sim_params,[n1,n2],[c])
    
        y=array(n2.saved_vars['y'])
        y_mat.append(y)

        
        
    plot(weight_mat,y_mat,'o')
    plot(weight_mat,20*weight_mat,'m--',linewidth=3)
    xlabel('Synaptic Weight')
    ylabel('Postsynaptic Variable, y')

    
    
def test_mean_variance():
    mn=[]
    sd=[]
    
    mag=linspace(0.1,2,20)
    
    for m in mag:
        
        n1=Constant_Poisson(169,rate=20)
        n2=Stochastic_Rate(1)
        n2.tau_activation=10/m
        n2.activation_magnitude=m
        
    
        c=Constant_Connection(n1,n2,initial_weight_range=[1,1])
    
        total_time=10000
        sim_params=default_sim_params(total_time)
        
        n1.save_spike_range=[0,total_time]
        n2.save_spike_range=[0,total_time]
        n2.savevar(['y'])
        c.savevar(['g'])
    
        run_sim(sim_params,[n1,n2],[c])
    
        clf()
        
##        subplot(2,1,1)
##        plot_spikes([n1,n2])
##        title('Stochastic Rate')
##        ylabel('neuron number')
##        
##        subplot(2,1,2)
        
        t=array(n2.saved_vars['t'])
        y=array(n2.saved_vars['y'])
        
        plot(t,y)
        xlabel('time (ms)')
        ylabel('instantaneous rate, y')
        show()

        
        mn.append(mean(y[7000:]))
        sd.append(std(y[7000:])**2)

        
    clf()

    subplot(1,2,1)
    plot(mag,mn,'-o')
    
    subplot(1,2,2)
    plot(mag,sd,'-o')
    

    # variance = mag * 10 * (tauact*m)/1000
    # tau act = 1000 / mag  ==> rate = input rate * weight * tauact/1000*mag
    
    # reducing tau_activation --> faster convergence, higher variance
    # 
    
def test_time_to_convergence():
    
    n1=Constant_Poisson(169,rate=30)
    n2=Stochastic_Rate(1)
    n2.tau_activation=50
    n2.activation_magnitude=0.1
    

    c=Constant_Connection(n1,n2,initial_weight_range=[1,1])

    total_time=[200,500,1500]
    
    clf()
    for i,tt in enumerate(total_time):
    
        sim_params=default_sim_params(tt)
        
        n1.save_spike_range=[0,tt]
        n2.save_spike_range=[0,tt]
        n2.savevar(['y'])
        c.savevar(['g'])
    
        run_sim(sim_params,[n1,n2],[c])
    
        
        subplot(len(total_time),1,i+1)
        
        t=array(n2.saved_vars['t'])
        y=array(n2.saved_vars['y'])
    
        
        plot(t,y)
        xlabel('time (ms)')
        ylabel('inst. rate, y')
        show()

        
    

if __name__=="__main__":
    
    n1=Constant_Poisson(169,rate=30)
    n2=Stochastic_Rate(1)
    n2.tau_activation=50
    n2.activation_magnitude=0.1
    

    c=Constant_Connection(n1,n2,initial_weight_range=[1,1])

    total_time=[200,500,1500]
    
    clf()
    for i,tt in enumerate(total_time):
    
        sim_params=default_sim_params(tt)
        
        n1.save_spike_range=[0,tt]
        n2.save_spike_range=[0,tt]
        n2.savevar(['y'])
        c.savevar(['g'])
    
        run_sim(sim_params,[n1,n2],[c])
    
        
        subplot(len(total_time),1,i+1)
        
        t=array(n2.saved_vars['t'])
        y=array(n2.saved_vars['y'])
    
        
        plot(t,y)
        xlabel('time (ms)')
        ylabel('inst. rate, y')
        show()

        
