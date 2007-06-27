import sys
import pdb
from copy import deepcopy
import time
import os

try:
    import splikes
except ImportError:
    sys.path.append('../../../splikes')

    
from splikes import *
from splikes.simplot import *

from tools import *

from numpy import *
import numpy


import sys
sys.path.append('/home/bblais/python/work/plasticity')

import plasticity

def splikes2plasticity(splikes_weights):

    weights=[numpy.zeros((1,prod(splikes_weights.shape)),numpy.float)]
    weights[0][:]=splikes_weights.ravel()
    
    return weights

def drawit(N,C):
    
    c=C[0]
    n1=N[0]
    
    th=c.saved_vars['th']
    t=c.saved_vars['t']
    
    params=n1.params
    weights=splikes2plasticity(c.weights)
    
    im=plasticity.utils.weights2image(params,weights)

    
    clf()
    subplot(1,2,1)
    pcolor(im,cmap=cm.gray)
    axis('scaled')
    title("[%.4e, %.4e]" % (min(c.weights),max(c.weights)))
    
    
    
##    mask=N[0].params['pattern_input'][0]['mask'].astype(int)
##    idx_on=nonzero(mask.ravel())
##    idx_off=nonzero((1-mask).ravel())
##    
##    weights[idx_off]=weights[idx_on].min()
##    
##    weights.shape=13,13
##    weights=flipud(weights)
##    im=weights
##
###    weights=[C[0].weights.T]
###    im=plasticity.utils.weights2image(params,weights)
##    clf()
##    subplot(1,2,1)
##    pcolor(im,cmap=cm.gray)
##    axis('scaled')
##    colorbar()
    
    subplot(1,2,2)
    plot(t,th,'.-')

    show()

def test_bcmlike_1D():
    
    rates=array([20.0,40.0])
    total_time=10000000
    
    n1=Variable_Poisson(1,rate=rates,time_between_rates=total_time/2)
    n1.save_spike_range=[0,total_time]
    n2=Stochastic_Rate(1)
    n2.activation_magnitude=10
    n2.tau_activation=100
    
    c=Spiking_Rate(n1,n2,initial_weight_range=[0.5,0.5])
    c.weights_min=0.0
    c.weights_max=numpy.inf
    c.activation_magnitude=1
    c.tau_activation=100
    c.eta=2e-6
    c.tau_thresh = 10000
    c.thresh_o = 2
    c.learning_rule=2  # charlie rule
    
        
    c.savevar(['x','y','th','weights'])
    
    sim_params=default_sim_params(total_time,1000,display_hash=True)
    
    run_sim(sim_params,[n1,n2],[c])
    
    t=array(c.saved_vars['t'])
    w=squeeze(array(c.saved_vars['weights']))

    
    x=array(c.saved_vars['x'])
    y=array(c.saved_vars['y'])
    th=array(c.saved_vars['th'])
    
    subplot(1,2,1)
    plot(t/1000.0/60.0,w,linewidth=2)
    plot([total_time/1000.0/60.0/2,total_time/1000/60.0/2],[0,1.2],'k--')
    ylabel('Weights')
    xlabel('time (min)')
    text(20,.75,'20 Hz Input')
    text(100,.75,'40 Hz Input')
    subplot(1,2,2)
    plot(t/1000/60.0,y,linewidth=2)
    plot([total_time/1000/60.0/2,total_time/1000/60.0/2],[0,10],'k--')
    ylabel('Postsynaptic Variable, y')
    xlabel('time (min)')
    text(20,6,'20 Hz Input')
    text(100,6,'40 Hz Input')
        
    return (n1,n2,c)

def test_bcmlike_2D():
     
    # rows - number of rates
    # columns - number of inputs
    
    
    pat1=[50.0,20.0]
    pat2=[20.0,50.0]
    rates=array([ pat1,pat2 ])
    
    
    n1=Variable_Poisson(2,rate=rates,time_between_rates=500)
    n2=Stochastic_Rate(1)
    n2.activation_magnitude=10
    n2.tau_activation=100
    
    c=Spiking_Rate(n1,n2,initial_weight_range=[0.5,0.5])
    c.weights_min=0.0
    c.weights_max=numpy.inf
    c.activation_magnitude=1
    c.tau_activation=100
    c.eta=2e-6
    c.tau_thresh = 10000
    c.thresh_o = 1
    c.learning_rule=2  # charlie rule
    
        
    c.savevar(['x','y','th','weights'])
    
    total_time=1e7
    sim_params=default_sim_params(total_time,100,display_hash=True)
    
    run_sim(sim_params,[n1,n2],[c])
    
    t=array(c.saved_vars['t'])

    weights=deepcopy(c.saved_vars['weights'])
    
    # testing
    
    c.eta=0.0  # stop learning
    y_pat1=[]
    y_pat2=[]
    sim_params=default_sim_params(total_time/1000,100)

    for w in weights:
        sys.stdout.write('.')
        sys.stdout.flush()
        
        c.initial_weights=w

        n1.rate=rates[0,:]  # pat1
        run_sim(sim_params,[n1,n2],[c])
        y_pat1.append(c.saved_vars['y'][-1])
        
        n1.rate=rates[1,:]  # pat2
        run_sim(sim_params,[n1,n2],[c])
        y_pat2.append(c.saved_vars['y'][-1])
    
    print
    
    
    
    w=squeeze(array(weights))
    subplot(1,2,1)
    for i in range(w.shape[1]):
        plot(t/1000.0/60.0,w[:,i],linewidth=2)
        
    ylabel('Weights')
    xlabel('Time (min)')

    subplot(1,2,2)
    plot(t/1000.0/60.0,y_pat1,'r-',linewidth=2)
    plot(t/1000.0/60.0,y_pat2,'m-',linewidth=2)
    xlabel('Time (min)')
    ylabel('Response')
    
    
    return (n1,n2,c,y_pat1,y_pat2,w)

def test_bcmlike_Nd():
    
    # rows - number of rates
    # columns - number of inputs
    
    N=5
    rates=numpy.eye(N)*30.0+20.0
    
    total_time=100000000
               
    n1=Variable_Poisson(N,rate=rates,time_between_rates=500)
    n2=Stochastic_Rate(1)
    n2.activation_magnitude=10
    n2.tau_activation=100
    n2.savevar(['y'])
    
    c=Spiking_Rate(n1,n2,initial_weight_range=[0.5,0.5])
#    c.weight_saturation=[-numpy.inf,numpy.inf]
    c.activation_magnitude=1
    c.tau_activation=100
    c.eta=2e-6
    c.tau_thresh = 10000
    c.thresh_o = 10
    c.learning_rule=2  # charlie rule
    
        
    c.savevar(['x','y','th','weights'])

    sim_params=default_sim_params(total_time,100,display_hash=True)
    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
    print "Time = ",time.time()-t1
    
    t=array(c.saved_vars['t'])

    weights=deepcopy(c.saved_vars['weights'])
    
    # testing
    
    c=Constant_Connection(n1,n2,initial_weight_range=[0.5,0.5])
    y_pats=[ [] for i in range(N)]
    T=10000
    sim_params=default_sim_params(T)
    

    print "weights[0]: ",weights[0]
    weights2=[]
    count=0
    for w in weights:
        print "%d/%d" % (count,len(weights))
        sys.stdout.flush()
        count=count+1
        
        weights2.append(w)
        c.initial_weights=w

        for i in range(N):
        
            n1.rate=rates[i,:]  # pat1
            run_sim(sim_params,[n1,n2],[c])
            y_pats[i].append(n2.saved_vars['y'][-1])
        
    
    
    
    
    w=squeeze(array(weights2))
    subplot(1,2,1)
    for i in range(w.shape[1]):
        plot(t/1000.0/60.0,w[:,i],linewidth=2)
        
    subplot(1,2,2)
    for i in range(N):
        plot(t/1000.0/60.0,y_pats[i],linewidth=2)
    

def test_natural_images():
    
    sfname='sim060107_1.dat'
    
    if not os.path.exists(sfname):
    
        # trying to get same params as doit110505.m
        
        n1=Natural_Images(169,time_between_rates=200,rate_range=[0.0,50.0])
        n2=Stochastic_Rate(1,activation_magnitude=10,tau_activation=100)
        
        c=Spiking_Rate(n1,n2,initial_weight_range=[0.5,0.5])
        c.weight_saturation=[-inf,inf]
        c.tau_thresh=10000.0
        c.tau_activation=100.0
        c.thresh_o=10.0
        c.eta=1e-7
        c.learning_rule=2  # 1 - bcm, 2 - lawcooper, 3 - hebb
        
        million=1e6
        billion=1e9
        trillion=1e12
        
#        total_time=2.5*billion  # same as doit110505.m in soctagon
        total_time=250*million
        
        #total_time=25000
        sim_params=default_sim_params(total_time,10000,
                            display=drawit,
                            display_hash=True,
                            display_step=total_time/10000,
                            hash_step=total_time/10000,
                            save_sim_file=sfname)
        
        n1.save_spike_range=[0,10000]
        n2.save_spike_range=[0,10000]
        c.savevar(['weights','th'])
    
        
        t1=time.time()
        run_sim(sim_params,[n1,n2],[c])
        print "Time = ",time.time()-t1
        
    
        del n1.plasticity
        for p in n1.params['pattern_input']:
            p['var']=None
        
        save_sim(sim_params,[n1,n2],[c])
        

        
    sim_params,N,C=load_sim(sfname)

    
    drawit(N,C)
    
    
def scale_shift_filter(image,scale=1.0,shift=0.0,truncate=False):

    
    mn=[]
    mx=[]
    for im in image:

        im[:]=scale*im+shift

        if truncate:
            im[im<0]=0.0
            
        mn.append(im.min())
        mx.append(im.max())

    print "Min: ",min(mn)
    print "Max: ",max(mx)

    
def test_natimg_ONOFF():
    sfname='sim060707_1.dat'
    
    # trying to get same params as doit110505.m, but with 2 eyes
    
    n1=Natural_Images([169,169],
            filter=[scale_shift_filter,scale_shift_filter],
            filter_params=[{'scale':4.0,'shift':25.0,'truncate':True},
                           {'scale':-4.0,'shift':25.0,'truncate':True}],
            time_between_rates=500)
    n2=Stochastic_Rate(1,activation_magnitude=10,tau_activation=100)
    
    c=Spiking_Rate(n1,n2,initial_weight_range=[0.5,0.5])
    c.weight_saturation=[-inf,inf]
    c.tau_thresh=10000.0
    c.tau_beta=10000.0
    c.tau_activation=100.0
    c.thresh_o=10.0
    c.eta=5e-8
    c.learning_rule=2  # 1 - bcm, 2 - lawcooper, 3 - hebb
    
    million=1e6
    billion=1e9
    trillion=1e12
    
#        total_time=2.5*billion  # same as doit110505.m in soctagon
    total_time=250*million
    
    #total_time=25000
    sim_params=default_sim_params(total_time,10000,
                        display=drawit,
                        display_hash=True,
                        display_step=total_time/10000,
                        hash_step=total_time/10000,
                        save_sim_file=sfname)
    
    n1.save_spike_range=[0,10000]
    n2.save_spike_range=[0,10000]
    c.savevar(['weights','th'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
    print "Time = ",time.time()-t1
    

    del n1.plasticity
    for p in n1.params['pattern_input']:
        p['var']=None
    
    save_sim(sim_params,[n1,n2],[c])
        
    
    
if __name__=="__main__":

    
    sfname='sim060807_1.dat'
    
    
        
    n1=Natural_Images(169,time_between_rates=200,rate_range=[0.0,50.0])
    n2=Stochastic_Rate(1,activation_magnitude=1,tau_activation=100)
    
    c=Spiking_Rate(n1,n2,initial_weight_range=[0.5,0.5])
    c.weight_saturation=[-inf,inf]
    c.tau_thresh=10000.0
    activation_magnitude=1
    c.tau_activation=100.0
    c.thresh_o=10.0
    c.eta=1e-7
    c.learning_rule=2  # 1 - bcm, 2 - lawcooper, 3 - hebb
    
    million=1e6
    billion=1e9
    trillion=1e12
    
#        total_time=2.5*billion  # same as doit110505.m in soctagon
    total_time=250*million
    
    #total_time=25000
    sim_params=default_sim_params(total_time,10000,
                        display=drawit,
                        display_hash=True,
                        display_step=total_time/10000,
                        hash_step=total_time/10000,
                        save_sim_file=sfname)
    
    n1.save_spike_range=[0,10000]
    n2.save_spike_range=[0,10000]
    c.savevar(['weights','th'])

    
    t1=time.time()
    run_sim(sim_params,[n1,n2],[c])
    print "Time = ",time.time()-t1
    

    del n1.plasticity
    for p in n1.params['pattern_input']:
        p['var']=None
    
    save_sim(sim_params,[n1,n2],[c])
        

    
