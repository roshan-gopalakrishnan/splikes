import numpy
import zpickle
from copy import deepcopy

import time

def wrapAsMethod(to_be_wrapped):

    def wrapped(*args, **kw):
        return to_be_wrapped(*args, **kw)
    
    return wrapped


class Neuron_Group(object):
    
    def __init__(self,qty=1):
        
        self.quantity=qty
        self.dt=1.0
        self.rest_voltage=-65
        
        self.debug=False
        
        self.spike=numpy.zeros(qty,numpy.bool)
        
        self.old_spike=numpy.zeros(qty,numpy.bool)
        self.V=numpy.ones(qty,numpy.float)*self.rest_voltage

        self.save_spike_range=[0,-1]
        self.saved_spikes={'t':[],'n':[]}
        
        self.min_spike_range=0.0
        self.max_spike_range=0.0

        self.saved_vars={'t':[]}
        self.save_var_times=[]
        self.save_count=0
        
        
        self.connections_from=[]
        self.connections_to=[]


    def update(self,t):
        pass

    def savevar(self,var):
        for v in var:
            self.saved_vars[v]=[]

            
    def save_spikes(self,t):
        if t>=self.save_spike_range[0] and t<=self.save_spike_range[1]:
            new_spikes=list(numpy.flatnonzero(self.spike))
            self.saved_spikes['n'].extend(new_spikes)
            self.saved_spikes['t'].extend([t]*len(new_spikes))

        
    def save(self,t):
        self.t=t
        for v in self.saved_vars:
            self.saved_vars[v].append(deepcopy(self.__getattribute__(v)))

        
    def swap_spikes(self):
        
        self.old_spike[:]=self.spike

#        self.spike,self.old_spike=self.old_spike,self.spike  # swap
        
    def _reset_(self):
        # reset for t=0
        
        self.min_spike_range=self.save_spike_range[0]
        self.max_spike_range=self.save_spike_range[1]
        
        for v in self.saved_vars.keys():
            self.saved_vars[v]=[]
            
        self.saved_spikes={'t':[],'n':[]}

        
        self.V=numpy.ones(self.quantity,numpy.float)*self.rest_voltage

