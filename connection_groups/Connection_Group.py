import numpy

import zpickle
from copy import deepcopy

from numpy import newaxis

import time
import sys


from neuron_groups import Silent_Neuron  # default value

def sig(x,beta):
    
  return((numpy.tanh(beta*x/2.0)+1.0)/2.0)

def wrapAsMethod(to_be_wrapped):

    def wrapped(*args, **kw):
        return to_be_wrapped(*args, **kw)
    
    return wrapped


class Connection_Group(object):
    
    def __init__(self,incell=Silent_Neuron(),
                      outcell=Silent_Neuron(),
                      sign=1,initial_weight_range=[0.0,1.0],
                      valid_weight_range=[0.0,1e500],
                      initial_weights=None):
        
        self.quantity=incell.quantity*outcell.quantity
        
        qty=self.quantity
        
        self.incell=incell
        self.outcell=outcell
        self.dt=1.0
        self.type=0
        
        incell.connections_from.append(self)
        outcell.connections_to.append(self)

        self.initial_weight_range=initial_weight_range
        self.initial_weights=initial_weights
        
        self.valid_weight_range=valid_weight_range
        
        self.saved_vars={'t':[]}
        self.save_var_times=[]
        self.save_count=0
        
        self.sign=sign
        self.g=numpy.zeros((incell.quantity,outcell.quantity),numpy.float)
        self.weights=numpy.random.uniform(initial_weight_range[0],
                    initial_weight_range[1],(incell.quantity,outcell.quantity))

        self.weight_saturation=numpy.array(valid_weight_range,numpy.float)
    
        
    def _reset_(self):
        if self.initial_weights is None:
            self.weights=numpy.random.uniform(self.initial_weight_range[0],
                                              self.initial_weight_range[1],
                                              (self.incell.quantity,
                                              self.outcell.quantity))
        
        else:
            self.weights=numpy.array(self.initial_weights).reshape(self.incell.quantity,self.outcell.quantity)

        self.g=numpy.zeros((self.incell.quantity,
                            self.outcell.quantity),numpy.float)
        

        self.weight_saturation=numpy.array(self.weight_saturation,float)
        
    def saturate(self):
        
        self.weights[self.weights<self.weight_saturation[0]]=self.weight_saturation[0]
        self.weights[self.weights>self.weight_saturation[1]]=self.weight_saturation[1]
        
        
#    saturate=wrapAsMethod(update_methods.saturate)
        
        
    def update(self,t):
        pass
        
    def savevar(self,var):
        for v in var:
            self.saved_vars[v]=[]
    
    def save(self,t):
        self.t=t
        for v in self.saved_vars:
            self.saved_vars[v].append(deepcopy(self.__getattribute__(v)))

