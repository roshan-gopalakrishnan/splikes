import sys
from Waitbar import Waitbar
import numpy
from copy import deepcopy

import time
import simutils

import zpickle


def save_sim(sim_params,N,C):
    
    sim_params=deepcopy(sim_params)
    sim_params['display']=None
    d={'sim_params':sim_params,'N':N,'C':C}
    
    
    zpickle.save(d,sim_params['save_sim_file'])
    
    
def load_sim(fname):
    
    d=zpickle.load(fname)

    return d['sim_params'],d['N'],d['C']

def default_sim_params(total_time=3000,times_to_save=[],
                        display=None,
                        display_step=None,
                        display_hash=False,
                        hash_step=None,
                        save_sim_file='untitled.dat'):
    
    sim_params={}
    sim_params['display']=display
    if display_step:
        sim_params['display_step']=display_step
    else:
        sim_params['display_step']=total_time/1000

    sim_params['total_time']=total_time
    
    sim_params['display_hash']=display_hash
    if hash_step:
        sim_params['hash_step']=hash_step
    else:
        sim_params['hash_step']=total_time/1000
        
    if sim_params['display_step']==0:
        sim_params['display_step']=1
    if sim_params['hash_step']==0:
        sim_params['hash_step']=1
        
    sim_params['save_sim_file']=save_sim_file
    
    if not times_to_save:
        T=sim_params['total_time']
        step=1.0
        sim_params['times_to_save']=numpy.arange(0,T,step,numpy.float)
    else:
        if isinstance(times_to_save,int):  # gave an integer
            T=sim_params['total_time']
            sim_params['times_to_save']=numpy.linspace(0,T-1,times_to_save)
        else:
            sim_params['times_to_save']=numpy.array(times_to_save,numpy.float)
    
    return sim_params
    


def run_sim(sim_params,n_list,c_list):
    
    wb = Waitbar(False)

    save_count=0
    num_iter=sim_params['total_time']
    
    hash_step=sim_params['hash_step']
    display_step=sim_params['display_step']
    
        
    # reset all of them
    for n in n_list:
        n._reset_()
    for c in c_list:
        c._reset_()
        

    
    display_hash=sim_params['display_hash']
    times_to_save=sim_params['times_to_save']
    
    for t in range(num_iter):
        
        for n in n_list:
            n.update(t)
        for c in c_list:
            c.update(t)
            c.saturate()
            
            
        for n in n_list:
            n.save_spikes(t)
        
        if save_count<len(times_to_save):
            
            if t>=times_to_save[save_count]:
                for n in n_list:
                    n.save(t)
                for c in c_list:
                    c.save(t)

                if save_count % display_step == 0:
                    if sim_params['display']:
                        sim_params['display'](n_list,c_list)
        
                save_count=save_count+1
        
        if display_hash:
            if t % hash_step == 0:
                wb.updated(t/float(num_iter))

                
# to speed it up                
run_sim=simutils.run_sim

run_sim_orig=simutils.run_sim_working

