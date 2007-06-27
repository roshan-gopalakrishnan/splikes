#search for 31 for maximum numbers
include "simutils.pxi"

include "c_connection_groups.pyx"

include "c_neuron_groups.pyx"


        
        
###===================== MAIN SIMULATION =======================###            
            
            
def run_sim(sim_params,n_list,c_list):
    
    cdef int time_to_save
    time_to_save=60*60*24
    
    wb = Waitbar(False)

    cdef long long t
    cdef int i

    cdef int ln,lc
    
    ln=len(n_list)
    lc=len(c_list)
    
    
    cdef long long num_iter
    num_iter=sim_params['total_time']
    
    
    for i from 0<=i<ln:
        n_list[i]._reset_()
    for i from 0<=i<lc:    
        c_list[i]._reset_()

    cdef long long hash_step
    hash_step=sim_params['hash_step']
    
    cdef long long display_step
    display_step=sim_params['display_step']
        
    cdef int display_hash
    display_hash=sim_params['display_hash']
    
    cdef int save_count,max_save_count
    cdef double *times_to_savep
    
    max_save_count=len(sim_params['times_to_save'])
    save_count=0
    
    if max_save_count>0:
        times_to_savep=DoubleData(sim_params['times_to_save'])
    else:
        times_to_savep=NULL

        
    cdef Connection_Group_struct Sc[31]
    for i from 0<=i<lc:    
        Sc[i]=init_Connection_Group(c_list[i])
    cdef Neuron_Group_struct Sn[31]
    for i from 0<=i<ln:    
        Sn[i]=init_Neuron_Group(n_list[i])
        
    t1=time.time()
    for t from 0<=t<num_iter:


        for i from 0<=i<ln:
            nupdate(&Sn[i],n_list[i],t)
            
        for i from 0<=i<lc:    
            cupdate(&Sc[i],c_list[i],t)
            saturate(&Sc[i])
            

        for i from 0<=i<ln:
            save_spikes(&Sn[i],n_list[i],t)
            
        if save_count<max_save_count:
            if t>=times_to_savep[save_count]:
                for i from 0<=i<ln:
                    n_list[i].save(t)
        
                for i from 0<=i<lc:
                    c_list[i].save(t)

                save_count=save_count+1
                
                t2=time.time()
                if (t2-t1)>(time_to_save):  # 20 minutes
                    print 'Preemptive save of ',sim_params['save_sim_file'],'...',
                    sys.stdout.flush()
                    sim.save_sim(sim_params,n_list,c_list)
                    print 'done.'
                    sys.stdout.flush()
                    t1=t2


                
            if t % display_step == 0:
                if sim_params['display']:
                    sim_params['display'](n_list,c_list)
                t2=time.time()
                if (t2-t1)>(time_to_save):  # 20 minutes
                    print 'Preemptive save of ',sim_params['save_sim_file'],'...',
                    sys.stdout.flush()
                    sim.save_sim(sim_params,n_list,c_list)
                    print 'done.'
                    sys.stdout.flush()
                    t1=t2
        
        
        if display_hash:
            if t % hash_step == 0:
                
                wb.updated(t/float(num_iter))
                t2=time.time()
                if (t2-t1)>(time_to_save):  # 20 minutes
                    print 'Preemptive save of ',sim_params['save_sim_file'],'...',
                    sys.stdout.flush()
                    sim.save_sim(sim_params,n_list,c_list)
                    print 'done.'
                    sys.stdout.flush()
                    t1=t2

                
                
def run_sim_working(sim_params,n_list,c_list):
    
    
    
    wb = Waitbar(False)

    cdef long long t
    cdef int i

    cdef int ln,lc
    
    ln=len(n_list)
    lc=len(c_list)
    
    cdef long long progress_step
    
    cdef long long num_iter
    num_iter=sim_params['total_time']
    
    # avoid divide by zero error
    if num_iter>=sim_params['num_display_hash']:
        progress_step=num_iter/sim_params['num_display_hash']
    else:
        progress_step=1
        
    
    for i from 0<=i<ln:
        n_list[i]._reset_()
    for i from 0<=i<lc:    
        c_list[i]._reset_()

    cdef int display_step
    display_step=sim_params['display_step']
        
    cdef int display_hash
    display_hash=sim_params['display_hash']
    
    cdef int save_count,max_save_count
    cdef double *times_to_savep
    
    max_save_count=len(sim_params['times_to_save'])
    save_count=0
    
    if max_save_count>0:
        times_to_savep=DoubleData(sim_params['times_to_save'])
    else:
        times_to_savep=NULL

        
        
    for t from 0<=t<num_iter:


        for i from 0<=i<ln:
            n_list[i].update(t)
            
        for i from 0<=i<lc:    
            c_list[i].update(t)
            c_list[i].saturate()
            

        for i from 0<=i<ln:
            n_list[i].save_spikes(t)
            
        if save_count<max_save_count:
            if t>=times_to_savep[save_count]:
                for i from 0<=i<ln:
                    n_list[i].save(t)
        
                for i from 0<=i<lc:
                    c_list[i].save(t)

                if save_count % display_step == 0:
                    if sim_params['display']:
                        sim_params['display'](n_list,c_list)
        
                save_count=save_count+1
        
        if display_hash:
            if t % progress_step == 0:
                
                wb.updated(t/float(num_iter))
                
        

                
                
    
