cdef struct Neuron_Group_struct:
    int quantity
    double dt
    char *spike,*old_spike
    int number_of_connection_groups_to
    int type
    double save_spike_range[2]
    
    char *c_old_spike[31]
    int c_num_incell[31],c_num_outcell[31]
    int c_sign[31]
    double *c_g[31],*c_weights[31]
    
    double *V
    
    # constant fixed
    double *time_to_next_spike
    double min_latency,max_latency
    
    # variable poisson
    double time_between_rates,time_to_next_rate
    double *rate    
    int number_of_rates
    int sequential
    int which_rate

    
    # integrate and fire
    double *time_to_stop_refract,*V_reset
    double tau_m,tau_ex,tau_in
    double reset_adaptation,t_refract
    double V_rev_exc,V_rev_inh,V_thresh
    double g_exc_max,g_inh_max
    
    # stochastic rate
    double *y
    double activation_magnitude,tau_activation

    # spike pattern
    double time_to_next_pattern,time_for_last_pattern
    int *pattern_count,*pattern_length
    double rate_val

    # Izh
    double *u,*epsp1,*epsp2,*ipsp1,*ipsp2
    double a,b,c,d,I,V_peak,tau_epsp1,tau_epsp2,epsp_scale
    
    
    
cdef copy_spikes(Neuron_Group_struct *s):
    
    cdef int i
    
    for i from 0<=i<s.quantity:
        s.old_spike[i]=s.spike[i]
    
#    cdef char *tmp
#    tmp=s.spike
#    s.spike=s.old_spike
#    s.old_spike=tmp
    
cdef Neuron_Group_struct init_Neuron_Group(object self):
    cdef Neuron_Group_struct s
    
    s.type=-1
    s.dt=self.dt
    s.quantity=self.quantity
    s.spike=CharData(self.spike)
    s.old_spike=CharData(self.old_spike)
    s.number_of_connection_groups_to=len(self.connections_to)
    s.save_spike_range[0]=self.save_spike_range[0]
    s.save_spike_range[1]=self.save_spike_range[1]
    s.V=DoubleData(self.V)
    
    cdef int cg
    
    for cg from 0<=cg<s.number_of_connection_groups_to:
        s.c_old_spike[cg]=CharData(self.connections_to[cg].incell.old_spike)
        s.c_num_incell[cg]=self.connections_to[cg].incell.quantity
        s.c_num_outcell[cg]=self.connections_to[cg].outcell.quantity
        
        s.c_sign[cg]=self.connections_to[cg].sign
        s.c_g[cg]=DoubleData(self.connections_to[cg].g)
        s.c_weights[cg]=DoubleData(self.connections_to[cg].weights)

    
    
    if self.__module__=='splikes.neuron_groups.Silent_Neuron':
        s.type=0
    elif self.__module__=='splikes.neuron_groups.Constant_Fixed':
        s.type=1
        s.time_to_next_spike=DoubleData(self.time_to_next_spike)
        s.min_latency=self.min_latency
        s.max_latency=self.max_latency
        
    elif self.__module__=='splikes.neuron_groups.Constant_Poisson':
        s.type=2
    elif self.__module__=='splikes.neuron_groups.Variable_Poisson':
        s.type=3
        s.time_between_rates=self.time_between_rates
        s.time_to_next_rate=self.time_to_next_rate
        s.rate=DoubleData(self.rate)
        s.number_of_rates=Dim0(self.rate)
        s.sequential=self.sequential
        s.which_rate=self.which_rate
    elif self.__module__=='splikes.neuron_groups.Integrate_and_Fire':
        s.type=4
        s.time_to_stop_refract=DoubleData(self.time_to_stop_refract)
        s.V_reset=DoubleData(self.V_reset)
        s.tau_m=self.tau_m
        s.tau_ex=self.tau_ex
        s.tau_in=self.tau_in
        s.reset_adaptation=self.reset_adaptation
        s.t_refract=self.t_refract
        s.V_rev_exc=self.V_rev_exc
        s.V_rev_inh=self.V_rev_inh
        s.V_thresh=self.V_thresh
        s.g_exc_max=self.g_exc_max
        s.g_inh_max=self.g_inh_max

    elif self.__module__=='splikes.neuron_groups.Stochastic_Rate':
        s.type=5
        s.y=DoubleData(self.y)
        s.activation_magnitude=self.activation_magnitude
        s.tau_activation=self.tau_activation
    elif self.__module__=='splikes.neuron_groups.Spike_Pattern':
        s.type=6
        s.rate_val=self.rate
        s.time_to_next_spike=DoubleData(self.time_to_next_spike)
        s.pattern_count=IntData(self.pattern_count)
        s.pattern_length=IntData(self.pattern_length)
        s.time_to_next_pattern=self.time_to_next_pattern
        s.time_for_last_pattern=self.time_for_last_pattern
    elif self.__module__=='splikes.neuron_groups.Izh':
        s.type=7
        s.u=DoubleData(self.u)
        s.epsp1=DoubleData(self.epsp1)
        s.epsp2=DoubleData(self.epsp2)
        s.ipsp1=DoubleData(self.ipsp1)
        s.ipsp2=DoubleData(self.ipsp2)
        s.a=self.a
        s.b=self.b
        s.c=self.c
        s.d=self.d
        s.I=self.I
        s.tau_epsp1=self.tau_epsp1
        s.tau_epsp2=self.tau_epsp2
        s.V_peak=self.V_peak
        s.epsp_scale=self.epsp_scale
        
        
    else:
        print "Unimplemented Neuron_Group Update",self.__module__
    
    return s
        

cdef save_spikes(Neuron_Group_struct *s,object self,double t):
    
    if t>=s.save_spike_range[0] and t<=s.save_spike_range[1]:
        self.save_spikes(t)
    

cdef Silent_Neuron_update(Neuron_Group_struct *s,object self,double t):
    pass
    
cdef Constant_Fixed_update(Neuron_Group_struct *s,object self,double t):
    
    cdef int i
    cdef int qty
    qty=s.quantity
    
    cdef double *ttn
    ttn=s.time_to_next_spike
    
    cdef double dt
    dt=self.dt
    
    if t==0:
        for i from 0<=i<qty:
            ttn[i]=floor(randu()*(s.max_latency-s.min_latency)+s.min_latency)
            
    
    cdef char *spike
    spike=s.spike
    
    for i from 0<=i<qty:
        if t>=ttn[i]:
            spike[i]=1
            ttn[i]=ttn[i]+1000.0*dt/self.rate
        else:
            spike[i]=0
    
cdef Constant_Poisson_update(Neuron_Group_struct *s,object self,double t):
    
    cdef int i
    cdef int qty
    
    qty=s.quantity

    cdef double rate,dt
    
    rate=self.rate
    dt=s.dt
    
    cdef char *spike
    spike=s.spike
    
    
    for i from 0<=i<qty:
        if randu()<(rate/1000.0/dt):
            spike[i]=1
        else:
            spike[i]=0
             
    
cdef Variable_Poisson_update(Neuron_Group_struct *s,object self,double t):

    cdef int i,offset
    cdef int qty,number_of_rates
    qty=s.quantity

    cdef double dt
    dt=s.dt
    cdef double time_between_rates,time_to_next_rate
    time_between_rates=s.time_between_rates
    time_to_next_rate=s.time_to_next_rate
    
    cdef double *rate    
    rate=s.rate
    number_of_rates=s.number_of_rates

    cdef char *spike
    spike=s.spike
    
    cdef int sequential
    sequential=s.sequential
    
    if t>=time_to_next_rate:
        s.time_to_next_rate=s.time_to_next_rate+s.time_between_rates
        if s.sequential:
            s.which_rate=s.which_rate+1
            if s.which_rate>=number_of_rates:
                s.which_rate=0
        else:
            s.which_rate=<int> (randu()*(number_of_rates))
        
        
    offset=qty*s.which_rate
    for i from 0<=i<qty:
        if randu()<(rate[i+offset]/1000.0/dt):
            spike[i]=1
        else:
            spike[i]=0
    

cdef Integrate_and_Fire_update(Neuron_Group_struct *s,object self,double t):
    
    cdef int i,n,c,cg,offset
    cdef int qty,number_of_rates
    qty=s.quantity

    cdef double dt
    dt=s.dt

    cdef double *time_to_stop_refract,*V,*V_reset
    
    time_to_stop_refract=s.time_to_stop_refract
    V=s.V
    V_reset=s.V_reset
    
    cdef int number_of_connection_groups_to
    number_of_connection_groups_to=s.number_of_connection_groups_to
    
    cdef char *old_spike,*spike
    spike=s.spike
    
    cdef double gg_ex,gg_in,gmax
    cdef int sign
    
    cdef double tau_m,tau_in,tau_ex,reset_adaptation,t_refract
    tau_m=s.tau_m
    tau_ex=s.tau_ex
    tau_in=s.tau_in
    reset_adaptation=s.reset_adaptation
    t_refract=s.t_refract
    
    
    
    cdef double V_rev_exc,V_rev_inh,V_thresh
    V_rev_exc=s.V_rev_exc
    V_rev_inh=s.V_rev_inh
    V_thresh=s.V_thresh
    
    
    cdef double g_exc_max,g_inh_max
    g_exc_max=s.g_exc_max
    g_inh_max=s.g_inh_max
    
    cdef double *g,*weights,sum
    cdef int num_incell,num_outcell
    cdef int ni,no
    
    gg_ex=0.0
    gg_in=0.0
    for cg from 0<=cg<number_of_connection_groups_to:
        old_spike=s.c_old_spike[cg]
        num_incell=s.c_num_incell[cg]
        num_outcell=s.c_num_outcell[cg]
        
        sign=s.c_sign[cg]
        g=s.c_g[cg]
        weights=s.c_weights[cg]


        if sign>0:
            gmax=g_exc_max
        else:
            gmax=g_inh_max

        for ni from 0<=ni<num_incell:
            if old_spike[ni]:
                for no from 0<=no<num_outcell:
                    g[no+ni*num_outcell]=g[no+ni*num_outcell]+gmax*weights[no+ni*num_outcell]
                    
        sum=0.0
        for c from 0<=c<num_incell*num_outcell:
            sum=sum+g[c]
            
        if sign>0:
            gg_ex=gg_ex+sum
        else:
            gg_in=gg_in+sum
            
                    
            
    for n from 0<=n<qty:
        
        if time_to_stop_refract[n]<=t:
            V[n]=(V[n]-(V[n]-V_reset[n])*dt/tau_m+
                      dt*gg_ex*(V_rev_exc-V[n])/tau_m +
                      dt*gg_in*(V_rev_inh-V[n])/tau_m)

        if V[n]>V_thresh:
            spike[n]=1
            if reset_adaptation:
                V_reset[n]=V_reset[n]-adaptation_step
            V[n]=V_reset[n]
            time_to_stop_refract[n]=time_to_stop_refract[n]+t_refract
        else:
            spike[n]=0
            
    
    # decay conductances
    for cg from 0<=cg<number_of_connection_groups_to:
        num_incell=s.c_num_incell[cg]
        num_outcell=s.c_num_outcell[cg]
        
        sign=s.c_sign[cg]
        g=s.c_g[cg]

        for c from 0<=c<num_incell*num_outcell:
            if sign>0:
                g[c]=g[c]-g[c]*dt/tau_ex
            else:
                g[c]=g[c]-g[c]*dt/tau_in
                
            
    if reset_adaptation:
        for n from 0<=n<qty:
            if spike[n]:
                V_reset[n]=V_reset[n]-(V_reset[n]-V_rest)*dt/tau_adaptation
    
            
cdef Stochastic_Rate_update(Neuron_Group_struct *s,object self,double t):
    
    cdef int i,n,c,cg
    cdef int number_of_connection_groups_to
    number_of_connection_groups_to=s.number_of_connection_groups_to
    
    cdef char *old_spike,*spike
    spike=s.spike

    
    cdef double *y
    y=s.y
    
    cdef int sign
    cdef double *g,*weights,sum
    cdef int num_incell,num_outcell
    cdef int ni,no
    cdef int count
    cdef double dt
    dt=s.dt
    
    cdef int qty
    qty=s.quantity
    
    cdef double activation_magnitude,tau_activation
    activation_magnitude=s.activation_magnitude
    tau_activation=s.tau_activation
    
    for cg from 0<=cg<number_of_connection_groups_to:
        old_spike=s.c_old_spike[cg]
        num_incell=s.c_num_incell[cg]
        num_outcell=s.c_num_outcell[cg]
        
        sign=s.c_sign[cg]
        g=s.c_g[cg]
        weights=s.c_weights[cg]

        for ni from 0<=ni<num_incell:
            if old_spike[ni]:
                for no from 0<=no<num_outcell:
                    g[no+ni*num_outcell]=g[no+ni*num_outcell]+activation_magnitude*weights[no+ni*num_outcell]*sign

        for no from 0<=no<num_outcell:
            sum=0.0
            for ni from 0<=ni<num_incell:
                sum=sum+g[no+ni*num_outcell]
                
            y[no]=y[no]+(1.0/tau_activation)*(sum-y[no])*dt
                
            
        # decay conductances
        for c from 0<=c<num_incell*num_outcell:
            g[c]=g[c]-g[c]*dt/tau_activation

                
            
    for n from 0<=n<qty:
        if randu()<(y[n]/1000.0/dt):
            spike[n]=1
        else:
            spike[n]=0
            
    

cdef Spike_Pattern_update(Neuron_Group_struct *s,object self,double t):
    
    cdef int i
    cdef int qty
    qty=s.quantity
    
    cdef double rate
    rate=s.rate_val
    
    cdef double dt
    dt=s.dt
    
    cdef double *ttn
    ttn=s.time_to_next_spike
    
    
    cdef char *spike
    spike=s.spike

    cdef int *pattern_count,*pattern_length
    pattern_count=s.pattern_count
    pattern_length=s.pattern_length
    
    if t>=s.time_to_next_pattern:
        s.time_for_last_pattern=t
        s.time_to_next_pattern=s.time_to_next_pattern+1000*dt/self.rate
        for i from 0<=i<qty:
            ttn[i]=self.pattern[i][0]+t
            pattern_count[i]=0
        
    for i from 0<=i<qty:
        if t>=ttn[i]:
            spike[i]=1
            
            pattern_count[i]=pattern_count[i]+1
            if pattern_count[i]==pattern_length[i]: # the end
                ttn[i]=1e500
            else:
                ttn[i]=self.pattern[i][pattern_count[i]]+s.time_for_last_pattern
        else:
            spike[i]=0
    
cdef Izh_update(Neuron_Group_struct *s,object self,double t):
    
    cdef int i,n,c,cg,offset
    cdef int qty,number_of_rates
    qty=s.quantity

    cdef double dt
    dt=s.dt

    cdef double *V,*u,*epsp1,*epsp2,*ipsp1,*ipsp2
    
    V=s.V
    u=s.u
    epsp1=s.epsp1
    epsp2=s.epsp2
    ipsp1=s.ipsp1
    ipsp2=s.ipsp2
    
    cdef int number_of_connection_groups_to
    number_of_connection_groups_to=s.number_of_connection_groups_to

    cdef double a,b,d,I,V_reset
    a=s.a
    b=s.b
    V_reset=s.c
    d=s.d
    I=s.I

    cdef char *old_spike,*spike
    spike=s.spike

    for n from 0<=n<qty:  # not sure if this works here...
        if spike[n]:
            V[n]=V_reset
    
    cdef epsp_scale
    cdef int sign
    
    epsp_scale=s.epsp_scale
    
    cdef double V_peak
    V_peak=s.V_peak
    
    cdef double g_exc_max,g_inh_max
    g_exc_max=s.g_exc_max
    g_inh_max=s.g_inh_max
    
    cdef double *g,*weights,sum
    cdef int num_incell,num_outcell
    cdef int ni,no
    
    gg_ex=0.0
    gg_in=0.0
    
    cdef double tm,norm_const,w
    
    cdef double tau_epsp1,tau_epsp2
    tau_epsp1=s.tau_epsp1
    tau_epsp2=s.tau_epsp2
    
    tm=(log(tau_epsp2)-log(tau_epsp1))/(1.0/tau_epsp1-1.0/tau_epsp2)
    norm_const=exp(-tm/tau_epsp1)-exp(-tm/tau_epsp2)
      
    # decay first, so that i don't decay new epsps on first iteration 

    for n from 0<=n<qty:
    
        epsp1[n]=epsp1[n]-epsp1[n]*dt/tau_epsp1
        epsp2[n]=epsp2[n]-epsp2[n]*dt/tau_epsp2
    
        ipsp1[n]=ipsp1[n]-ipsp1[n]*dt/tau_epsp1
        ipsp2[n]=ipsp2[n]-ipsp2[n]*dt/tau_epsp2

    for cg from 0<=cg<number_of_connection_groups_to:

        old_spike=s.c_old_spike[cg]
        num_incell=s.c_num_incell[cg]
        num_outcell=s.c_num_outcell[cg]
        
        sign=s.c_sign[cg]
        weights=s.c_weights[cg]

        for ni from 0<=ni<num_incell:
            if old_spike[ni]:
                
                w=weights[no+ni*num_outcell]
                if sign>0:
                    for no from 0<=no<num_outcell:
                        epsp1[no]=epsp1[no]+epsp_scale*w/norm_const
                        epsp2[no]=epsp2[no]+epsp_scale*w/norm_const
                else:
                    for no from 0<=no<num_outcell:
                        ipsp1[no]=ipsp1[no]+epsp_scale*w/norm_const
                        ipsp2[no]=ipsp2[no]+epsp_scale*w/norm_const
                    
            
    for n from 0<=n<qty:
        
        V[n]=(V[n] +(0.04*V[n]*V[n]+5.0*V[n]+140.0-u[n])*dt+
                      ((epsp1[n]-epsp2[n])-(ipsp1[n]-ipsp2[n]))*dt
             )

        u[n]=u[n]+a*(b*V[n]-u[n])*dt

        if V[n]>V_peak:
            spike[n]=1
            V[n]=V_peak
        else:
            spike[n]=0
            
    

cdef nupdate(Neuron_Group_struct *s,object self,double t):
    
    if s.type==-1:
        self.update(t)
        return
    
    copy_spikes(s)

    if s.type==0: # Silent Neuron
        pass
    elif s.type==1: # Constant Fixed
        Constant_Fixed_update(s,self,t)
    elif s.type==2:
        Constant_Poisson_update(s,self,t)
    elif s.type==3:
        Variable_Poisson_update(s,self,t)
    elif s.type==4:
        Integrate_and_Fire_update(s,self,t)
    elif s.type==5:
        Stochastic_Rate_update(s,self,t)
    elif s.type==6:
        Spike_Pattern_update(s,self,t)
    elif s.type==7:
        Izh_update(s,self,t)
    else:
        raise TypeError,"Unknown type: %d" % (s.type)
        

        
