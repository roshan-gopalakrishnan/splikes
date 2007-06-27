cdef struct Connection_Group_struct:
    int output_num,input_num
    double *weights,sign,dt
    char *outspike,*inspike
    char *old_outspike,*old_inspike
    double *in_V,*out_V
    int type
    double *weight_saturation
    
    # Spiking Rate
    double *x,*xtmp,*y,*ytmp,*th,*beta
    int use_beta
    double tau_beta
    double activation_magnitude,tau_activation,tau_thresh
    double thresh_o,xo,yo,eta,lambda_decay
    int learning_rule

    # calcium
    double *v_backspike_slow, *v_backspike_fast
    double *v_total
    
    double *g_nmda
    double *I_nmda,*I_nmda_slow,*I_nmda_fast
    double *Ca


    double peak_backspike_slow,peak_backspike_fast
    double backspike_amplitude
    double k_minus,k_plus,g_t,Vo,Vp
    double i_nmda_mu,i_nmda_s,i_nmda_f,tau_nmda_s,tau_nmda_f
    double mg1,mg2,tau_ca
    double omega_offset,alpha1,alpha2,beta1,beta2,theta_o,omega_max
    double eta_down,eta_up
    double eta_p1,eta_p2,eta_p3,eta_p4,eta_gamma0
    double tau_backspike_fast,tau_backspike_slow

    
    # STDP
    
    double tau_minus,tau_plus,a_minus,a_plus,g_max
    int nearest_neighbor
    double *M,*P
    

    # Gerstner06
    double tau_x,tau_y
    double A2_minus,A2_plus,A3_minus,A3_plus

    double *r1,*r2,*o1,*o2
    
    
    
cdef saturate(Connection_Group_struct *s):
    
    cdef int output_num,input_num
    cdef int count
    cdef int o,i

    cdef double *weights,*weight_saturation
    
    weights=s.weights
    output_num=s.output_num
    input_num=s.input_num
    weight_saturation=s.weight_saturation
    
    count=0    
    for i from 0<=i<input_num:
        for o from 0<=o<output_num:

            if weights[count]<weight_saturation[0]:
                weights[count]=weight_saturation[0]
                
            if weights[count]>weight_saturation[1]:
                weights[count]=weight_saturation[1]
                
            count=count+1
    
 

cdef double sig(double x,double beta):
    #   return(exp(beta*x)/(1.0+exp(beta*x))); 
    return ((tanh(beta*x/2.0)+1.0)/2.0)


    
cdef Connection_Group_struct init_Connection_Group(object self):
    cdef Connection_Group_struct s

    
    # this is here, so that we don't have to keep calling C-API when
    # we get constant parameters...this gives us about a factor of 10-20!
    
    s.type=-1
    s.sign=self.sign
    s.dt=self.dt
    s.weights=DoubleData(self.weights)
    s.outspike=CharData(self.outcell.spike)
    s.inspike=CharData(self.incell.spike)
    
    s.in_V=DoubleData(self.incell.V)
    s.out_V=DoubleData(self.outcell.V)
    s.old_outspike=CharData(self.outcell.old_spike)
    s.old_inspike=CharData(self.incell.old_spike)
    s.output_num=Dim0(self.outcell.spike)
    s.input_num=Dim0(self.incell.spike)
    s.weight_saturation=DoubleData(self.weight_saturation)
    
    if self.__module__=="splikes.connection_groups.Constant_Connection":
        s.type=0
    elif self.__module__=="splikes.connection_groups.Spiking_Rate":
        s.type=1

        s.activation_magnitude=self.activation_magnitude
        s.tau_activation=self.tau_activation
        s.tau_thresh=self.tau_thresh
        s.thresh_o=self.thresh_o
        s.learning_rule=self.learning_rule
        s.xo=self.xo
        s.yo=self.yo
        s.eta=self.eta
        s.use_beta=self.use_beta
        s.tau_beta=self.tau_beta
        s.lambda_decay=self.lambda_decay
    
        s.x=DoubleData(self.x)
        s.xtmp=DoubleData(self.xtmp)
        s.y=DoubleData(self.y)
        s.ytmp=DoubleData(self.ytmp)
        s.th=DoubleData(self.th)
        s.beta=DoubleData(self.beta)
        
    elif self.__module__=="splikes.connection_groups.Calcium":
        s.type=2
        
        s.v_backspike_slow=DoubleData(self.v_backspike_slow)
        s.v_backspike_fast=DoubleData(self.v_backspike_fast)
        s.v_total=DoubleData(self.v_total)

        s.g_nmda=DoubleData(self.g_nmda)
        s.I_nmda=DoubleData(self.I_nmda)
        s.I_nmda_slow=DoubleData(self.I_nmda_slow)
        s.I_nmda_fast=DoubleData(self.I_nmda_fast)
        s.Ca=DoubleData(self.Ca)


        s.peak_backspike_slow=self.peak_backspike_slow
        s.peak_backspike_fast=self.peak_backspike_fast
        s.backspike_amplitude=self.backspike_amplitude
        s.k_minus=self.k_minus 
        s.k_plus=self.k_plus
        s.g_t=self.g_t
        s.Vo=self.Vo
        s.Vp=self.Vp
        s.i_nmda_mu=self.i_nmda_mu
        s.i_nmda_s=self.i_nmda_s
        s.i_nmda_f=self.i_nmda_f
        s.tau_nmda_s=self.tau_nmda_s
        s.tau_nmda_f=self.tau_nmda_f
        s.mg1=self.mg1
        s.mg2=self.mg2
        s.g_t=self.g_t
        s.tau_ca=self.tau_ca
        s.learning_rule=self.learning_rule
        s.omega_offset=self.omega_offset
        s.alpha1=self.alpha1
        s.alpha2=self.alpha2
        s.beta1=self.beta1
        s.beta2=self.beta2
        s.theta_o=self.theta_o
        s.omega_max=self.omega_max
        s.eta_down=self.eta_down
        s.eta_up=self.eta_up
        s.eta_p1=self.eta_p1
        s.eta_p2=self.eta_p2
        s.eta_p3=self.eta_p3
        s.eta_p4=self.eta_p4
        s.eta_gamma0=self.eta_gamma0
        s.lambda_decay=self.lambda_decay
        s.tau_backspike_fast=self.tau_backspike_fast
        s.tau_backspike_slow=self.tau_backspike_slow
        
    elif self.__module__=="splikes.connection_groups.STDP":
        s.type=3
        s.tau_minus=self.tau_minus
        s.tau_plus=self.tau_plus
        s.a_minus=self.a_minus
        s.a_plus=self.a_plus
        s.g_max=self.g_max
        s.nearest_neighbor=self.nearest_neighbor
        s.M=DoubleData(self.M)
        s.P=DoubleData(self.P)
        
    elif self.__module__=="splikes.connection_groups.Gerstner06":
        s.type=4
        
        s.tau_minus=self.tau_minus
        s.tau_plus=self.tau_plus
        s.tau_x=self.tau_x
        s.tau_y=self.tau_y
        s.A2_minus=self.A2_minus
        s.A2_plus=self.A2_plus
        s.A3_minus=self.A3_minus
        s.A3_plus=self.A3_plus
        s.nearest_neighbor=self.nearest_neighbor

        s.r1=DoubleData(self.R1)
        s.r2=DoubleData(self.R2)
        s.o1=DoubleData(self.O1)
        s.o2=DoubleData(self.O2)
    else:
        print "Unimplemented Connection_Group Update",self.__module__
    
    return s

cdef Constant_Connection_update(Connection_Group_struct *s,object self,double t):

    pass

    
cdef Spiking_Rate_update(Connection_Group_struct *s,object self,double t):
    cdef int output_num,input_num
    cdef double *weights,*x,*xtmp,*y,*ytmp,*th,*beta
    cdef char *outspike,*inspike
    
    cdef double activation_magnitude,sign,dt,tau_activation,tau_thresh
    cdef double thresh_o,xo,yo,eta,lambda_decay,tau_beta
    cdef int learning_rule,use_beta
    cdef double tmp_y
    
    activation_magnitude=s.activation_magnitude
    sign=s.sign
    dt=s.dt
    tau_activation=s.tau_activation
    tau_thresh=s.tau_thresh
    thresh_o=s.thresh_o
    learning_rule=s.learning_rule
    xo=s.xo
    yo=s.yo
    eta=s.eta
    lambda_decay=s.lambda_decay
    tau_beta=s.tau_beta
    use_beta=s.use_beta
    
    
    x=s.x
    xtmp=s.xtmp
    y=s.y
    ytmp=s.ytmp
    th=s.th
    beta=s.beta
    
    weights=s.weights
    outspike=s.outspike
    inspike=s.inspike
    output_num=s.output_num
    input_num=s.input_num
    
    cdef int o,i
    
    for i from 0<=i<input_num:
        if inspike[i]:
            xtmp[i]=xtmp[i]+activation_magnitude*sign
        x[i]=x[i]+(1.0/tau_activation)*(xtmp[i]-x[i])*dt
        xtmp[i]=xtmp[i]-1.0/tau_activation*xtmp[i]*dt
        
            
    for o from 0<=o<output_num:
        if outspike[o]:
            ytmp[o]=ytmp[o]+activation_magnitude
        y[o]=y[o]+(1.0/tau_activation)*(ytmp[o]-y[o])*dt
        ytmp[o]=ytmp[o]-1.0/tau_activation*ytmp[o]*dt

    if use_beta:
        for o from 0<=o<output_num:
            beta[o]=beta[o]+(1.0/tau_beta)*(y[o]-beta[o])*dt
            
    
        
    if use_beta:
        for o from 0<=o<output_num:
            tmp_y=y[o]-beta[o]
            if tmp_y<0.0:
                tmp_y=0.0
            th[o]=th[o]+(1.0/tau_thresh)*((tmp_y-yo)*(tmp_y-yo)/thresh_o-tmp_y)

        count=0    
        for i from 0<=i<input_num:
            for o from 0<=o<output_num:
                
                tmp_y=y[o]-beta[o]
                if tmp_y<0.0:
                    tmp_y=0.0
                
                # matches plasticity, not soctagon
                if learning_rule==1: # BCM
                    weights[count]=weights[count]+eta*((x[i]-xo)*(tmp_y-yo)*((tmp_y-yo)-th[o])-lambda_decay*weights[count])
                elif learning_rule==2: # law and cooper
                    weights[count]=weights[count]+eta*((x[i]-xo)*(tmp_y-yo)*((tmp_y-yo)-th[o])/th[o]-lambda_decay*weights[count])
                elif learning_rule==3: # hebb w/oja
                    weights[count]=weights[count]+eta*((x[i]-xo)*(tmp_y-yo)-(tmp_y-yo)*(tmp_y-yo)*weights[count]-lambda_decay*weights[count])
                
                count=count+1
        
    else:
        for o from 0<=o<output_num:
            th[o]=th[o]+(1.0/tau_thresh)*((y[o]-yo)*(y[o]-yo)/thresh_o-th[o])
        
        count=0    
        for i from 0<=i<input_num:
            for o from 0<=o<output_num:
                
                # matches plasticity, not soctagon
                if learning_rule==1: # BCM
                    weights[count]=weights[count]+eta*((x[i]-xo)*(y[o]-yo)*((y[o]-yo)-th[o])-lambda_decay*weights[count])
                elif learning_rule==2: # law and cooper
                    weights[count]=weights[count]+eta*((x[i]-xo)*(y[o]-yo)*((y[o]-yo)-th[o])/th[o]-lambda_decay*weights[count])
                elif learning_rule==3: # hebb w/oja
                    weights[count]=weights[count]+eta*((x[i]-xo)*(y[o]-yo)-(y[o]-yo)*(y[o]-yo)*weights[count]-lambda_decay*weights[count])
                
                count=count+1
        

cdef Calcium_update(Connection_Group_struct *s,object self,double t):
    
    cdef double *v_backspike_slow, *v_backspike_fast
    cdef double *v_total
    cdef double B
    cdef double h,eta,omega
    cdef double *g_nmda
    cdef double *I_nmda,*I_nmda_slow,*I_nmda_fast
    cdef double *Ca
    
    v_backspike_slow=s.v_backspike_slow
    v_backspike_fast=s.v_backspike_fast
    v_total=s.v_total
    g_nmda=s.g_nmda
    I_nmda=s.I_nmda
    I_nmda_slow=s.I_nmda_slow
    I_nmda_fast=s.I_nmda_fast
    Ca=s.Ca
        

    cdef int output_num,input_num
    cdef double *weights
    cdef char *outspike,*inspike
    
    weights=s.weights
    outspike=s.outspike
    inspike=s.inspike
    output_num=s.output_num
    input_num=s.input_num
    
    
    ##    from yeung's code 11/03 
    ##  
    ##      if(is_meta)
    ##                meta_dyn(g, v, t, p, ic); 
    ##    nmda_dyn(i_nmda, exc_spikes, g, v, t, p, ic);
    ##    ca_dyn(ca, i_nmda, t, p, ic);
    ##    learn_ca(w, ca, t, p, ic); 
    ##  
    ##    
    ##    yeung does: input, then learning, then output.
    ##    I do: input and output, then learning
    ##    
    ##    she also does her metaplasticity before the rest of the learning. 
    ##    
    ##    bkspike_dyn(&v_bk, post_spike_train, t, p, ic);
    ##    v = v_rest + v_bk; 
    ##
    ##    and she just does rest + v_bk
    ##  

    cdef int o,i,c
    
    # set active backspikes to peak value
    for o from 0<=o<output_num:
        if outspike[o]:
            v_backspike_slow[o]=s.peak_backspike_slow
            v_backspike_fast[o]=s.peak_backspike_fast
            
    
    c=0    
    for i from 0<=i<input_num:    # erase v_total
        for o from 0<=o<output_num:
            v_total[c]=0.0
            c=c+1
            

    cdef double *V
    V=s.out_V # voltage of the out cell
    
    c=0    
    for i from 0<=i<input_num:    
        for o from 0<=o<output_num:
            v_total[c]=v_total[c]+V[o]
            v_total[c]=v_total[c]+(
                s.backspike_amplitude*(v_backspike_fast[o]+
                                       v_backspike_slow[o]))
            c=c+1
        
        

    # meta plasticity
    
    if s.k_minus !=0 or s.k_plus!=0:  # do the metaplasticity

        ##      /* from yeung
        ##       tau_meta = 1.0/(k_minus*pow(v+65,2)+ k_plus);
        ##       g[i] = g[i] + t_scale/tau_meta * (g_t*k_plus*tau_meta - g[i]);
        ##      */

        c=0    
        for i from 0<=i<input_num:    
            for o from 0<=o<output_num:
                
                g_nmda[c]=g_nmda[c]+s.dt*(s.k_plus*s.g_t-
                    (s.k_plus+s.k_minus*(v_total[c]-s.Vo)**s.Vp)*g_nmda[c])
                c=c+1
        

    # NMDA dynamics
    
    c=0    
    for i from 0<=i<input_num:

        if inspike[i]:
            I_nmda_slow[i]=I_nmda_slow[i]+s.i_nmda_mu*(s.i_nmda_s-I_nmda_slow[i])
            I_nmda_fast[i]=I_nmda_fast[i]+s.i_nmda_mu*(s.i_nmda_f-I_nmda_fast[i])

        for o from 0<=o<output_num:
        
    
            B= 1.0/(1.0 + (exp(s.mg1 * v_total[c]) / s.mg2))
	
            # this h is negative
            h=B*(v_total[c]-self.v_reversal)

            # this g_nmda is negative 
            I_nmda[c]= g_nmda[c]*(I_nmda_fast[c] + I_nmda_slow[c]) * h

            Ca[c]=Ca[c]+(I_nmda[c] - Ca[c]/s.tau_ca)*s.dt
    
            # p1=0.1; p2=p1*1e-4; p3=3; p4=1; % new paper % p2 ~=p1/1e-4 */
            # p1=0.5; p2=1e-4; p3=3; p4=0.05; % old paper */
    
            #  tau=p1[i]/(p2[i]+pow(Ca[i],p3[i]))+p4[i];   tau in seconds 
	    # eta=1.0/(tau*1000.0);   % eta in 1/ms */
    
	        
            #  eta=p1[i]*Ca[i]/1000; 
	
            if s.learning_rule==1: # difference of sigmoids, omega 
                omega=s.omega_offset+(sig(Ca[c]-s.alpha2,s.beta2)-
                                          0.5*sig(Ca[c]-s.alpha1,s.beta1))
                eta=s.eta_gamma0*Ca[c]
            elif s.learning_rule==2: # quadratic
                omega=s.omega_offset+ Ca[c]*(Ca[c]-s.theta_o)
                if omega>s.omega_max:
                    omega=s.omega_max
                if omega<0.0:
                    eta=s.eta_down
                else:
                    eta=s.eta_up
            elif s.learning_rule==3 or s.learning_rule==4: 
                # harel's hill function or  # biocyb omega and eta
    
                ##       p1: global magnitude
                ##       p4: basal calcium 
                ##       p3: exponent 
                ##       p2: half max 
                ##       want saturation, otherwise too much LTP 
                ##      
                ##       if p4 is 0, and Ca = p2 then we have half max 
                
                  #omega=P.omega_offset+sig(Ca(c)-P.alpha2,P.beta2)-
                      #0.25*sig(Ca(c)-P.alpha1,P.beta1);
    
                # shouldn't it be with .5?
                omega=s.omega_offset+(sig(Ca[c]-s.alpha2,s.beta2)-
                                          0.5*sig(Ca[c]-s.alpha1,s.beta1))
                eta=s.eta_gamma0*Ca[c]
    
                if s.learning_rule==3: # harels 
                    eta=s.eta_p1*pow(Ca[c]+s.eta_p4,s.eta_p3)/(pow(Ca[c]+s.eta_p4,s.eta_p3)+pow(s.eta_p2,s.eta_p3))
                else:                     # biocyb 
                    eta=(s.eta_p2+pow(Ca[c],s.eta_p3))/(s.eta_p1+s.eta_p4*(s.eta_p2+pow(Ca[c],s.eta_p3)))
            else:
                raise ValueError, "Unknown Learning Rule"
            
            weights[c]=weights[c]+eta*(omega-s.lambda_decay*weights[c])
            
            c=c+1
            
            
    # decay the backspike, NMDA currents, etc.
            
    for i from 0<=i<input_num:    # erase v_total
        I_nmda_slow[i]=I_nmda_slow[i]-I_nmda_slow[i]/s.tau_nmda_s
        I_nmda_fast[i]=I_nmda_fast[i]-I_nmda_fast[i]/s.tau_nmda_f
            
    for o from 0<=o<output_num:
        v_backspike_slow[o]=v_backspike_slow[o]-v_backspike_slow[o]/s.tau_backspike_slow
        v_backspike_fast[o]=v_backspike_fast[o]-v_backspike_fast[o]/s.tau_backspike_fast
       
    
                
cdef STDP_update(Connection_Group_struct *s,object self,double t):

    
    cdef int output_num,input_num
    cdef double *M,*P,*weights
    cdef char *outspike,*inspike
    
    cdef double tau_minus,tau_plus,a_minus,a_plus,sign,g_max,dt
    cdef int count,nearest_neighbor
    
    tau_minus=s.tau_minus
    tau_plus=s.tau_plus
    a_minus=s.a_minus
    a_plus=s.a_plus
    sign=s.sign
    g_max=s.g_max
    dt=s.dt
    nearest_neighbor=s.nearest_neighbor
    
    M=s.M
    P=s.P

    weights=s.weights
    outspike=s.outspike
    inspike=s.inspike
    output_num=s.output_num
    input_num=s.input_num
    
    
    
    cdef int o,i

    for o from 0<=o<output_num:
        M[o]=M[o]-M[o]*dt/tau_minus
        if outspike[o]:
            if nearest_neighbor:
                M[o]=-a_minus*sign
            else:
                M[o]=M[o]-a_minus*sign
            
    for i from 0<=i<input_num:
        P[i]=P[i]-P[i]*dt/tau_plus
        if inspike[i]:
            if nearest_neighbor:
                P[i]=a_plus*sign
            else:
                P[i]=P[i]+a_plus*sign
            
    count=0    
    for i from 0<=i<input_num:
        for o from 0<=o<output_num:

            weights[count]=weights[count]+(M[o]*inspike[i]*g_max+
              P[i]*outspike[o]*g_max)*sign

            count=count+1
            
cdef Gerstner06_update(Connection_Group_struct *s,object self,double t):
    
    cdef int output_num,input_num
    cdef double *r1,*r2,*o1,*o2,*P,*weights
    cdef char *outspike,*inspike
    
    cdef double tau_minus,tau_plus,tau_x,tau_y
    cdef double A2_minus,A2_plus,A3_minus,A3_plus
    cdef double sign,dt
    
    tau_minus=s.tau_minus
    tau_plus=s.tau_plus
    tau_x=s.tau_x
    tau_y=s.tau_y
    A2_minus=s.A2_minus
    A3_minus=s.A3_minus
    A2_plus=s.A2_plus
    A3_plus=s.A3_plus
    
    sign=s.sign
    dt=s.dt
    
    r1=s.r1
    r2=s.r2
    o1=s.o1
    o2=s.o2
    

    weights=s.weights
    outspike=s.outspike
    inspike=s.inspike
    output_num=s.output_num
    input_num=s.input_num
    
    
    
    cdef int o,i,c

    # update r1 and o1
    for o from 0<=o<output_num:
        o1[o]=o1[o]-o1[o]*dt/tau_minus
        if outspike[o]:
            if s.nearest_neighbor:
                o1[o]=1.0
            else:
                o1[o]=o1[o]+1.0
            
    for i from 0<=i<input_num:
        r1[i]=r1[i]-r1[i]*dt/tau_plus
        if inspike[i]:
            if s.nearest_neighbor:
                r1[i]=1.0
            else:
                r1[i]=r1[i]+1.0
            
                
                
    # update the weights
    c=0    
    for i from 0<=i<input_num:
        for o from 0<=o<output_num:
            
            
            weights[c]=weights[c] + (
                - inspike[i]*o1[o]*(A2_minus+A3_minus*r2[i])+
                outspike[o]*r1[i]*(A2_plus+A3_plus*o2[o]))
                
            c=c+1
            
            
            
            
    # update r2 and o2 afterwards
    for o from 0<=o<output_num:
        o2[o]=o2[o]-o2[o]*dt/tau_y
        if outspike[o]:
            if s.nearest_neighbor:
                o2[o]=1.0
            else:
                o2[o]=o2[o]+1.0
            
    for i from 0<=i<input_num:
        r2[i]=r2[i]-r2[i]*dt/tau_x
        if inspike[i]:
            if s.nearest_neighbor:
                r2[i]=1.0
            else:
                r2[i]=r2[i]+1.0

            
            
                
cdef cupdate(Connection_Group_struct *s,object self,double t):
    
    if s.type==-1:
        self.update(t)
        return
    
    if s.type==0: # Constant_Connection
        pass
    elif s.type==1: 
        Spiking_Rate_update(s,self,t)
    elif s.type==2: 
        Calcium_update(s,self,t)
    elif s.type==3: 
        STDP_update(s,self,t)
    elif s.type==4: 
        Gerstner06_update(s,self,t)
    else:
        raise TypeError,"Unknown type: %d" % (s.type)

