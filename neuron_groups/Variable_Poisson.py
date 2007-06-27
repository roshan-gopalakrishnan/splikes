from Neuron_Group import *

class Variable_Poisson(Neuron_Group):
    
    def __init__(self,qty,rate=None,time_between_rates=1000,sequential=True):
        
        super(Variable_Poisson,self).__init__(qty)

        if rate is None:
            self.rate=50*numpy.ones((1,self.quantity),numpy.float)
        elif isinstance(rate,basestring):  # filename
            d=zpickle.load(rate)
            self.rate=d['var']
        else:
            self.rate=rate.reshape((numpy.prod(rate.shape)/qty,qty))
            
        self.time_between_rates=time_between_rates
        self.sequential=sequential
        
        self.time_to_next_rate=0
        self.which_rate=-1
        
    def _reset_(self):
        super(Variable_Poisson,self)._reset_()
        
        self.time_to_next_rate=0
        self.which_rate=-1
        
        if len(self.rate.shape)==1:
            self.rate.shape=(1,len(self.rate))
        
        
    def update(self,t):
        
        self.swap_spikes()

        number_of_rates=self.rate.shape[0]
        
        if t>=self.time_to_next_rate:
            self.time_to_next_rate+=self.time_between_rates
            if self.sequential:
                self.which_rate+=1
                if self.which_rate>=number_of_rates:
                    self.which_rate=0
            else:
                self.which_rate=numpy.random.randint(number_of_rates)

        self.spike[:]=numpy.random.rand(self.quantity)<(self.rate[self.which_rate,:]/1000.0/self.dt)
    
    # replacing the update method
#    update=wrapAsMethod(update_methods.Variable_Poisson_update)
     
