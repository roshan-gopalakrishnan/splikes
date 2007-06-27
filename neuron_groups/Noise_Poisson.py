from Neuron_Group import *

class Noise_Poisson(Neuron_Group):
    
    def __init__(self,C=None,M=None,offset=0,
                time_between_rates=1000,type='gaussian',save_rates=False):
        
        if ( (not C is None) and (not M is None) ):
            raise ValueError,"Cannot specify both C and M"
        
        elif not C is None:
            M=self.get_mixing_matrix(C)
            
            
        if not M.shape[0]==M.shape[1]:
            raise ValueError,"Mixing matrix size isn't square: %dx%d" % (M.shape[0],M.shape[1],qty)
                     

        qty=M.shape[0]
        
        super(Noise_Poisson,self).__init__(qty)

            
        self.rate=0*numpy.ones((1,self.quantity),numpy.float)
        self.C=C
        self.M=numpy.matrix(M)
        self.time_between_rates=time_between_rates
        self.offset=numpy.array(offset)
        
        self.save_rates=save_rates
        
        self.type=type
        

        self.time_to_next_rate=0
        self.saved_rates=[]
        
    def _reset_(self):
        super(Noise_Poisson,self)._reset_()
        
        self.time_to_next_rate=0
        self.rate=0*numpy.ones((1,self.quantity),numpy.float)
        self.rate=self.rate.ravel()
        self.saved_rates=[]
        

    def get_mixing_matrix(self,C):
        from numpy import diag
        from numpy.linalg import eig,cholesky
        from warnings import warn
        
        D,V=eig(C)
        
        if any(D<0):  # negative eigenvalues
            warn("""The correlation function given has negative eigenvalues, which
                means that it is not a proper correlation function.  An
                approximate version will be used to create the environment
                by using only those vectors with positive eigenvalue.""")
        

            D[D<0]=1e-3*max(D)
            D=diag(D)
            C=V*D*V.T
            D=diag(D)
            
    
                
        M=cholesky(C)
        
        return M
                
        
        
    def update(self,t):
        
        self.swap_spikes()

        if t>=self.time_to_next_rate:
            self.time_to_next_rate+=self.time_between_rates
            
            if self.type=='gaussian':
	        self.rate=self.M.T*numpy.random.randn(self.quantity,1)
            elif self.type=='uniform':
	        self.rate=self.M.T*(numpy.random.rand(self.quantity,1)-0.5)*2*numpy.sqrt(3)
            elif self.type=='laplace':
	        self.rate=self.M.T*numpy.random.laplace(size=(self.quantity,1))/numpy.sqrt(2)
            else:
                raise TypeError,"Illegal noise type %s" % self.type
                
                
            self.rate=self.rate.ravel()+self.offset
            
            if self.save_rates:
                self.saved_rates.append(numpy.array(self.rate))
            
            #print self.rate
            
    
        self.spike[:]=numpy.random.rand(self.quantity)<(self.rate/1000.0/self.dt)
            
