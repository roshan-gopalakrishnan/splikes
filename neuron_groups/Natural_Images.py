from Neuron_Group import *

default_pics_path='/home/bblais/python/work/plasticity/plasticity/pics/'


def scale_shift_filter(image,scale=1.0,shift=0.0):

    for im in image:

        im[:]=scale*im+shift



def min_max_filter(image,min=0.0,max=1.0):

    for im in image:

        im[:]= (im-im.min())/(im.max()-im.min())*(max-min)+min


class Natural_Images(Neuron_Group):
    
    def __init__(self,qty,time_between_rates=1000,
                        pics=['new_images12_dog.pics'],
                        filter=[min_max_filter],
                        filter_params=[{'min':0.0,'max':50.0}],
                        rate_range=None):
        
        if isinstance(qty,list):
            total_qty=sum(qty)
            num_channels=len(qty)
        else:
            total_qty=qty
            qty=[qty]
            num_channels=1
        
            
        if rate_range:
            for f in filter_params:
                f['min']=rate_range[0]
                f['max']=rate_range[1]
            
            
        if len(pics)!=num_channels:
            pics=pics*num_channels
                            
        if len(filter)!=num_channels:
            filter=filter*num_channels
            
        if len(filter_params)!=num_channels:
            filter_params=filter_params*num_channels
            
            
        super(Natural_Images,self).__init__(total_qty)
        
        import sys
        sys.path.append('/home/bblais/python/work/plasticity')
        
        import plasticity

        
        self.plasticity=plasticity
        
        self.time_between_rates=time_between_rates
        self.time_to_next_rate=0
        self.num_channels=num_channels
        
        params=plasticity.utils.default_params()
        self.params={}
        
        keys=['neuron_offsets','num_neurons',
              'pattern_input','noise_input']
        
        for k in keys:
            self.params[k]=params[k]
            
        self.params['pattern_input']=[]
        self.params['noise_input']=[]
        for c in range(num_channels):
            self.params['pattern_input'].append(plasticity.utils.get_pattern())
            self.params['noise_input'].append(plasticity.utils.get_noise())

        for c in range(num_channels):
            self.params['pattern_input'][c]['num_inputs']=qty[c]
            self.params['pattern_input'][c]['filter']=filter[c]
            self.params['pattern_input'][c]['filter_params']=filter_params[c]

            if pics[c]:
                self.params['pattern_input'][c]['filename']=default_pics_path+pics[c]
            else:
                self.params['pattern_input'][c]['type']=0

        num_neurons=1
        num_inputs=total_qty
        
        self.X=numpy.zeros((num_neurons,num_inputs),numpy.float64)
        
        
        plasticity.train.init_params(self.params)
        
        
            
        
    def _reset_(self):
        
        # reset for time t=0
        
        super(Natural_Images,self)._reset_()
        self.time_to_next_rate=0
        self.plasticity.train.init_params(self.params)
        
        
    def update(self,t):
        
        
        self.swap_spikes()

        
        if t>=self.time_to_next_rate:
            self.plasticity.train.get_input_vector(self.X,self.params)
            
            self.time_to_next_rate+=self.time_between_rates

            
            
        self.spike[:]=numpy.random.rand(self.quantity)<(self.X[:]/1000.0/self.dt)

