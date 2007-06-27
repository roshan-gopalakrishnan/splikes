from pylab import *
ion()

def plot_spikes(n_list):
    
    markers=['o','+','x','D','^','s','v','h','p']
    colors=['b','r','g','m','k']
    
    color_count=0
    marker_count=0
    
    qty=0
    for n in n_list:
        s=n.saved_spikes
        
        style=colors[color_count]+markers[marker_count]
        
        plot(s['t'],[a+qty for a in s['n']],style)
        gca().hold(True)
        qty=qty+n.quantity
        
    
        color_count+=1
        if color_count>=len(colors):
            color_count=0
            marker_count+=1
            if marker_count>=len(markers):
                marker_count=0


    gca().set_ylim([-0.5,qty-1+0.5])
    draw()
    
