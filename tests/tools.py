    
    
def print_rates(n1,total_time):
    ts=n1.saved_spikes['t']
    ns=n1.saved_spikes['n']
    
    
    try: # assume list
        total_time_val=total_time[1]-total_time[0]
        time_min=total_time[0]
        time_max=total_time[1]
    except TypeError:
        total_time_val=total_time
        time_min=0
        time_max=total_time
    
    
    
    for nn in range(n1.quantity):
        tt=[t for n,t in zip(ns,ts) if n==nn and (time_min<=t<=time_max)]
        print "  %d: %.1f Hz %d spikes" % (nn,
                                           len(tt)/float(total_time_val)*1000,
                                           len(tt))
    
