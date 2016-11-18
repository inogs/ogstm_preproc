
def dom_dec(ji, jj, jpiglo, jpjglo, nproc_i,nproc_j):

    jpi  = jpiglo/nproc_i
    jpj  = jpjglo/nproc_j
    i_start=ji*jpi
    j_start=jj*jpj
    i_end = jpi*(ji+1)
    j_end = jpj*(jj+1)
    # last one takes rest of division
    if (ji==(nproc_i-1)): 
        i_end = jpiglo
        jpi = i_end - i_start
    if (jj==(nproc_j-1)): 
        j_end = jpjglo   
        jpj = j_end - j_start
    return i_start,i_end, j_start, j_end
#             print ji,jj
#             print i_start, i_end, i_end-i_start, jpi
#             print j_start, j_end, j_end-j_start, jpj

if __name__== "__main__":
    jpiglo=3308
    jpjglo=1580
    
    nproc_i=2
    nproc_j=2
    import numpy as np
    mglo=np.zeros((jpjglo,jpiglo),np.bool)
    
    for jj in range(nproc_j):
        for ji in range(nproc_i):
            i_start,i_end, j_start, j_end = dom_dec(ji, jj, jpiglo, jpjglo, nproc_i, nproc_j)
            Mlog = mglo[j_start:j_end,i_start:i_end]
            outfile = "clim_%d_%d" %(jj,ji)
            print outfile
        
