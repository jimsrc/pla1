#!/usr/bin/env ipython
#comm.Barrier() # I need everyone to stop touching the files
# let's unify all files into one
#if rank==0:
fout = h5(fname_out, 'w')
for r in range(wsize):
    fnm_inp = fname_out+'_%02d'%r
    finp = h5(fnm_inp, 'r')
    cont = finp.keys()       # list of groups
    for c in cont:           # iterate over each group
        finp.copy(c, fout)
    finp.close()
    os.system('rm {fname}'.format(fname=fnm_inp))

print " ----> We generated: "+fout.filename
fout.close()
# clean backup
os.system('rm {fname}_'.format(fname=fname_out))
print " [r:%d] I'm finished!" % rank
#EOF
