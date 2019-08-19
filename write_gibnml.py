MODELVARS=[]
MODELDIR="/gpfs/scratch/userexternal/gbolzon0/GIB_allvars/wrkdir/MODEL/"
F=file(MODELDIR + "namelist.passivetrc")
for line in F:
    if line.find("ctrcnm") != -1:
        apex_1 = line.find("\"")
        MODELVARS.append(line[apex_1+1:apex_1+4])

print "&VARS_DIMENSION"
print "    n_vars = %d " % len(MODELVARS)
print "/"

print "\n&CORE\n"

for i, var in enumerate(MODELVARS):
    s = "    vars(%d) = \"%s\"" %(i+1, var)
    print s
print ""

for i, var in enumerate(MODELVARS):
    s ="    var_names_idx(%d) = %d" %(i+1,i+1)
    print s
print ""

print "    alpha = 4.0d0"
print "    reduction_value_t = 1.0d-6"
print "    length = -7.5d0"

print "\n/ \n &NUDGING_VARS_DIMENSION\n"
print "    n_vars = %d " % len(MODELVARS)
print "/\n\n&NUDGING_CORE\n"
print "data_file = \"bounmask.nc\" \n"

for i, var in enumerate(MODELVARS):
    s = "    vars(%d) = \"%s\"" %(i+1, var)
    print s
print ""

for i, var in enumerate(MODELVARS):
    s ="    var_names_idx(%d) = %d" %(i+1,i+1)
    print s
print ""

for i, var in enumerate(MODELVARS):
    s ="    rst_corr(%d) = 1.0d0" %(i+1)
    print s
print "/"
