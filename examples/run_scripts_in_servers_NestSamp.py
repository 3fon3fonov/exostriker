import dill, os, sys


_orig_cwd = os.getcwd()
import exostriker.lib.RV_mod as rv
# Restore working directory
os.chdir(_orig_cwd)

#sys.path.append('../exostriker/lib/') #RV_mod directory must be in your path
#import RV_mod as rv
import time


file_name = "Eta_Ceti_py3.ses"

file_ = open(file_name,"rb")
fit = dill.load(file_)
file_.close() 

 
 
# Options for the NS sampler (see gui.py how these are passed from the GUI) 
#fit.ns_maxiter[1] = 6000000
#fit.ns_maxcall[1] = 6000000

# When options are done, either on the GUI side or here, just run the NS and 
# colect the the outout .ses file 

fit = rv.run_nestsamp(fit)


file_name = "Eta_Ceti_py3_NS_out.ses"

time.sleep(60) # to avoud some memory problems !


### save the session. Then you shuld be able to open it from the Exo-Striker
file_ = open(file_name,"wb")
dill.dump(fit,file_)
file_.close()


print("writing done!")
time.sleep(10)  # again, give it some time to avoud memory problems !

