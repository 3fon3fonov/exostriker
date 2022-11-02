import dill, os, sys
sys.path.append('../exostriker/lib/') #RV_mod directory must be in your path
import RV_mod as rv
import time


file_name = "Your_saved_session_with_NS_setup.ses"

file_ = open(file_name,"rb")
fit = dill.load(file_)
file_.close() 

fit.cwd = '../exostriker'


#fit.ns_maxiter[1] = 6000000
#fit.ns_maxcall[1] = 6000000

fit = rv.run_nestsamp(fit)


file_name = "output_session_with_NS_samples.ses"

time.sleep(60) # to avoud some memory problems !


### save the session. Then you shuld be able to open it from the Exo-Striker
file_ = open(file_name,"wb")
dill.dump(fit,file_)
file_.close()


print("writing done!")
time.sleep(10)  # again, give it some time to avoud memory problems !

