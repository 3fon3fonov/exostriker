#!/usr/bin/env python3
"""
@author: Trifon Trifonov
"""

import sys,  os
sys.path.insert(0, './lib')
sys.path.append('./exostriker/lib/RV_mod/') #RV_mod directory must be in your path

#lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')
#sys.path.insert(0,lib_path)
 
#os.chdir(os.path.dirname(os.path.abspath(__file__)))

import RV_mod as rv
import dill
 

arguments = len(sys.argv) - 1

fit=rv.signal_fit(name='mcmc_ses')


if arguments==3 and sys.argv[1] == '-ses' and os.path.exists(sys.argv[2]):
    try:
        file_pi = open(sys.argv[2], 'rb')
        fit_ses = dill.load(file_pi)
        file_pi.close() 
        fit = rv.check_for_missing_instances(fit,fit_ses)

 
    except (ImportError, KeyError, AttributeError) as e:
        print("%s cannot be recognaized"%sys.argv[2])
 
    target_name = sys.argv[3] 
    #fit_ses.cwd = '/home/trifonov/git/exostriker_TIC149601126/exostriker'

    fit = rv.run_mcmc(fit)
      
    file_ses = open("%s_out.ses"%target_name, 'wb')
    dill.dump(fit, file_ses)
    file_ses.close()    
    
    


     
