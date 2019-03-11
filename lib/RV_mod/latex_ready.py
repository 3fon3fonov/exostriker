 #!/usr/bin/python
__author__ = 'Trifon Trifonov'


#import sys 
#sys.path.insert(0, '../lib')
import numpy as np
from functions import *       


# class for signal to run gls and plotting functions on (it can be both the entire signal and O-C, depends what we pass as rvs in the __ini__ function)
# This class must completely refurbished!
class Latex_out(object): 

    def __init__(self,jd,rvs,rv_error, sig_for_gls=np.array([0.1,0.01,0.001])):
        self.jd=np.array(list(map(float,jd))) # for some reason sometimes we get strings instead of floats, so...
        self.rvs=np.array(list(map(float,rvs)))
        self.rv_error=np.array(list(map(float,rv_error)))
        self.sig=sig_for_gls # list of significances for gls, so the algorithm will do the computation for all three of them and include relevant levels in the z array outputed by lomb_scargle
        self.sig=np.array(self.sig)
        # It is convenient to store significance levels in decreasing order, so let's sort them
       
        sort = convert_array_to_int(np.array(sorted(range(len(self.sig)), key=lambda k: self.sig[k], reverse=True))) # find the permutation in which they are sorted in decreasing order 
        # sorting based on the permutation sort calculated above
        self.sig = self.sig[sort]             
        
  
    def pl_param_table(self, ind_sig=2,sig_for_gls=np.array([0.1,0.01,0.001])): # ind_sig will be the index of the significance level which we assume to be sufficient to declare a planet candidate
 
        self.sig = self.sig[sort] 
        
        
text == '''       
\begin{table}[ht]
% \begin{adjustwidth}{-4.0cm}{} 
% \resizebox{0.69\textheight}{!}
% {\begin{minipage}{1.1\textwidth}

\centering   
\caption{{}}   
\label{table:mod_stable}      

\begin{tabular}{lrrr}     % 2 columns 


\hline\hline  \noalign{\vskip 0.7mm}      

% 
% \multicolumn{4}{c}{ Updated two-planet coplanar edge-on  }}           \\
% \multicolumn{4}{c}{{\bf configuration based upon stability}}           \\


\hline \noalign{\vskip 0.7mm}  

Parameter &\hspace{0.0 mm} &  GJ\,1148 b  &  GJ\,1148~c  \\
\hline\noalign{\vskip 0.5mm}

Semi-amplitude $K$  [m\,s$^{-1}$]                        &  & 34.1$_{-2.9}^{+2.3}$          &  25.3$_{-3.3}^{+2.1}$        \\
Period $P$ [days]   		        	          &  & 517.8$_{-3.9}^{+8.9}$        &  946.6$_{-20.9}^{+20.7}$      \\  
Eccentricity $e$                                       &  & ~~0.011$_{-0.011}^{+0.078}$   &  ~~0.028$_{-0.012}^{+0.065}$ \\
Arg. of periastron $\omega$ [deg]                            &  & ~~201.6$_{-98.5}^{+97.8}$    &  ~~136.8$_{-62.7}^{+143.0}$       \\   
Mean anomaly $M_0$ [deg]                               & &  ~~79.2$_{-28.1}^{+198.6}$    &  180.0$_{-108.9}^{+98.9}$    \\  \noalign{\vskip 0.9mm} 
Mean longditude $\lambda$ [deg]                           & &  ~~298.8$_{-14.6}^{+70.0}$    &  262.8$_{-53.3}^{+37.8}$    \\ 
Semi-major axis $a$ [AU]                                  &  & 1.566$_{-0.007}^{+0.016}$     &  2.342$_{-0.035}^{+0.034}$                     \\  
Dyn. mass $m$ [$M_{\mathrm{Jup}}$]                  & &  1.996$_{-0.100}^{+0.220}$     &  1.864$_{-0.227}^{+0.177}$                    \\ \noalign{\vskip 0.9mm}

Inclination $i$ [deg]                                 & & 90.0 (fixed)       & 90.0 (fixed)              \\ 
Node $\Delta \Omega$ [deg]                     & &  \multicolumn{2}{c}{0.0 (fixed) }                    \\
%$\Delta i$ [deg]                          & 180.0                &                            \\ \noalign{\vskip 0.9mm} 

$\gamma_{\rm HIRES}$~[m\,s$^{-1}$]        &  & \multicolumn{2}{c}{--14.42$_{-2.22}^{+1.70}$}  \\ \noalign{\vskip 0.9mm} 
$\gamma_{\rm CARM.}$~[m\,s$^{-1}$]        &  & \multicolumn{2}{c}{--14.42$_{-2.22}^{+1.70}$}  \\  


Jitter$_{\rm HIRES}$ [m\,s$^{-1}$]        &  & \multicolumn{2}{c}{9.36$_{-0.57}^{+2.03}$}   \\ \noalign{\vskip 0.9mm}
Jitter$_{\rm CARM.}$ [m\,s$^{-1}$]        &  & \multicolumn{2}{c}{9.36$_{-0.57}^{+2.03}$}   \\   
 

%$r.m.s. $ [m\,s$^{-1}$]                   & 19.46               &                             \\ 
%$-\ln\mathcal{L}$                       &   &  \multicolumn{2}{c}{154.51$_{-3.74}^{+1.40}$} \\                                     
%$\chi^2$                                 &   & 41.6830              &                            \\  
%$\chi_{\nu}^2$                           &   & 1.3894              &                            \\
\noalign{\vskip 0.5mm}

\hline\hline\noalign{\vskip 1.2mm}  

% 
% \multicolumn{3}{c}{Mutually inclined}            \\
% 
% % &  \multicolumn{1}{c}{ ($i$~=~90$^\circ$, $\Delta\Omega$~=~180$^\circ$) } \\ 
% 
% \hline \noalign{\vskip 0.7mm}  
 


\end{tabular}  

%\tablefoot{\small same comments as in Table~\ref{table:orb_par_stable} }

\end{table}
'''  

