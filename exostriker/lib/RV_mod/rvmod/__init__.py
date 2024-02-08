#!/usr/bin/env python3

#####################################################
# RVMOD Class                                       #
# Matheus J. Castro & Trifon Trofonov               #
# v0.8                                              #
# Last Modification: 09/21/2023                     #
# Contact: matheusdejesuscastro@gmail.com           #
#####################################################

from pathlib import Path
import numpy as np
import sys
import os

from .Warning_log import Warning_log


class GP_parameters(object): # class for GP process parameters which allows for different kernels with different number of parameters

    def __init__(self,npar,parameters,kernel_id=0):
        gpparameterswarnings=Warning_log([],'generating GP_parameters object')
        self.gp_par=parameters
        if not (npar==len(parameters)):
            npar=len(parameters)
            gpparameterswarnings.update_warning_list('Different number of parameters declared than the number of parameters actually provided! Redefined.')
        self.npar=npar
        self.kernel_id=kernel_id
        #self.rot_kernel=rot_kernels.choose_kernel(kernel_id,parameters)
        gpparameterswarnings.print_warning_log()


# RvModelWrap Class
# Wrap to create rv_model variable with the right parameters
class RvModelWrap:
    def __init__(self, jd, rvs, rv_err, o_c):
        self.jd = jd
        self.rvs = rvs
        self.rv_err = rv_err
        self.o_c = o_c


# ParamsWrap Class
# Wrap to create stat variable with the right parameters
class ParamsWrap:
    def __init__(self, params, param_errors, dof=0):
        self.params = params
        self.param_errors = param_errors
        self.dof = dof


# PlanetParamsWrap Class
# Wrap to create planet_params variable with the right parameters
class PlanetParamsWrap:
    def __init__(self, offsets, jitters, planet_params, linear_trend, stellar_mass,
                 GP_params=None, GP_kernel_id=0):
        if GP_params is None:
            GP_params = [10.0] * 4

        self.offsets = offsets
        self.jitters = jitters
        self.planet_params = planet_params
        self.linear_trend = linear_trend
        self.stellar_mass = stellar_mass

        self.GP_params = GP_parameters(len(GP_params), GP_params,
                                       kernel_id=GP_kernel_id)
        # we always want to have this attribute, but we only use it if we call GP, and then we update it anyway

    def update_offset(self, dataset, offset):  # change offset of one dataset
        self.offsets[dataset] = offset
        return

    def update_offsets(self, offsets):  # change all offsets
        self.offsets = offsets
        return

    def update_jitter(self, dataset, jitter):  # change jitter of one dataset
        self.jitters[dataset] = jitter
        return

    def update_jitters(self, jitters):  # change all jitters
        self.jitters = jitters
        return

    def update_K(self, planet, newK):  # update K for a given planet
        self.planet_params[7 * planet] = newK
        return

    def update_P(self, planet, newP):  # update P for a given planet
        self.planet_params[7 * planet + 1] = newP
        return

    def update_e(self, planet, newe):  # update e for a given planet
        self.planet_params[7 * planet + 2] = newe
        return

    def update_w(self, planet, neww):  # update w for a given planet
        self.planet_params[7 * planet + 3] = neww
        return

    def update_M0(self, planet, newM0):  # update M0 for a given planet
        self.planet_params[7 * planet + 4] = newM0
        return

    def update_inclination(self, planet, newi):  # update inclination info for one planet
        self.planet_params[7 * planet + 5] = newi
        return

    def update_lineofnodes(self, planet, newcap):  # update lineofnodes info for one planet
        self.planet_params[7 * planet + 6] = newcap
        return

    def update_planet_params_one_planet(self, planet, K, P, e, w, M0, i,
                                        cap):  # update all planet_params for one planet
        self.update_K(planet, K)
        self.update_P(planet, P)
        self.update_e(planet, e)
        self.update_w(planet, w)
        self.update_M0(planet, M0)
        self.update_inclination(planet, i)
        self.update_lineofnodes(planet, cap)
        return

    def update_planet_params(self, planet_params):  # update all planet_params in one go
        self.planet_params = planet_params
        return

    def update_linear_trend(self, linear_trend):  # update linear trend
        self.linear_trend = linear_trend
        return

    def update_GP_param_value(self, i, newpar):
        self.GP_params.gp_par[i] = newpar
        return

    def update_GP_params(self, newparams, kernel_id=-1):
        # redefine entire GP_params object, if kernel_id=-1 then we do not wish to change kernel type
        if (kernel_id == -1):
            kernel_id = self.GP_params.kernel_id
        # self.GP_params=GP_parameters(newparams,newparams,kernel_id=kernel_id)
        self.GP_params = GP_parameters(len(newparams), newparams, kernel_id=kernel_id)
        return

    def update_stellar_mass(self, stellar_mass):
        self.stellar_mass = stellar_mass
        return

    def getinit_input(self, masses, semimajor, fileinput=False, filename='geninit_j_input'):

        '''Prepares input for a fortran code which calculates Jacobi coordinates based on orbital parameters'''

        if not (fileinput):  # if we want to save input in a file we don't want this line in the input string
            ppp = './geninit_j3_in_days << EOF\n'
        else:
            ppp = ''

        npl = len(self.planet_params) / 7

        ppp += '1 \n%s \n%s \n1.d0 \npl.in\n' % (str(self.stellar_mass), str(npl))

        for i in range(npl):
            ppp += '%s \n%s %s %s %s %s %s\n' % (
            str(masses[i]), str(semimajor[i]), str(self.planet_params[7 * i + 2]),
            str(self.planet_params[7 * i + 5]), str(self.planet_params[7 * i + 3]),
            str(self.planet_params[7 * i + 6]), str(self.planet_params[7 * i + 4]))

        if not (fileinput):
            ppp += 'EOF'  # end of the command to run in the case of saving input directly
        else:  # here's what we do if we want to generate a file as well
            # first we save the ppp string in a file (by default 'Kep_input')
            file_getin = open(filename, 'wb')
            file_getin.write('%s' % ppp)
            file_getin.close()
            # then we overwrite ppp with the command to pass this file as input for the fortran code
            ppp = './geninit_j3_in_days < %s' % (program, filename)

        return ppp


# PlanetParamsWrap Class
# Wrap to create planet_params variable with the right parameters
class PlanetParamsErrorsWrap:
    def __init__(self, offsets, jitters, planet_params, linear_trend, stellar_mass,
                 GP_params_errors=None):
        if GP_params_errors is None:
            GP_params_errors = [[0.0, 0.0]] * 4

        self.offset_errors = offsets
        self.jitter_errors = jitters
        self.planet_params_errors = planet_params
        self.linear_trend_error = linear_trend
        self.stellar_mass_error = stellar_mass

        self.GP_params_errors = GP_params_errors

        '''In all functions below 'error' should in fact be a [-error,+error] 2 element array'''

    def update_offset_error(self, dataset, offset_error):  # change offset error of one dataset
        self.offset_errors[dataset] = offset_error
        return

    def update_offset_errors(self, offset_errors):  # change all offset errors
        self.offset_errors = offset_errors
        return

    def update_jitter_error(self, dataset, jitter_error):  # change jitter error for one dataset
        self.jitter_errors[dataset] = jitter_error
        return

    def update_jitter_errors(self, jitter_errors):  # change all jitters
        self.jitter_errors = jitter_errors
        return

    def update_Kerror(self, planet, newKerror):  # update K error for a given planet
        self.planet_params_errors[7 * planet] = newKerror
        return

    def update_Perror(self, planet, newPerror):  # update P error for a given planet
        self.planet_params_errors[7 * planet + 1] = newPerror
        return

    def update_eerror(self, planet, neweerror):  # update e error for a given planet
        self.planet_params_errors[7 * planet + 2] = neweerror
        return

    def update_werror(self, planet, newwerror):  # update w error for a given planet
        self.planet_params_errors[7 * planet + 3] = newwerror
        return

    def update_M0error(self, planet, newM0error):  # update M0 error for a given planet
        self.planet_params_errors[7 * planet + 4] = newM0error
        return

    def update_inclination_error(self, planet, newierror):  # update inclination error for one planet
        self.planet_params_errors[7 * planet + 5] = newierror
        return

    def update_lineofnodes_error(self, planet, newcaperror):  # update lineofnodes error for one planet
        self.planet_params_errors[7 * planet + 6] = newcaperror
        return

    def update_planet_param_errors_one_planet(self, planet, Kerror, Perror, eerror, werror, M0error, ierror,
                                              lineofnodeserror):  # update all planet_params_errors for one planet
        self.update_Kerror(planet, Kerror)
        self.update_Perror(planet, Perror)
        self.update_eerror(planet, eerror)
        self.update_werror(planet, werror)
        self.update_M0error(planet, M0error)
        self.update_inclination_error(planet, ierror)
        self.update_lineofnodes_error(planet, lineofnodeserror)
        return

    def update_planet_param_errors(self, planet_params_errors):  # update all planet_param_errors in one go
        self.planet_params_errors = planet_params_errors
        return

    def update_linear_trend_error(self, linear_trend_error):  # update linear trend error
        self.linear_trend_error = linear_trend_error
        return

    def update_GP_param_errors(self, i, newerror):
        self.GP_params_errors[i] = newerror
        return

    def update_GP_params_errors(self, GP_params_errors):
        self.GP_params_errors = GP_params_errors
        return

    def update_stellar_mass_error(self, stellar_mass_error):
        self.stellar_mass_error = stellar_mass_error
        return


# Rvfit Class
# This class uses the rvmod_for.f95 as a Python library compiled with F2PY
class Rvfit:
    def __init__(self):
        self.arguments = None
        self.res_dict = {}

        self.res, self.model_rvs, self.rvs, self.o_c, self.rv_err, \
            self.t0, self.RR, self.aR, self.m, self.tw, self.h, self.k, \
            self.lamb, self.V0, self.jitt, self.lin_trend_out, self.quad_trend_out, \
            self.Ndata, self.mfit, self.rms2, \
            self.reduced_chi22, self.jacobi_major_ax, self.fit, self.npl, \
            self.model_jd, self.model, self.rv_lintr, self.rv_lintr_err, \
            self.rv_quadtr, self.rv_quadtr_err, self.stat, self.to_sort, self.npl_pos, \
            self.stat_array_saved, self.stat, self.stellar_mass, self.ind_orig, self.ind_used = [None] * 38

        self.rms, self.wrms, self.chi2, self.reduced_chi2, self.dof, self.loglik, \
            self.P, self.K, self.e, self.w, self.M0, self.incl, self.cap0m, self.epoch,\
            self.ndset, self.ndata, self.coplar_inc = [0] * 17
        self.a = [0] * 20
        self.jd, self.idset = [[]]*2
        self.rv_model = RvModelWrap(self.jd, self.rvs, self.rv_err, self.o_c)

        self.omega_dot, self.mass = [list(np.zeros(9))] * 2
        self.omega_dot_err = np.zeros((9, 2))

        self.offsets, self.jitters, self.planet_params = np.zeros(20), np.zeros(20), np.zeros(70)
        self.offset_errors, self.jitter_errors, self.planet_params_errors = \
            np.zeros((20, 2)), np.zeros((20, 2)), np.zeros((70, 2))

        self.params = PlanetParamsWrap(self.offsets, self.jitters, self.planet_params,
                                       self.lin_trend_out, self.stellar_mass)
        self.params_errors = PlanetParamsErrorsWrap(self.offset_errors, self.jitter_errors, self.planet_params_errors,
                                                    self.lin_trend_out, self.stellar_mass)
        self.stat = ParamsWrap(self.params, self.params_errors, dof=self.dof)

    # Print the documentation for the class and all callable subroutines
    @staticmethod
    def __doc__():
        return rvmod_for.__doc__ + "\n" + rvmod_for.kepfit_amoeba.__doc__ + "\n" + rvmod_for.kepfit_lm.__doc__

    def __str__(self):
        return "Code Initialized. Call functions to run and get results."

    @staticmethod
    def __fix_index(data):
        data = np.array(data)
        for i, j in enumerate(np.unique(data.T[3])):
            for k in range(len(data)):
                if data[k, 3] == j:
                    data[k, 3] = i + 1
        return data

    @staticmethod
    def __sort_planets(data):
        data = np.array(data)
        to_sort = data[:, 1, 0].argsort()
        data = data[to_sort]
        return data, to_sort

    # Initialization of the arguments that should be passed to run the Fortran code
    def init_args(self, data_array, planets_param, epoch, hkl, dyn_eps=1e-10, dyn_dt=864000., amoeba_iter=0,
                  timeout=10, rv_model_npoints=5000, rv_model_max=50, rv_model_min=0, rv_gr_flag=0, stellar_mass=1.,
                  get_best_par=0, get_RV=0, get_fit_model=0,
                  ndset=None, ndata=None, nplanet=None,
                  rv_ofset=None, rv_jitt=None,
                  lin_trend_in=None, quad_trend_in=None,
                  dyn_planets=None, coplar_inc=0, npl_pos=None):

        self.ind_used = np.unique(np.array(data_array).T[3]).astype(int) if len(data_array) != 0 else [0]

        # Set default arguments for optional entries
        if ndset is None:
            ndset = len(self.ind_used) if len(data_array) != 0 else 0
        if ndata is None:
            ndata = len(data_array)
        if list(planets_param) == []:
            planets_param = [[[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                              [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                              [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]]
            nplanet = 0
        elif nplanet is None:
            nplanet = len(planets_param)
        if rv_ofset is None or ndset == 0:
            rv_ofset = [[0, 0]]
        if rv_jitt is None or ndset == 0:
            rv_jitt = [[0, 0]]
        if lin_trend_in is None:
            lin_trend_in = [0., 0]
        if quad_trend_in is None:
            quad_trend_in = [0., 0]
        if dyn_planets is None:
            if nplanet != 0:
                dyn_planets = np.ones(nplanet)
            else:
                dyn_planets = [1]

        # Sort the input array respectively from the JD collumn or set empty array
        if len(data_array) != 0:
            data_array = data_array[data_array[:, 0].argsort()]
            self.ind_orig = np.array(data_array).T[3]
        else:
            data_array = [[0, 0, 0, 0]]
            self.ind_orig = np.array([])

        data_array = self.__fix_index(data_array)

        if nplanet != 0:
            planets_param, self.to_sort = self.__sort_planets(planets_param)
        else:
            self.to_sort = []

        if npl_pos is not None:
            ocurr = 0
            for i in npl_pos:
                ocurr += 1 if i == 1 else 0
            if ocurr != len(self.to_sort):
                sys.exit("Error: positions of the planets does not match the number of planets.")

        self.npl = nplanet
        self.ndset = int(ndset)
        self.ndata = ndata
        self.stellar_mass = stellar_mass
        self.coplar_inc = coplar_inc
        self.npl_pos = npl_pos

        # Concatenate some of the arguments to be in the right form for Fortran code
        files_param = np.concatenate([rv_ofset, rv_jitt], axis=1)
        final_params = np.concatenate([lin_trend_in, quad_trend_in, [epoch], [hkl]])

        # Create the final multidimensional array accepted in Fortran code
        self.arguments = [dyn_eps, dyn_dt, amoeba_iter, timeout, rv_model_npoints,
                          rv_model_max, rv_model_min, rv_gr_flag, stellar_mass,
                          get_best_par, get_RV, get_fit_model,
                          ndset, ndata, data_array, files_param,
                          nplanet, planets_param, final_params,
                          dyn_planets, coplar_inc]
        return self.arguments

    # Internal subroutine. Used to auto update the arguments respectively from the last run
    def __update_args(self, res):
        data_array = self.arguments[14]
        planets_param = res[2]
        epoch = self.res_dict["epoch"]
        hkl = self.arguments[18][-1]
        dyn_eps = self.arguments[0]
        dyn_dt = self.arguments[1]
        amoeba_iter = self.arguments[2]
        timeout = self.arguments[3]
        rv_model_npoints = self.arguments[4]
        rv_model_max = self.arguments[5]
        rv_model_min = self.arguments[6]
        rv_gr_flag = self.arguments[7]
        stellar_mass = self.arguments[8]
        get_best_par = self.arguments[9]
        get_RV = self.arguments[10]
        get_fit_model = self.arguments[11]
        rv_ofset = np.concatenate([[self.res_dict["V0"][:, 0]], [self.arguments[15][:, 1]]], axis=0).T
        rv_jitt = np.concatenate([[self.res_dict["jitt"][:, 0]], [self.arguments[15][:, 3]]], axis=0).T
        lin_trend_in = self.res_dict["lin_trend_out"]
        quad_trend_in = self.res_dict["quad_trend_out"]
        dyn_planets = self.arguments[19]

        self.init_args(data_array, planets_param, epoch, hkl, dyn_eps=dyn_eps, dyn_dt=dyn_dt,
                       amoeba_iter=amoeba_iter, timeout=timeout, rv_model_npoints=rv_model_npoints,
                       rv_model_max=rv_model_max, rv_model_min=rv_model_min, rv_gr_flag=rv_gr_flag,
                       stellar_mass=stellar_mass, get_best_par=get_best_par, get_RV=get_RV,
                       get_fit_model=get_fit_model, rv_ofset=rv_ofset, rv_jitt=rv_jitt,
                       lin_trend_in=lin_trend_in, quad_trend_in=quad_trend_in,
                       dyn_planets=dyn_planets)

    # Compile the Fortran code
    @staticmethod
    def f_compile(path=None):
        if path is None:
            path = Path(__file__).parts[:-1]
            path = Path(path[0]).joinpath(*path[1:])

        # Change directory
        old_path = os.getcwd()
        os.chdir(path)

        # Get the actual Python version
        vers = str(sys.version_info.major) + "." + str(sys.version_info.minor)

        # Recursive Flag only on linux
        #recur = "" if "win" in sys.platform[0:3] else "-frecursive"
        recur = "-fmax-stack-var-size=2147483646" if "win" in sys.platform[0:3] else "-fmax-stack-var-size=2147483646"        
#        recur = "-frecursive" if "win" in sys.platform[0:3] else "-frecursive"                 

        # Compile it using the Numpy F2PY
        os.system("python{} -m numpy.f2py -c --opt=\"-O3 -std=legacy {}\" -m rvmod_for rvmod_for.f95".format(vers, recur))

        # If Windows, move the created DLL
        if "win" in sys.platform[0:3]:
            lib_path = path.joinpath("rvmod_for", ".libs")
            fls = os.listdir(lib_path)
            for fl in fls:
                os.rename(lib_path.joinpath(fl), path.joinpath(fl))
            os.rmdir(lib_path)
            os.rmdir(path.joinpath("rvmod_for"))

        # Return to the old directory
        os.chdir(old_path)

    # Check if there is the right compiled Fortran code for the Python version
    # If there is not, compile it
    # This function runs automatically in every initialization of the RVMOD module
    def check_compiled_version(self, path=None):
        # Rename the old executables, so Python will not try to import the wrong one
        def rename(fl_names):
            for exec_fl in fl_names:
                if "old_" not in exec_fl:
                    num = 1
                    while True:
                        if "old_{}_".format(num) + exec_fl in os.listdir():
                            num += 1
                        else:
                            os.rename(exec_fl, "old_{}_".format(num) + exec_fl)
                            break

        if path is None:
            path = Path(__file__).parts[:-1]
            path = Path(path[0]).joinpath(*path[1:])

        # Change directory
        old_path = os.getcwd()
        os.chdir(path)

        # Check for others executables in the path
        execs = []
        for names in os.listdir():
            if "rvmod_for" in names:
                if (".so" in names and "linux" in sys.platform) or (".so" in names and "darwin" in sys.platform) or (".pyd" in names and "win" in sys.platform[0:3]):
                    execs.append(names)
        # If there is no other executable, call the compile function, otherwise rename them
        if execs is not []:
            # Get python version and check for executables for this version
            vers = str(sys.version_info.major) + str(sys.version_info.minor)
            compile_flag = any(vers in exec_fl for exec_fl in execs)
            if not compile_flag:
                rename(execs)
                self.f_compile()
            else:
                execs_true = np.array([vers in exec_fl for exec_fl in execs])
                execs_ver = np.array(execs)[np.where(execs_true == True)[0]]
                if all("old_" in exec_fl for exec_fl in execs_ver):
                    ind, max_val = -1, 1
                    for i in range(len(execs_ver)):
                        if max_val < int(execs_ver[i][4]):
                            max_val = int(execs_ver[i][4])
                            ind = i
                    os.rename(execs_ver[ind], execs_ver[ind][6:])
                execs_ren = np.array(execs)[np.where(execs_true == False)[0]]
                rename(execs_ren)
        else:
            self.f_compile()

        # Return to the old directory
        os.chdir(old_path)

    # Run the amoeba code in Fortran
    # The options for mtype defines the type of run between Keplerian and N-body
    def run_amoeba(self, mtype, auto_update=False):
        if mtype == "kep":
            try:
                res = rvmod_for.kepfit_amoeba(*self.arguments)
            except Exception as e:
                print(repr(e))
                return 0
        elif mtype == "dyn":
            try:
                res = rvmod_for.dynfit_amoeba(*self.arguments)
            except Exception as e:
                print(repr(e))
                return 0
        else:
            exit("Model type not recognized!")

        self.__create_dict(res)
        self.__update_args(res) if auto_update else None
        return 1

    # Run the Levenberg-Marquardt code in Fortran
    # The options for mtype defines the type of run between Keplerian and N-body
    def run_lm(self, mtype, auto_update=False):
        if mtype == "kep":
            try:
                res = rvmod_for.kepfit_lm(*self.arguments)
            except Exception as e:
                print(repr(e))
                return 0
        elif mtype == "dyn":
            try:
                res = rvmod_for.dynfit_lm(*self.arguments)
            except Exception as e:
                print(repr(e))
                return 0
        else:
            exit("Model type not recognized!")

        self.__create_dict(res)
        self.__update_args(res) if auto_update else None
        return 1

    @staticmethod
    def __sort_array_by_index(array, index):
        if len(array) == 0:
            return array
        array_temp = np.zeros((max(index), 2))
        for i, j in enumerate(index-1):
            array_temp[j, :] = array[i]
        return array_temp

    def __unsort_planets_args(self, param):
        if len(param) == 0:
            if self.npl_pos is not None:
                return np.zeros(len(self.npl_pos))
            else:
                return np.array([])

        unsort = np.argsort(self.to_sort)
        param = np.array(param)[unsort]
        param_out = []
        k = 0
        if self.npl_pos is not None:
            for i in self.npl_pos:
                if i == 1:
                    param_out.append(param[k])
                    k += 1
                else:
                    param_out.append(0)
        else:
            return param.copy()

        try:
            return np.array(param_out)
        except ValueError:
            return param_out

    # Internal subroutine. Take the output of Fortran and format it as a dict
    def __create_dict(self, res):
        self.res_dict.clear()
        self.res_dict = {"jd": res[0].T[0], "model_rvs": res[0].T[5], "rvs": res[0].T[1],
                         "o_c": res[0].T[4], "rv_err": res[0].T[2], "idset": res[0].T[3]-1,
                         "loglik": res[1][0], "reduced_chi2": res[1][1],
                         "chi2": res[1][2], "rms": res[1][3],
                         "K": res[2][:, 0], "P": res[2][:, 1], "e": res[2][:, 2],
                         "w": res[2][:, 3], "M0": res[2][:, 4], "incl": res[2][:, 5],
                         "cap0m": res[2][:, 6], "w_dot": res[2][:, 7], "t0": res[2][:, 8],
                         "RR": res[2][:, 9], "aR": res[2][:, 10], "a_new": res[2][:, 11],
                         "m": res[2][:, 12], "tw": res[2][:, 13], "h": res[2][:, 14],
                         "k": res[2][:, 15], "lamb": res[2][:, 16],
                         "V0": res[3], "jitt": res[4],
                         "lin_trend_out": res[5][0:2], "quad_trend_out": res[5][2:4],
                         "Ndata": res[5][4], "mfit": res[5][5], "rms2": res[5][6],
                         "reduced_chi22": res[5][7], "epoch": res[5][8],
                         "jup_mass": res[5][9:9 + (len(res[5]) - 9) // 2],
                         "jacobi_major_ax": res[5][9 + (len(res[5]) - 9) // 2:],
                         "fit": res[6],"model_data":res[0]}

        # All variables below saves specific results of Fortran run
        self.res = self.res_dict.copy()
        self.jd = self.res_dict["jd"].copy()
        self.model_rvs = self.res_dict["model_rvs"].copy()
        self.rvs = self.res_dict["rvs"].copy()
        self.o_c = self.res_dict["o_c"].copy()
        self.rv_err = self.res_dict["rv_err"].copy()
        self.idset = self.ind_orig.astype(int) - 1
        self.loglik = self.res_dict["loglik"].copy()
        self.reduced_chi2 = self.res_dict["reduced_chi2"].copy()
        self.chi2 = self.res_dict["chi2"].copy()
        self.rms = self.res_dict["rms"].copy()

        if self.arguments[9]:
            self.K = self.__unsort_planets_args(self.res_dict["K"])
            self.P = self.__unsort_planets_args(self.res_dict["P"])
            self.e = self.__unsort_planets_args(self.res_dict["e"])
            self.w = self.__unsort_planets_args(self.res_dict["w"])
            self.M0 = self.__unsort_planets_args(self.res_dict["M0"])
            self.incl = self.__unsort_planets_args(self.res_dict["incl"])
            self.cap0m = self.__unsort_planets_args(self.res_dict["cap0m"])
            len_size = len(self.npl_pos) if self.npl_pos is not None else len(self.res_dict["w_dot"])
            self.omega_dot[0:len_size] = self.__unsort_planets_args(self.res_dict["w_dot"].T[0])
            self.omega_dot_err[0:len_size, :] = np.array([self.__unsort_planets_args(
                                                                self.res_dict["w_dot"].T[1])] * 2).T
            self.t0 = self.__unsort_planets_args(self.res_dict["t0"])
            self.RR = self.__unsort_planets_args(self.res_dict["RR"])
            self.aR = self.__unsort_planets_args(self.res_dict["aR"])
            self.a_new = self.__unsort_planets_args(self.res_dict["a_new"])
            self.m = self.__unsort_planets_args(self.res_dict["m"])
            self.tw = self.__unsort_planets_args(self.res_dict["tw"])
            self.h = self.__unsort_planets_args(self.res_dict["h"])
            self.k = self.__unsort_planets_args(self.res_dict["k"])
            self.lamb = self.__unsort_planets_args(self.res_dict["lamb"])
            self.V0 = self.__sort_array_by_index(self.res_dict["V0"].copy(), self.ind_used)
            self.jitt = self.__sort_array_by_index(self.res_dict["jitt"].copy(), self.ind_used)
            self.lin_trend_out = self.res_dict["lin_trend_out"].copy()
            self.quad_trend_out = self.res_dict["quad_trend_out"].copy()
            self.Ndata = self.res_dict["Ndata"].copy()
            self.mfit = self.res_dict["mfit"].copy()
            self.rms2 = self.res_dict["rms2"].copy()
            self.reduced_chi22 = self.res_dict["reduced_chi22"].copy()
            self.epoch = self.res_dict["epoch"].copy()
            self.mass = self.__unsort_planets_args(list(self.res_dict["jup_mass"]))
            self.a = self.__unsort_planets_args(list(self.res_dict["jacobi_major_ax"]))

            self.planet_params[0:len(res[2]) * 7] = res[2][:, :7, 0].reshape(len(res[2]) * 7)
            self.planet_params_errors[0:len(res[2]) * 7, :] = np.array(
                [res[2][:, :7, 1].reshape(len(res[2]) * 7)] * 2).T
        else:
            self.npl = 1 if self.npl == 0 else self.npl

            self.K = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 0, 0], [0]*self.npl]).T)
            self.P = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 1, 0], [0]*self.npl]).T)
            self.e = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 2, 0], [0]*self.npl]).T)
            self.w = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 3, 0], [0]*self.npl]).T)
            self.M0 = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 4, 0], [0]*self.npl]).T)
            self.incl = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 5, 0], [0]*self.npl]).T)
            self.cap0m = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 6, 0], [0]*self.npl]).T)
            len_size = len(self.npl_pos) if self.npl_pos is not None else len(self.res_dict["w_dot"])
            self.omega_dot[0:len_size] = self.__unsort_planets_args(np.array(
                                                                        self.arguments[17])[:, 7, 0])
            self.omega_dot_err[0:len_size, :] = 0 if self.npl_pos is None else np.zeros((len(self.npl_pos), 2))
            self.t0 = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 8, 0], [0]*self.npl]).T)
            self.RR = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 9, 0], [0]*self.npl]).T)
            self.aR = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 10, 0], [0]*self.npl]).T)
            self.a_new = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 11, 0], [0]*self.npl]).T)
            self.m = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 12, 0], [0]*self.npl]).T)
            self.tw = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 13, 0], [0]*self.npl]).T)
            self.h = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 14, 0], [0]*self.npl]).T)
            self.k = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 15, 0], [0]*self.npl]).T)
            self.lamb = self.__unsort_planets_args(np.array([np.array(self.arguments[17])[:, 16, 0], [0]*self.npl]).T)
            self.V0 = self.__sort_array_by_index(np.array([np.array(self.arguments[15])[:, 0], [0]*self.ndset]).T,
                                                 self.ind_used)
            self.jitt = self.__sort_array_by_index(np.array([np.array(self.arguments[15])[:, 2], [0]*self.ndset]).T,
                                                   self.ind_used)
            self.lin_trend_out = np.array([self.arguments[18][0], 0])
            self.quad_trend_out = np.array([self.arguments[18][2], 0])
            self.Ndata = self.ndata
            self.mfit = 0
            self.rms2 = 0
            self.reduced_chi22 = 0
            self.epoch = self.arguments[18][4]

            if self.npl_pos is not None:
                self.mass = np.zeros(len(self.npl_pos))
                self.a = np.zeros(len(self.npl_pos))
            else:
                self.mass = np.zeros(self.npl)
                self.a = np.zeros(self.npl)

            self.planet_params[0:self.npl * 7] = np.array(self.arguments[17])[:, :7, 0].reshape(self.npl * 7)
            self.planet_params_errors[0:self.npl * 7, :] = 0

        self.model_data = self.res_dict["model_data"].copy().T
        self.fit = self.res_dict["fit"].copy()
        self.model_jd = self.res_dict["fit"].copy().T[0]
        self.model = self.res_dict["fit"].copy().T[1]
        self.rv_quadtr = self.res_dict["quad_trend_out"].copy().T[0]
        self.rv_quadtr_err = np.array([self.res_dict["quad_trend_out"].copy().T[1]]*2).T
        self.rv_lintr = self.res_dict["lin_trend_out"].copy().T[0]
        self.rv_lintr_err = np.array([self.res_dict["lin_trend_out"].copy().T[1]]*2).T

        self.rv_model = RvModelWrap(self.jd, self.rvs, self.rv_err, self.o_c)
        self.stat_array_saved = True
        self.wrms = np.sqrt(np.average(self.o_c ** 2, weights=1 / self.rv_err)) if self.ndset != 0 else 0

        self.offsets[0:len(self.V0.T[0])] = self.V0.T[0]
        self.jitters[0:len(self.jitt.T[0])] = self.jitt.T[0]
        self.params = PlanetParamsWrap(self.offsets, self.jitters, self.planet_params,
                                              self.lin_trend_out[0], self.stellar_mass)

        self.offset_errors[0:len(self.V0.T[1]), :] = np.array([self.V0.T[1]]*2).T
        self.jitter_errors[0:len(self.jitt.T[1]), :] = np.array([self.jitt.T[1]]*2).T
        self.params_errors = PlanetParamsErrorsWrap(self.offset_errors, self.jitter_errors,
                                                           self.planet_params_errors,
                                                           np.array([self.lin_trend_out[1]]*2).T, 0.)

        self.dof = self.Ndata - self.mfit
        self.stat = ParamsWrap(self.params, self.params_errors, dof=self.dof)


# If it is not the main calling, check for compiled versions and import the Fortran
if __name__ != '__main__':
    re = Rvfit()
    re.check_compiled_version()
    del re

    from . import rvmod_for
