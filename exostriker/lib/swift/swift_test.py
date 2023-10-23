import numpy as np

import swift

iuflg = 1
mstar0 = 1.0
nbod = 2
frho3 = 1.0
planets_params = [[0.010956056701692769, 1.286350183286512, 0.12495909431420299, 90.0, 8.761185986070139,
                   0.0, 236.93745955638647],
                  [0.01207470374287773, 4.267216547236886, 0.19126053191013156, 90.0, 10.005317383453992,
                   0.0, 221.47361992077785]]

print("results")
res = swift.geninit_j3_in_days(iuflg, mstar0, nbod, frho3, planets_params)
# print(res)

t0 = 0.00
tstop = 3652500.0
dt = 5.0
dtout = 36525.0
dtdump = 365250.0
lflg = [0, 1, 1, 1, 1, 0]
rmin = 0.0001
rmax = 50.0
rmaxu = 50.0
qmin = -1.
lclose = "T"
outfile = "bin.dat"
fopenstat = "unknown"
coordinates = res[1]

mtiny = 1e-40
ll = 1000

res = swift.swift_symba5_j(t0, tstop, dt, dtout, dtdump, lflg,
                            rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat,
                            nbod, coordinates, mtiny, (tstop-t0)/(dtdump-t0))
# res = swift.swift_mvs_j(t0, tstop, dt, dtout, dtdump, lflg,
#                         rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat,
#                         nbod, coordinates, (tstop-t0)/(dtdump-t0))
#res = swift.swift_mvs_j_gr(t0, tstop, dt, dtout, dtdump, lflg,
#                           rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat,
#                           nbod, coordinates, ll, (tstop-t0)/(dtdump-t0))
print(res)

nbody = 2

res = swift.follow_symba2(t0, tstop, dt, dtout, dtdump, lflg,
                          rmin, rmax, rmaxu, qmin, lclose, outfile, fopenstat,
                          nbod, coordinates, nbody, (tstop-t0)/dtout+1)

print(res)
