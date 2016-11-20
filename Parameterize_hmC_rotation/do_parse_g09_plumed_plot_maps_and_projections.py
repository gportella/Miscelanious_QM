#! /usr/bin/env python

import numpy as np
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import math

EV2H = 1./27.211385050
H2KJMOL = 2625.5

def lastscf(fn):
    opt_energy = 0.0
    with open(fn) as f:
        for l in f:
            if "SCF Done:" in l.rstrip():
                opt_energy = float(l.rstrip("\n").\
                                   split("=")[1].split()[0])

    if opt_energy >= 0.0:
        return None
    else:
        return opt_energy


def change_range(ang_degrees):
    new_range_deg = ang_degrees if ang_degrees <= 180 else ang_degrees - 360
    return new_range_deg


def e_vs_d(dih_freq):
    d1 = []
    d2 = []
    Edd = []
    dih_values = []
    for i, u in enumerate(xrange(len(dih_freq))):
        dih_values.append([x for x in xrange(-180, 180, dih_freq[i])])
    for dih_comb in itertools.product(*dih_values):
        # I want the order switched wrt itertools product, a bit hacky
        # this is due to the way the plotting works
        # data remains the same
        new_d1 = dih_comb[1] % 360
        new_d2 = dih_comb[0] % 360
        fn = "dh" + "_" + str(new_d1) +\
        "_" + str(new_d2) + ".out"
        E = lastscf(fn)*H2KJMOL
        d1.append(dih_comb[1])
        d2.append(dih_comb[0])
        Edd.append(E)
    d1d2 = np.column_stack((np.array(d1), np.array(d2)))
    g = np.column_stack((d1d2, (np.array(Edd))- min(Edd)))
    return g


def heat_plot(data, title_plot="DFT", do_write=True, outpng="dih_hmC.png"):
    vecfunc = np.vectorize(change_range)
    d1 = vecfunc(data[:,0]).reshape(len(data[:,0]),1)
    d2 = vecfunc(data[:,1]).reshape(len(data[:,0]),1)
    g = np.hstack((d1, d2, data[:,2].reshape((len(data[:,0]),1))))
    N = int(len(data[:,2])**.5)
    z = g[:,2].reshape(N, N)
    fig, ax = plt.subplots(figsize=(9,9), dpi=120)
    im = plt.imshow(z, extent=(min(d1)[0], max(d1)[0], min(d2)[0], max(d2)[0]), origin="lower", cmap=cm.CMRmap)
    ax.set_xlabel("C4-C5-C55-OH5 (deg)")
    ax.set_ylabel("C5-C55-OH5-HO5 (deg)")
    ax.set_title(title_plot, fontsize=14)
    plt.clim(0,45)
    cb = fig.colorbar(im, fraction=0.046, pad=0.04, ax=ax)
    cb.set_label("Energy (kJ/mol)")
    #cs = plt.contour(z, 25, extent=(-180, 180, -180, 180), origin="lower",  cmap=cm.CMRmap)
    cs = plt.contour(z, 25, extent=(min(d1)[0], max(d1)[0], min(d2)[0], max(d2)[0]),
                     origin="lower", linewidths=.5,  cmap=cm.CMRmap)
    plt.clim(0,1000)
    if do_write:
        plt.savefig(outpng)
    else:
        plt.show()


def project(data, axis=1):
    gg_tup = tuple(map(tuple, data))
    d = {}
    if axis != 0 and axis != 1:
        print "Axis can only be 0 or 1"
        return None
    for el in gg_tup:
        d.setdefault(el[axis], []).append((el[2]))
    deg = []
    ener = []
    for k in d:
        norm = 0
        for e_s in d[k]:
            norm += math.exp(-e_s)
        s = 0
        for e_s in d[k]:
            s += e_s*math.exp(-e_s)/norm
        deg.append(k)
        ener.append(s)
    emin = min(ener)
    ener = [x - emin for x in ener]
    deg = np.array(deg)
    ener = np.array(ener)
    projected = np.column_stack((deg, ener))
    return projected


def do_plot(data, do_write=False, outpng="ff_vs_qm.png"):
    count = 0
    f, ax = plt.subplots(len(data),1, sharex=True,
                         figsize=(9,9), dpi=120)
    f.subplots_adjust(hspace=0.2)
    labels = ["FF", "DFT"]
    for datasets in data:
        count += 1
        ax[count-1].grid(b=True, which='major',
                         color='gray', linestyle='--', linewidth=0.5)
        ax[count-1].set_axisbelow(True)
        for c_counter, dataset in enumerate(datasets[0]):
            thelabel = labels[c_counter]
            dihedral = datasets[1]
            ax[count-1].set_title(dihedral)
            for j in range(dataset.shape[1]-1):
                ax[count-1].plot(dataset[:, 0],dataset[:,j+1],
                     linewidth=1.2, label=thelabel)

        ax[count-1].legend(ncol=4).get_frame().set_linewidth(0.1)
        plt.xlabel("Dihedral (deg)")
        ax[count-1].set_ylabel("Energy (kJ/mol)")

    if do_write:
        plt.savefig(outpng)
    else:
        plt.show()


def mod_360(val):
    return val % 360


def print_xy(p, label, fname="out.dat"):
    p[:, 0] = map(mod_360, p[:, 0])
    p.view("f8, f8").sort(axis=0)
    with open(fname, "w") as fh:
        np.savetxt(fh, p, fmt="%f")


dfreq = [10, 10]
gg = e_vs_d(dfreq)
gg[:,0] = map(mod_360, gg[:, 0])
gg[:,1] = map(mod_360, gg[:, 1])
gg.view("f8, f8, f8").sort(axis=0)
np.savetxt("2d_dft_deg.dat", gg, fmt="%f")

# heat_plot(gg, title_plot="DFT", do_write=True, outpng="dih_hmC_dft.png")

ff_2d_rad = np.genfromtxt("fe_2d.dat", usecols=(0, 1, 2))
d1 = np.rad2deg(ff_2d_rad[:, 0])
d1[:] = map(mod_360, d1[:])
d2 = np.rad2deg(ff_2d_rad[:, 1])
d2[:] = map(mod_360, d2[:])
d1d2 = np.column_stack((d1, d2))
ff_2d = np.column_stack((d1d2, ff_2d_rad[:, 2]))
ff_2d.view("f8, f8, f8").sort(axis=0)
np.savetxt("2d_fe_deg.dat", ff_2d, fmt="%f")

p0 = project(gg, axis=0)
p1 = project(gg, axis=1)
p0.view("f8, f8").sort(axis=0)
p1.view("f8, f8").sort(axis=0)
fe_one = np.genfromtxt("fe_one.dat")
fe_one[:,0] *= 57.295779513
fe_two = np.genfromtxt("fe_two.dat")
fe_two[:,0] *= 57.295779513

print_xy(p0, "DFT C4-C5-C55-OH5", "dft_cc.dat")
print_xy(p1, "DFT C5-C55-OH5-HO5", "dft_co.dat")
print_xy(fe_one, "FF C4-C5-C55-OH5", "ff_cc.dat")
print_xy(fe_two, "FF C5-C55-OH5-HO5", "ff_co.dat")

# data = [([fe_one, p0], "C4-C5-C55-OH5"),
#        ([fe_two, p1], "C5-C55-OH5-HO5")]
# do_plot(data, do_write=True, outpng="ff_vs_qm.png")
