#! /usr/bin/env python

"""
    Docstring
"""
# import matplotlib.pyplot as plt
# from matplotlib import rc
import pmx
import gromacs
import glob
import os
import shutil
import g09
import sys
import tarfile
import itertools

HEAD_TOP_FILE = """
#include "parmbsc0-star-ildn.ff/forcefield.itp"
#include "dhn_gmx.itp"
#ifdef DIHE_RES
"""
BOTTOM_TOP_FILE = """
#endif
#include "parmbsc0-star-ildn.ff/spce.itp"

[ system ]
dafunk
[ molecules ]
DNA           1
"""


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
        return


def mdruncleanup():
    to_delete = ["*.edr", "*.xtc", "*.log", "mdout.mdp",
                 "*.trr", "em.tpr", "#*"]
    for exp in to_delete:
        for filename in glob.glob(exp):
            os.remove(filename)


def ipt_pdb_em_cleanup():
    to_delete = ["dih_restr_*.itp", "out_em_*.pdb"]
    for exp in to_delete:
        for filename in glob.glob(exp):
            os.remove(filename)


def get_dihe_from_atname(model, at):
    """
        Return the dihedral formed by four atoms
    """
    if len(at) != 4:
        print >> sys.stderr, "Need 4 atoms to compute a dihedral"
        return None
    at_d = []
    for a in at:
        at_d.append(model.residues[0][a])

    return 22


def get_dihe(model, at):
    """
        Return the dihedral formed by four atoms
    """
    if len(at) != 4:
        print >> sys.stderr, "Need 4 atoms to compute a dihedral"
        return None

    return at[0].dihedral(at[1], at[2], at[3], degree=True)


def new_dihe_at_id(redmod, at_quadruplets):
    at_ind = []
    for i, ats in enumerate(at_quadruplets):
        if len(ats) != 4:
            print >> sys.stderr, "Need 4 atoms to compute a dihedral"
            return None
        l_at_ind = []
        for j, a in enumerate(ats):
            l_at_ind.append(redmod.residues[0][a].id)
        at_ind.append(l_at_ind)

    return at_ind


def write_dih_restr(mod, at_quadruplets, dihedrals, itpf,
                    fct=500.0):
    """ Here I should iterate first over the number
        of different dihedrals
    """
    if len(at_quadruplets) != len(dihedrals):
        print >> sys.stderr, "Number of atom quadruplets and dihedrals",
        print >> sys.stderr, " does not match",\
                             len(at_quadruplets), len(dihedrals)
        sys.exit(1)
    # make it indep of number of dih
    suff = ""
    for d in dihedrals:
        suff += "_" + str(d)
    itp_outfn = itpf + suff + ".itp"
    with open(itp_outfn, "w") as fh:
        fh.write("[dihedral_restraints]\n")
        for i, ats in enumerate(at_quadruplets):
            if len(ats) != 4:
                print >> sys.stderr, "Need 4 atoms to compute a dihedral"
                return None
            at_ind = []
            for j, a in enumerate(ats):
                at_ind.append(mod.residues[0][a].id)
            fh.write("{0:3d} {1:3d} {2:3d} {3:3d}  1"
                     "{4: 8.4f}  0  {5:f}\n".format(
                      at_ind[0], at_ind[1], at_ind[2], at_ind[3],
                      float(dihedrals[i]), float(fct)))
    return itp_outfn


def do_energy_min(pdb, itpfn, outpdb="out.pdb"):
    # write the topology file
    with open("topol_em.top", "w") as fh_t:
        fh_t.write(HEAD_TOP_FILE)
        fh_t.write("#include " + '"' + itpfn + '"')
        fh_t.write(BOTTOM_TOP_FILE)
    # TODO function that checks if files are present
    gromacs.grompp(f="em.mdp", c=pdb, p="topol_em.top",
                   o="em.tpr", maxwarn="99")
    gromacs.mdrun(s="em.tpr", c=outpdb)


def reduce_system(pdb_file):
    """ Get rid of the sugar and place a H instead
    """
    full_model = pmx.Model(pdb_file)
    red_model = pmx.Model()
    at_counter = 1
    for c in full_model.chains:
        for r in c.residues:
            for a in r.atoms:
                if "'" not in a.name and "H5T" != a.name and "H3T" != a.name:
                    a.id = at_counter
                    red_model.atoms.append(a)
                    at_counter += 1
                if "C1'" == a.name:
                    a.name = "HN1"
                    a.symbol = "H"
                    a.id = at_counter
                    red_model.atoms.append(a)
                    at_counter += 1

    red_model.make_chains()
    red_model.make_residues()

    return red_model


def do_gen_ff_prof(inpdb, ats_dih, dih_freq=None, fct=500.0,
                   itpf="dih_rest"):

    if dih_freq is None:
        # use a default of 20
        for at_dih in ats_dih:
                dih_freq.append(20)

    model = pmx.Model(inpdb)

    directory = "g09_input_files"
    if not os.path.exists(directory):
        os.makedirs(directory)

    dih_values = []
    for i, u in enumerate(xrange(len(dih_freq))):
        dih_values.append([x for x in xrange(0, 360, dih_freq[i])])

    for dih_combination in itertools.product(*dih_values):
        itpf_d = write_dih_restr(model, ats_dih, dih_combination, itpf, fct)
        suff = ""
        for d in dih_combination:
            suff += "_" + str(d)
        out_pdb_em = "out_em" + suff + ".pdb"
        fout_gau = "dh" + suff + ".com"
        if itpf_d is not None:
            do_energy_min(inpdb, itpf_d, outpdb=out_pdb_em)
            mdruncleanup()
            reduced_model = reduce_system(out_pdb_em)
            frozen_dih = new_dihe_at_id(reduced_model, ats_dih)
            gau_obj = g09.G09_calc(reduced_model, method="B3LYP",
                                   bset="6-31+G(d,p)",
                                   calc="Opt+dihfreeze", nmol=1,
                                   outfn=fout_gau, title="dafunk",
                                   frozen=frozen_dih,
                                   frozen_dih = list(dih_combination),
                                   extrakey="SCRF=(SMD,solvent=water)")
            if gau_obj.write():
                dst_file = os.path.join(directory, fout_gau)
                if os.path.exists(dst_file):
                    os.remove(dst_file)
                shutil.move(fout_gau, directory)
            else:
                print >> sys.stderr, "Problems writing gaussian"

        ipt_pdb_em_cleanup()

    # compress the folder
    make_tarfile(directory + ".tar.bz2", directory)
    shutil.rmtree(directory)

if __name__ == "__main__":
    at_dih1 = ["C4", "C5", "C55", "OH5"]
    at_dih2 = ["C5", "C55", "OH5", "HO5"]
    atoms_dih = [at_dih1, at_dih2]
    # dfreq is the frequency to scan each dihedral
    dfreq = [10, 10]
    inpdb = "dhn_centered_box.pdb"
    do_gen_ff_prof(inpdb, atoms_dih, dih_freq=dfreq,
                   fct=500.0, itpf="dih_restr")
