#! /usr/bin/env python
import re
import sys
__doc__ = """
Early implementation of a class to setup G09 calculations.
Right now it should be able to write simple input for SP, Opt
and Counterpoise=2 calculations. Very few G09 keywords
(non related to Opt) are implemented,
make it better as you do more G09.

Depends on pmx, as the input for the coords, etc... is not a pdb
but a pmx object. You could lift this requirement by changing
it to biophython, such that we rely on something "robust", but
you definetly need to have an automatized way to deal with the
pdb files.

TODO : Test for input consistencies, e.g. what happen if you ask
for CP calculation but then define two fragn2, does the code break?

"""


CALC_T = {
    "SP": "single point",
    "CP": "Counterpoise",
    "Opt": "Opt",
    "Opt+dihfreeze": "Opt=ModRedun"
}


class G09_calc(object):

    def __init__(self, model, method="MP2", bset="cc-pVTZ",
                 calc="SP", nmol=1, fragn1="DEF", outfn="gau_inp.com",
                 fragn2="DEF", charge=[0], spin=[1], title="Run4Fun",
                 frozen=[], frozen_dih=[], extrakey=""):

        self.memory = "12GB"
        self.sharedproc = "8"
        self.scf = "Tight"
        self.method = method
        self.bset = bset
        self.title = title
        self.calctype = calc
        self.molnum = nmol
        self.outfn = outfn
        self.model = model
        self.resname_1 = fragn1
        self.resname_2 = fragn2
        self.frozen = frozen
        self.frozen_dih = frozen_dih
        self.b_pre_pm3 = False
        self.b_chekpoint = False
        self.extrakey = extrakey
        if calc == "CP":
            self.charge = 3 * [0]
            self.spin = 3 * [1]
        else:
            self.charge = charge[0]
            self.spin = spin[0]
        if calc == "Opt+dihfreeze":
            # self.b_pre_pm3 = True # I don't force it by default
            self.b_chekpoint = True
            # check that we have frozen coords, although
            # we could just simply warn or ignore
            if len(self.frozen) == 0:
                print >> sys.stderr, "Please pass me a list of frozen coords"
                return None

    def _write_header(self, fh, force_pm3=False, geom_check=False):

        try:
            _calc = CALC_T[self.calctype]
        except:
            print >> sys.stderr, "Calculation type {0} " \
                "is unknown to me, what did you mean?".format(
                    self.calctype)
            print >> sys.stderr, "Allowed keywords are the keys "\
                "in the following dictionary"
            print >> sys.stderr, CALC_T
            return 0
        if CALC_T[self.calctype] == "Counterpoise":
            if self.molnum != 2:
                print >> sys.stderr, "I only do counterpoise " \
                    "with 2 molecules."
                return 0
            else:
                _calc = _calc + "=" + str(self.molnum)
        if self.calctype == "Opt":
            _calc = CALC_T[self.calctype]
        elif self.calctype == "Opt+dihfreeze":
            _calc = CALC_T[self.calctype]
        else:
            _calc = ""
        # Header bit
        if force_pm3 or self.calctype != "Opt+dihfreeze":
            fh.write('%NProcShared={0}\n'.format(self.sharedproc))
        fh.write('%chk=molecule\n')
        if force_pm3:
            fh.write('# {0} '.format("PM3"))
        else:
            fh.write('# {0}/{1} '.format(self.method, self.bset))
        fh.write('SCF={0} '.format(self.scf))
        fh.write('{0} '.format(_calc))
        if geom_check:
            fh.write('Geom=Check ')
        if len(self.extrakey) != 0:
            fh.write('{0}'.format(self.extrakey))
        fh.write('\n\n')
        # Title
        fh.write('{0}\n\n'.format(self.title))

        return 1

    def _write_chgspin(self, fh):
        if self.calctype == "CP":
            if len(self.charge) == 3 and len(self.spin) == 3:
                chg_supermol = self.charge[0]
                spin_supermol = self.spin[0]
                chg_f1 = self.charge[1]
                spin_f1 = self.spin[1]
                chg_f2 = self.charge[2]
                spin_f2 = self.spin[2]
                fh.write('{0},{1} {2},{3} {4},{5}\n'
                         .format(chg_supermol,
                                 spin_supermol, chg_f1,
                                 spin_f1, chg_f2, spin_f2))
            else:
                print >> sys.stderr, "I require three pairs "\
                    "of integer charge,spin values for a CP calculation"
                return 0
        else:
            fh.write('{0},{1}\n'.format(self.charge, self.spin))
        return 1

    def _write_coords_CP(self, fh):
        p = re.compile('[a-z]', re.IGNORECASE)
        for c in self.model.chains:
            for res in c.residues:
                if res.resname == self.resname_1:
                    frag_num = 1
                elif res.resname == self.resname_2:
                    frag_num = 2
                for a in res.atoms:
                    atom = p.search(a.name).group()
                    frag = atom + "(Fragment=" + str(frag_num) + ")"
                    fh.write('{0}    {1: 7f} {2: 7f} {3: 7f}\n'.format
                             (frag, a.x[0], a.x[1], a.x[2]))
        # write a blank line afert coords
        fh.write('\n')
        return 1

    def _write_coords(self, fh):
        p = re.compile('[a-z]', re.IGNORECASE)
        for c in self.model.chains:
            for res in c.residues:
                for a in res.atoms:
                    atom = p.search(a.name).group()
                    frag = atom
                    fh.write('{0}    {1: 7f} {2: 7f} {3: 7f}\n'
                             .format(frag, a.x[0], a.x[1], a.x[2]))
            # write a blank line afert coords
        fh.write('\n')
        return 1

    def _write_freeze_dih(self, fh):

        for i,f in enumerate(self.frozen):
            if len(f) != 4:
                print >> sys.stderr, "Freeze_dihedrals requires a quadruplet"
                print >> sys.stderr, "I got: ",  f
                sys.exit(1)

            at_list = []
            for ind in f:
                at_list.append(self.model.atoms[ind-1])
            if len(self.frozen_dih) == 0 :
                dih = at_list[0].dihedral(at_list[1], at_list[2], at_list[3],
                                          degree=True)
            elif len(self.frozen_dih) == len(self.frozen):
                dih = self.frozen_dih[i]
            else:
                print sys.stderr, "Unmatched number of dihedrals to be frozen",
                print sys.stderr, "and matching reference positions"
                print sys.stderr, self.frozen_dih, self.frozen
                sys.exit(1)

            fh.write("{0:3d} {1:3d} {2:3d} {3:3d} {4:8.4f} F\n".format(
                f[0], f[1], f[2], f[3], dih))

    def write(self, fn="default"):
        if fn == "default":
            fn = self.outfn
        with open(fn, 'w') as fh:
            if self.b_pre_pm3 and self.calctype == "Opt+dihfreeze":
                self._write_header(fh, force_pm3=True)
                self._write_chgspin(fh)
                self._write_coords(fh)
                self._write_freeze_dih(fh)
                fh.write("\n\n")
                fh.write("--Link1--\n")
                self._write_header(fh, geom_check=True)
                self._write_chgspin(fh)
                fh.write("\n")
                self._write_freeze_dih(fh)
                fh.write("\n\n")
                return 1

            elif self.calctype == "Opt+dihfreeze":
                self._write_header(fh, force_pm3=False)
                self._write_chgspin(fh)
                self._write_coords(fh)
                self._write_freeze_dih(fh)
                fh.write("\n\n")
                return 1

            else:
                self._write_header(fh)
                self._write_chgspin(fh)

            if self.calctype == "CP":
                self._write_coords_CP(fh)
                return 1
            elif self.calctype != "Opt+dihfreeze":
                self._write_coords(fh)
                return 1
