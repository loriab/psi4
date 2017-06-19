#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Module with utility functions that act on molecule objects."""
from __future__ import absolute_import
import re
import math

from psi4 import core
from psi4.driver.p4util import constants
from psi4.driver.inputparser import process_pubchem_command, pubchemre
from psi4.driver import qcdb
from psi4.driver.qcdb import periodictable

#NUMBER = r'((?:[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?)|(?:[-+]?\d+\.\d*(?:[DdEe][-+]?\d+)?))'
NUMBER = r'((?:[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?)|(?:[-+]?\d+\.\d*(?:[DdEe][-+]?\d+)?)|(?:[-+]?\d+(?:[DdEe][-+]?\d+)?))'
SEP = r'[\t ,]+'


def extract_clusters(mol, ghost=True, cluster_size=0):
    """Function to return all subclusters of the molecule *mol* of
    real size *cluster_size* and all other atoms ghosted if *ghost*
    equals true, all other atoms discarded if *ghost* is false. If
    *cluster_size* = 0, returns all possible combinations of cluster size.

    """
    # How many levels of clusters are possible?
    nfrag = mol.nfragments()

    # Initialize the cluster array
    clusters = []

    # scope the arrays
    reals = []
    ghosts = []

    # counter
    counter = 0

    # loop over all possible cluster sizes
    for nreal in range(nfrag, 0, -1):

        # if a specific cluster size size is requested, only do that
        if (nreal != cluster_size and cluster_size > 0):
            continue

        # initialize the reals list
        reals = []

        # setup first combination [3,2,1] lexical ordering
        # fragments indexing is 1's based, bloody hell
        for index in range(nreal, 0, -1):
            reals.append(index)

        # start loop through lexical promotion
        while True:

            counter = counter + 1

            # Generate cluster from last iteration
            if (ghost):
                ghosts = []
                for g in range(nfrag, 0, -1):
                    if (g not in reals):
                        ghosts.append(g)
                clusters.append(mol.extract_subsets(reals, ghosts))
            else:
                clusters.append(mol.extract_subsets(reals))

            # reset rank
            rank = 0

            # look for lexical promotion opportunity
            # i.e.: [4 2 1] has a promotion opportunity at
            #   index 1 to produce [4 3 1]
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            # do the promotion
            reals[rank] = reals[rank] + 1

            # demote the right portion of the register
            val = 1
            for k in range(nreal - 1, rank, -1):
                reals[k] = val
                val = val + 1

            # boundary condition is promotion into
            # [nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break

    return clusters


def extract_cluster_indexing(mol, cluster_size=0):
    """Function to returns a LIST of all subclusters of the molecule *mol* of
    real size *cluster_size*. If *cluster_size* = 0, returns all possible
    combinations of cluster size.

    """
    import copy

    # How many levels of clusters are possible?
    nfrag = mol.nfragments()

    # Initialize the cluster array
    clusters = []

    # scope the arrays
    reals = []

    # counter
    counter = 0

    # loop over all possible cluster sizes
    for nreal in range(nfrag, 0, -1):

        # if a specific cluster size size is requested, only do that
        if (nreal != cluster_size and cluster_size > 0):
            continue

        # initialize the reals list
        reals = []

        # setup first combination [3,2,1] lexical ordering
        # fragments indexing is 1's based, bloody hell
        for index in range(nreal, 0, -1):
            reals.append(index)

        # start loop through lexical promotion
        while True:

            counter = counter + 1

            # Generate cluster from last iteration
            clusters.append(copy.deepcopy(reals))

            # reset rank
            rank = 0

            # look for lexical promotion opportunity
            # i.e.: [4 2 1] has a promotion opportunity at
            #   index 1 to produce [4 3 1]
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1
                    break

            # do the promotion
            reals[rank] = reals[rank] + 1

            # demote the right portion of the register
            val = 1
            for k in range(nreal - 1, rank, -1):
                reals[k] = val
                val = val + 1

            # boundary condition is promotion into
            # [nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break

    return clusters


def molecule_set_attr(self, name, value):
    """Function to redefine __setattr__ method of molecule class."""
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)
    if isvar:
        fxn = object.__getattribute__(self, "set_variable")
        fxn(name, value)
        return

    object.__setattr__(self, name, value)


def molecule_get_attr(self, name):
    """Function to redefine __getattr__ method of molecule class."""
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)

    if isvar:
        fxn = object.__getattribute__(self, "get_variable")
        return fxn(name)

    return object.__getattribute__(self, name)


def BFS(self):
    """Perform a breadth-first search (BFS) on the real atoms
    in molecule, returning an array of atom indices of fragments.
    Relies upon van der Waals radii and so faulty for close
    (esp. hydrogen-bonded) fragments. Original code from
    Michael S. Marshall.

    """
    vdW_diameter = {
        'H':  1.001 / 1.5,
        'HE': 1.012 / 1.5,
        'LI': 0.825 / 1.5,
        'BE': 1.408 / 1.5,
        'B':  1.485 / 1.5,
        'C':  1.452 / 1.5,
        'N':  1.397 / 1.5,
        'O':  1.342 / 1.5,
        'F':  1.287 / 1.5,
        'NE': 1.243 / 1.5,
        'NA': 1.144 / 1.5,
        'MG': 1.364 / 1.5,
        'AL': 1.639 / 1.5,
        'SI': 1.716 / 1.5,
        'P':  1.705 / 1.5,
        'S':  1.683 / 1.5,
        'CL': 1.639 / 1.5,
        'AR': 1.595 / 1.5}

    Queue = []
    White = range(self.natom())  # untouched
    Black = []  # touched and all edges discovered
    Fragment = []  # stores fragments

    start = 0  # starts with the first atom in the list
    Queue.append(start)
    White.remove(start)

    # Simply start with the first atom, do a BFS when done, go to any
    #   untouched atom and start again iterate until all atoms belong
    #   to a fragment group
    while White or Queue:                # Iterates to the next fragment
        Fragment.append([])

        while Queue:                     # BFS within a fragment
            for u in Queue:              # find all white neighbors to vertex u
                for i in White:
                    dist = constants.bohr2angstroms * math.sqrt(
                           (self.x(i) - self.x(u)) ** 2 +
                           (self.y(i) - self.y(u)) ** 2 +
                           (self.z(i) - self.z(u)) ** 2)
                    if dist < vdW_diameter[self.symbol(u)] + \
                       vdW_diameter[self.symbol(i)]:
                        Queue.append(i)  # if you find you, put in the queue
                        White.remove(i)  # & remove it from the untouched list
            Queue.remove(u)              # remove focus from Queue
            Black.append(u)
            Fragment[-1].append(int(u))  # add to group (0-indexed)
            Fragment[-1].sort()          # preserve original atom ordering

        if White:                        # can't move White -> Queue if empty
            Queue.append(White[0])
            White.remove(White[0])

    return Fragment


def dynamic_variable_bind(cls):
    """Function to dynamically add extra members to
    the core.Molecule class.

    """
    cls.__setattr__ = molecule_set_attr
    cls.__getattr__ = molecule_get_attr

    cls.BFS = BFS


dynamic_variable_bind(core.Molecule)  # pass class type, not class instance


#
# Define geometry to be used by PSI4.
# The molecule created by this will be set in options.
#
# geometry("
#   O  1.0 0.0 0.0
#   H  0.0 1.0 0.0
#   H  0.0 0.0 0.0
#
def geometry(geom, name="default"):
    """Function to create a molecule object of name *name* from the
    geometry in string *geom*. Permitted for user use but deprecated
    in driver in favor of explicit molecule-passing. Comments within
    the string are filtered.

    """
    geom = pubchemre.sub(process_pubchem_command, geom)
    molecule = init_psi4_molecule_from_any_string(geom, name=name)
    molecule.set_name(name)

    activate(molecule)

    return molecule


def activate(mol):
    """Function to set molecule object *mol* as the current active molecule.
    Permitted for user use but deprecated in driver in favor of explicit
    molecule-passing.

    """
    core.set_active_molecule(mol)


def init_psi4_molecule_from_any_string(mol_str, **kwargs):



#    #if not has_efp:

    # Figure out how and if we will parse the Molecule adata
    mname = kwargs.pop("name", "default")
    dtype = kwargs.pop("dtype", "psi4").lower()
    if mol_str is not None:

        print('<<< QMMOL', mol_str, '>>>')
        mol_str = filter_comments(mol_str)

        if dtype == 'psi4':

            # Psi4 defaults
            mol_init = {'name': mname,
                        'units': 'Angstrom',
                        'input_units_to_au': 1.0 / constants.bohr2angstroms,
                        'fragments': [],
                        'fragment_types': [],
                        'fragment_charges': [],
                        'fragment_multiplicities': []}

            # handle units, com, orient, symm
            mol_str, univ = filter_universals(mol_str)
            mol_init.update(univ)
            #mol_str, univ = filter_universals(mol_str, seed=mol_init)

            mol_str, efp_init = filter_libefp(mol_str, confer=mol_init)
            if efp_init:
                print('<<< core.EFP INTO', efp_init, '>>>')
                core.efp_init()
                efp = core.get_active_efp()
                efp.construct_from_pydict(efp_init)
                mol_init['fix_com'] = True
                mol_init['fix_orientation'] = True
                mol_init['fix_symmetry'] = 'c1'

            mol_str, mints_init = filter_mints(mol_str)
            mol_init.update(mints_init)

            mol_init = reconcile_cgmp(mol_init)
            print('\nINTO core.Molecule.from_dict <<<', mol_init, '>>>\n')
            molecule = core.Molecule.from_dict(mol_init)
            molecule.update_geometry()
            molecule.print_out()
            return molecule

            #has_efp, geom = filter_libefp(geom)

        if dtype == "oldpsi4":
            qcdbmol = qcdb.Molecule(mol_str)
            qcdbmol.update_geometry()
            qcdbmol.everything()
            pymol = qcdbmol.export_for_libmints()
            molecule = core.Molecule.construct_from_pydict(pymol)
            #molecule.set_name(mname)
            molecule.print_out()

            activate(molecule)

            print('RRRR'
            ,molecule.natom()
            ,molecule.nfragments()
            ,molecule.com_fixed()
            ,molecule.orientation_fixed()
            ,molecule.charge_specified()
            ,molecule.multiplicity_specified()
            ,molecule.symmetry_specified()
            )

            return molecule

        #    self._molecule_from_string_psi4(mol_str)
        #elif dtype == "numpy":
        #    frags = kwargs.pop("frags", [])
        #    self._molecule_from_numpy(mol_str, frags)
        #elif dtype == "json":
        #    self._molecule_from_json(mol_str)
        else:
            raise KeyError("Molecule: dtype of %s not recognized.")
    else:
        # In case a user wants to build one themselves
        pass


def reconcile_cgmp(pydict):
    def apply_default(lst, default):
        return [default if (c is None) else c for c in lst]

    def reconcile_charges(chg, fchg):
        """expects None as placeholder for not specified"""

        naive_fchg = apply_default(fchg, 0)

        if chg is None:
            # Sys: None, Frag: [-1, 1] --> Sys: 0, Frag: [-1, 1]
            # Sys: None, Frag: [None, -2] --> Sys: -2, Frag: [0, -2]
            return sum(naive_fchg), naive_fchg
        else:
            if None in fchg:
                # Sys: 0, Frag: [None, 1] --> Sys: 0, Frag: [-1, 1]
                # Sys: 2, Frag: [1, None, None] --> Sys: 2, Frag: [1, 1, 0]
                missing_frag_chg = chg - sum(naive_fchg)
                first = fchg.index(None)
                tmp = fchg[:]
                tmp[first] = missing_frag_chg
                return chg, apply_default(tmp, 0)
            else:
                if chg == sum(fchg):
                    # Sys: 0, Frag: [-1, 1] --> Sys: 0, Frag: [-1, 1]
                    return chg, fchg
                else:
                    # Sys: 0, Frag: [-2, 1] --> Irreconcilable!
                    raise ValidationError('System charge: {} not reconcilable with fragment charges: {}'.format(chg, fchg))

    def reconcile_multiplicities(mult, fmult):

        def high_spin_sum(mult_list):
            mm = 1
            for m in mult_list:
                mm += m - 1
            return mm

        naive_fmult = apply_default(fmult, 1)

        if mult is None:
            return high_spin_sum(naive_fmult), naive_fmult
        else:
            if None in fmult:
                missing_frag_mult = mult - high_spin_sum(naive_fmult) + 1
                first = fmult.index(None)
                tmp = fmult[:]
                tmp[first] = missing_frag_mult
                return mult, apply_default(tmp, 1)
            else:
                if mult == high_spin_sum(fmult):
                    return mult, fmult
                else:
                    # TODO: could be ok, could be not
                    raise ValidationError('Spin Arithmetic')

    chg = pydict.get('system_charge', None)
    mult = pydict.get('system_multiplicity', None)
    fchg = pydict['fragment_charges']
    fmult = pydict['fragment_multiplicities']
    # TODO look at the atoms and see if chgmult appropriate

    if (chg is not None) or (len(fchg) == 1 and fchg[0] is not None):
        chg_specified = True
    else:
        chg_specified = False

    if (mult is not None) or (len(fmult) == 1 and fmult[0] is not None):
        mult_specified = True
    else:
        mult_specified = False

    chg, fchg = reconcile_charges(pydict.get('system_charge'), pydict['fragment_charges'])
    mult, fmult = reconcile_multiplicities(pydict.get('system_multiplicity'), pydict['fragment_multiplicities'])

    if chg_specified:
        pydict['system_charge'] = chg
    pydict['fragment_charges'] = fchg
    if mult_specified:
        pydict['system_multiplicity'] = mult
    pydict['fragment_multiplicities'] = fmult

    return pydict



def filter_pubchem(mol_str):
    pass
    pubchemerror = re.compile(r'^\s*PubchemError\s*$', re.IGNORECASE)
    pubcheminput = re.compile(r'^\s*PubchemInput\s*$', re.IGNORECASE)


def filter_comments(mol_str):
    comment = re.compile(r'(^|[^\\])#.*')
    mol_str = re.sub(comment, '', mol_str)
    return mol_str


def filter_universals(mol_str, seed=None):

    com = re.compile(r'\A(no_com|nocom)\Z', re.IGNORECASE)
    orient = re.compile(r'\A(no_reorient|noreorient)\Z', re.IGNORECASE)
    bohrang = re.compile(r'\Aunits?[\s=]+((?P<ubohr>(bohr|au|a.u.))|(?P<uang>(ang|angstrom)))\Z', re.IGNORECASE)
    symmetry = re.compile(r'\Asymmetry[\s=]+(?P<pg>\w+)\Z', re.IGNORECASE)

    def process_com(matchobj):
        univ['fix_com'] = True
        return ''

    def process_orient(matchobj):
        univ['fix_orientation'] = True
        return ''

    def process_bohrang(matchobj):
        if matchobj.group('uang'):
            univ['units'] = 'Angstrom'
            univ['input_units_to_au'] = 1.0 / constants.bohr2angstroms
        elif matchobj.group('ubohr'):
            univ['units'] = 'Bohr'
            univ['input_units_to_au'] = 1.0
        return ''

    def process_symmetry(matchobj):
        univ['fix_symmetry'] = matchobj.group('pg').lower()
        return ''

    reconstitute = []
    univ = {}
    if seed:
        univ = seed

    for line in mol_str.split('\n'):
        line = re.sub(com, process_com, line.strip())
        line = re.sub(orient, process_orient, line)
        line = re.sub(bohrang, process_bohrang, line)
        line = re.sub(symmetry, process_symmetry, line)
        if line:
            reconstitute.append(line)

    return '\n'.join(reconstitute), univ


def filter_mints(mol_str, seed=None):

    fragment_re = re.compile(r'^\s*--\s*$', re.MULTILINE)
    ATOM = r'(?:(?P<gh1>@)|(?P<gh2>Gh\())?(?P<label>(?P<symbol>[A-Z]{1,3})(?:(_\w+)|(\d+))?)(?(gh2)\))(?:@(?P<mass>\d+\.\d+))?'
    CHG = r'(?P<chg>-?\d+)'
    MULT = r'(?P<mult>\d+)'
    cgmp = re.compile(r'\A' + CHG + SEP + MULT + r'\Z')
    atom_cart = re.compile(r'\A' + ATOM + SEP + NUMBER + SEP + NUMBER + SEP + NUMBER + r'\Z', re.IGNORECASE)
    atom_zmat1 = re.compile(r'\A' + ATOM + r'\Z', re.IGNORECASE)
    atom_zmat2 = re.compile(r'\A' + ATOM + SEP + r'(\d+)' + SEP + NUMBER + r'\Z', re.IGNORECASE)
    atom_zmat3 = re.compile(r'\A' + ATOM + SEP + r'(\d+)' + SEP + NUMBER
                                         + SEP + r'(\d+)' + SEP + NUMBER + r'\Z', re.IGNORECASE)
    atom_zmat4 = re.compile(r'\A' + ATOM + SEP + r'(\d+)' + SEP + NUMBER
                                         + SEP + r'(\d+)' + SEP + NUMBER
                                         + SEP + r'(\d+)' + SEP + NUMBER + r'\Z', re.IGNORECASE)
    #frag = re.compile(r'^\s*--\s*$')
    #variable = re.compile(r'^\s*(\w+)\s*=\s*(-?\d+\.\d+|-?\d+\.|-?\.\d+|-?\d+|tda)\s*$', re.IGNORECASE)
    #ghost = re.compile(r'@(.*)|Gh\((.*)\)', re.IGNORECASE)

    def process_system_cgmp(matchobj):
        mints_init['system_charge'] = int(matchobj.group('chg'))
        mints_init['system_multiplicity'] = int(matchobj.group('mult'))
        return ''

    def filter_fragment(mol_str):

        def process_fragment_cgmp(matchobj):
            mints_init['fragment_charges'].append(int(matchobj.group('chg')))
            mints_init['fragment_multiplicities'].append(int(matchobj.group('mult')))
            return ''

        def process_atom_spec(matchobj):
            atom_symbol = matchobj.group('symbol').upper()
            if atom_symbol not in periodictable.el2z:
                raise ValidationError("""Illegal atom symbol in geometry specification: {}""".format(atom_symbol))

            if matchobj.group('gh1') or matchobj.group('gh2'):
                atom_Z = 0.0
                atom_charge = 0.0
                atom_ghosted = True
            else:
                atom_Z = float(periodictable.el2z[atom_symbol])
                atom_charge = float(atom_Z)
                atom_ghosted = False

            if matchobj.group('mass'):
                atom_mass = float(matchobj.group('mass'))
            else:
                atom_mass = periodictable.el2mass[atom_symbol]

            atom_init = {}
            atom_init['Z'] = atom_Z
            atom_init['charge'] = atom_charge
            atom_init['ghosted'] = atom_ghosted
            atom_init['symbol'] = atom_symbol
            atom_init['label'] = matchobj.group('label').upper()
            atom_init['mass'] = atom_mass
            return atom_init

        def process_atom_cart(matchobj):

            atom_init = process_atom_spec(matchobj)
            atom_init['qm_type'] = 'qmcart'
            atom_init['x'] = float(matchobj.group(8))
            atom_init['y'] = float(matchobj.group(9))
            atom_init['z'] = float(matchobj.group(10))

            mints_init['full_atoms'].append(atom_init)
            return ''

        #moldict['fragment_charges'] = self.fragment_charges
        #moldict['fragment_multiplicities'] = self.fragment_multiplicities

        #frag_reconstitute = []
        #start_atom = len(mints_init["full_atoms"])

        #for iln, line in enumerate(mol_str.split('\n')):
        #    line = line.strip()
        #    if iln == 0:
        #        line = re.sub(cgmp, process_fragment_cgmp, line)
        #    line = re.sub(atom_cart, process_atom_cart, line)
        #    #line = re.sub(atom_zmat1, process_atom_zmat1, line)
        #    #line = re.sub(atom_zmat2, process_atom_zmat2, line)
        #    #line = re.sub(atom_zmat3, process_atom_zmat3, line)
        #    #line = re.sub(atom_zmat4, process_atom_zmat4, line)
        #    if line:
        #        frag_reconstitute.append(line)

        #end_atom = len(mints_init["full_atoms"])
        #if end_atom > 0:
        #    #mints_init['fragments'].append([start_atom, end_atom + 1])
        #    mints_init['fragments'].append([start_atom, end_atom])
        #    mints_init['fragment_types'].append('Real')
        #    if len(mints_init['fragment_charges']) < len(mints_init['fragments']):
        #        mints_init['fragment_charges'].append(0)
        #        mints_init['fragment_multiplicities'].append(1)

        frag_reconstitute = []
        start_atom = len(mints_init["full_atoms"])

        for iln, line in enumerate(mol_str.split('\n')):
            line = re.sub(atom_cart, process_atom_cart, line.strip())
            #line = re.sub(atom_zmat1, process_atom_zmat1, line)
            #line = re.sub(atom_zmat2, process_atom_zmat2, line)
            #line = re.sub(atom_zmat3, process_atom_zmat3, line)
            #line = re.sub(atom_zmat4, process_atom_zmat4, line)
            if line:
                frag_reconstitute.append(line)

        end_atom = len(mints_init["full_atoms"])
        if (end_atom - start_atom) > 0:
            if re.sub(cgmp, process_fragment_cgmp, mol_str.split('\n')[0].strip()) != '':
                mints_init['fragment_charges'].append(None)
                mints_init['fragment_multiplicities'].append(None)
            #mints_init['fragments'].append([start_atom, end_atom + 1])
            mints_init['fragments'].append([start_atom, end_atom])
            mints_init['fragment_types'].append('Real')

        return '\n'.join(frag_reconstitute), mints_init

    reconstitute = []
    mints_init = {}
    mints_init['fragments'] = []
    mints_init['fragment_types'] = []
    mints_init['fragment_charges'] = []
    mints_init['fragment_multiplicities'] = []
    mints_init['full_atoms'] = []
    if seed:
        mints_init = seed

    # handle `--`-demarcated blocks
    for ifr, frag in enumerate(re.split(fragment_re, mol_str)):
        frag = frag.strip()
        if ifr == 0:
            frag = re.sub(cgmp, process_system_cgmp, frag)
        frag, mints_init = filter_fragment(frag)
        if frag:
            reconstitute.append(frag)

    return '\n--\n'.join(reconstitute), mints_init


def filter_libefp(mol_str, confer=None):
    """confer is read-only"""

    ENDL = r'[\t ,]*$'

    fragment = re.compile(r'^\s*--\s*$', re.MULTILINE)
    efpxyzabc = re.compile(
        r'\A' + r'efp' + SEP + r'(\w+)' +
        SEP + NUMBER + SEP + NUMBER + SEP + NUMBER +
        SEP + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL + r'\Z',
        re.IGNORECASE)
    efppoints = re.compile(
        r'\A' + r'efp' + r'\s+' + r'(\w+)' + ENDL +
        r'[\s,]*' + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL +
        r'[\s,]*' + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL +
        r'[\s,]*' + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL + r'\Z',
        re.IGNORECASE | re.MULTILINE)

    def process_efpxyzabc(matchobj):
        efp_init.append({'efp_type': 'xyzabc',
                         'fragment_file': matchobj.group(1),
                         'coordinates_hint': [float(matchobj.group(2)) * input_units_to_au,
                                              float(matchobj.group(3)) * input_units_to_au,
                                              float(matchobj.group(4)) * input_units_to_au,
                                              float(matchobj.group(5)),
                                              float(matchobj.group(6)),
                                              float(matchobj.group(7))]})
        return ''

    def process_efppoints(matchobj):
        efp_init.append({'efp_type': 'points',
                         'fragment_file': matchobj.group(1),
                         'coordinates_hint': [float(matchobj.group(2)) * input_units_to_au,
                                              float(matchobj.group(3)) * input_units_to_au,
                                              float(matchobj.group(4)) * input_units_to_au,
                                              float(matchobj.group(5)) * input_units_to_au,
                                              float(matchobj.group(6)) * input_units_to_au,
                                              float(matchobj.group(7)) * input_units_to_au,
                                              float(matchobj.group(8)) * input_units_to_au,
                                              float(matchobj.group(9)) * input_units_to_au,
                                              float(matchobj.group(10)) * input_units_to_au]})
        return ''

    reconstitute = []
    efp_init = []  # list of dicts, one for each efp fragment in mol_str

    # NOTE: applying libefp default of AU
    if confer and ('input_units_to_au' in confer):
        input_units_to_au = confer['input_units_to_au']
    else:
        input_units_to_au = 1.0


    # handle `--`-demarcated blocks
    for frag in re.split(fragment, mol_str):
        frag = re.sub(efpxyzabc, process_efpxyzabc, frag.strip())
        frag = re.sub(efppoints, process_efppoints, frag)
        if frag:
            reconstitute.append(frag)

    return '\n--\n'.join(reconstitute), efp_init


