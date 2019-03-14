#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

import copy
import pprint
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
from typing import Dict, List, Any

import numpy as np
import pydantic

from psi4 import core
from psi4.driver import p4util
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver import qcdb
from psi4.driver.task_base import BaseTask, SingleResult, unnp, plump_qcvar

# CONVENTIONS:
# n_ at the start of a variable name is short for "number of."
# _pi at the end of a variable name is short for "per irrep."
# h is the index of an irrep.

array_format = {"precision": 10}


def _displace_cart(mass, geom, salc_list, i_m, step_size):
    """Displace a geometry along the specified displacement SALCs.

    Parameters
    ----------
    mass : ndarray
        (nat, ) masses [u] of atoms of the molecule (const).
    geom : ndarray
        (nat, 3) reference geometry [a0] of the molecule (const).
    salc_list : :py:class:`~psi4.core.CdSalcList`
        A list of Cartesian displacement SALCs
    i_m : iterator of tuples
        An iterator containing tuples. Each tuple has the index of a salc in
        salc_list and the number of steps (positive or negative) to displace
        the salc at that index.
    step_size : float
        The size of a single "step," i.e., the stencil size.

    Returns
    -------
    disp_geom : ndarray
        (nat, 3) Displaced geometry.
    label : str
        Displacement label for the metadata dictionary.

    """
    label = []
    disp_geom = np.copy(geom)
    # This for loop and tuple unpacking is why the function can handle
    # an arbitrary number of SALCs.
    for salc_index, disp_steps in i_m:
        # * Python error if iterate through `salc_list`
        for i in range(len(salc_list[salc_index])):
            salc = salc_list[salc_index][i]
            disp_geom[salc.atom, salc.xyz] += disp_steps * step_size * salc.coef / np.sqrt(mass[salc.atom])
        label.append(f"{salc_index}: {disp_steps}")

    # salc_index is in descending order. We want the label in ascending order, so...
    # ...add the new label part from the left of the string, not the right.
    label = ', '.join(reversed(label))
    return disp_geom, label


def _initialize_findif(mol, freq_irrep_only, mode, stencil_size, step_size, initialize_string, verbose=0):
    """Perform initialization tasks needed by all primary functions.

    Parameters
    ----------
    mol : qcdb.molecule or :py:class:`~psi4.core.Molecule`
        The molecule to displace
    freq_irrep_only : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    mode : {"1_0", "2_0", "2_1"}
         The first number specifies the derivative level determined from
         displacements, and the second number is the level determined at.
    stencil_size : {3, 5}
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size : float
        [a0]

    initialize_string : function
         A function that returns the string to print to show the caller was entered.
         The string is both caller-specific and dependent on values determined
         in this function.
    verbose : int
         Set to 0 to silence extra print information, regardless of the print level.
         Used so the information is printed only during geometry generation, and not
         during the derivative computation as well.

    Returns
    -------
    data : dict
        Miscellaneous information required by callers.
    """

    core.print_out("\n         ----------------------------------------------------------\n")
    core.print_out("                                   FINDIF\n")
    core.print_out("                     R. A. King and Jonathon Misiewicz\n")
    core.print_out("         ---------------------------------------------------------\n\n")

    print_lvl = core.get_option("FINDIF", "PRINT")

    data = {"print_lvl": print_lvl, "stencil_size": stencil_size, "step_size": step_size}

    if print_lvl:
        core.print_out(initialize_string(data))

    # Get settings for CdSalcList, then get the CdSalcList.
    method_allowed_irreps = 0x1 if mode == "1_0" else 0xFF
    t_project = not core.get_global_option("EXTERN") and (not core.get_global_option("PERTURB_H"))
    # core.get_option returns an int, but CdSalcList expect a bool, so re-cast
    r_project = t_project and bool(core.get_option("FINDIF", "FD_PROJECT"))
    salc_list = core.CdSalcList(mol, method_allowed_irreps, t_project, r_project)

    n_atom = mol.natom()
    n_irrep = salc_list.nirrep()
    n_salc = salc_list.ncd()

    if print_lvl and verbose:
        core.print_out(f"    Number of atoms is {n_atom}.\n")
        if method_allowed_irreps != 0x1:
            core.print_out(f"    Number of irreps is {n_irrep}.\n")
        core.print_out("    Number of {!s}SALCs is {:d}.\n".format(
            "" if method_allowed_irreps != 0x1 else "symmetric ", n_salc))
        core.print_out("    Translations projected? {:d}. Rotations projected? {:d}.\n".format(t_project, r_project))

    # TODO: Replace with a generator from a stencil to a set of points.
    # Diagonal displacements differ between the totally symmetric irrep, compared to all others.
    # Off-diagonal displacements are the same for both.
    pts_dict = {
        3: {
            "sym_irr": ((-1, ), (1, )),
            "asym_irr": ((-1, ), ),
            "off": ((1, 1), (-1, -1))
        },
        5: {
            "sym_irr": ((-2, ), (-1, ), (1, ), (2, )),
            "asym_irr": ((-2, ), (-1, )),
            "off": ((-1, -2), (-2, -1), (-1, -1), (1, -1), (-1, 1), (1, 1), (2, 1), (1, 2))
        }
    }

    try:
        disps = pts_dict[stencil_size]
    except KeyError:
        raise ValidationError(f"FINDIF: Number of points ({stencil_size}) not among {pts_dict.keys()}!")

    # Convention: x_pi means x_per_irrep. The ith element is x for irrep i, with Cotton ordering.
    salc_indices_pi = [[] for h in range(n_irrep)]

    # Validate that we have an irrep matching the user-specified irrep, if any.
    try:
        salc_indices_pi[freq_irrep_only]
    except (TypeError, IndexError):
        if freq_irrep_only != -1:
            raise ValidationError(
                f"FINDIF: 0-indexed Irrep value ({freq_irrep_only}) not in valid range: <{len(salc_indices_pi)}.")

    # Populate salc_indices_pi for all irreps.
    # * Python error if iterate through `salc_list`
    for i in range(len(salc_list)):
        salc_indices_pi[salc_list[i].irrep_index()].append(i)

    # If the method allows more than one irrep, print how the irreps partition the SALCS.
    if print_lvl and method_allowed_irreps != 0x1 and verbose:
        core.print_out("    Index of SALCs per irrep:\n")
        for h in range(n_irrep):
            if print_lvl > 1 or freq_irrep_only in {h, -1}:
                tmp = (" {:d} " * len(salc_indices_pi[h])).format(*salc_indices_pi[h])
                core.print_out("     {:d} : ".format(h + 1) + tmp + "\n")
        core.print_out("    Number of SALCs per irrep:\n")
        for h in range(n_irrep):
            if print_lvl > 1 or freq_irrep_only in {h, -1}:
                core.print_out("     Irrep {:d}: {:d}\n".format(h + 1, len(salc_indices_pi[h])))

    # Now that we've printed the SALCs, clear any that are not of user-specified symmetry.
    if freq_irrep_only != -1:
        for h in range(n_irrep):
            if h != freq_irrep_only:
                salc_indices_pi[h].clear()

    n_disp_pi = []

    for irrep, indices in enumerate(salc_indices_pi):
        n_disp = len(indices) * len(disps["asym_irr" if irrep != 0 else "sym_irr"])
        if mode == "2_0":
            # Either len(indices) or len(indices)-1 is even, so dividing by two is safe.
            n_disp += len(indices) * (len(indices) - 1) // 2 * len(disps["off"])
        n_disp_pi.append(n_disp)

    # Let's print out the number of geometries, the displacement multiplicity, and the CdSALCs!
    if print_lvl and verbose:
        core.print_out("    Number of geometries (including reference) is {:d}.\n".format(sum(n_disp_pi) + 1))
        if method_allowed_irreps != 0x1:
            core.print_out("    Number of displacements per irrep:\n")
            for i, ndisp in enumerate(n_disp_pi, start=1):
                core.print_out(f"      Irrep {i}: {ndisp}\n")

    if print_lvl > 1 and verbose:
        for i in range(len(salc_list)):
            salc_list[i].print_out()

    data.update({
        "n_disp_pi": n_disp_pi,
        "n_irrep": n_irrep,
        "n_salc": n_salc,
        "n_atom": n_atom,
        "salc_list": salc_list,
        "salc_indices_pi": salc_indices_pi,
        "disps": disps,
        "project_translations": t_project,
        "project_rotations": r_project
    })

    return data


def _geom_generator(mol, freq_irrep_only, mode, stencil_size, step_size):
    """
    Generate geometries for the specified molecule and derivative levels.
    You probably want to instead use one of the convenience functions:
    gradient_from_energies_geometries, hessian_from_energies_geometries,
    hessian_from_gradients_geometries.

    Parameters
    ----------
    mol : qcdb.molecule or :py:class:`~psi4.core.Molecule`
        The molecule on which to perform a finite difference calculation.
    freq_irrep_only : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    mode : {"1_0", "2_0", "2_1"}
        The first number specifies the targeted derivative level. The
        second number is the compute derivative level. E.g., "2_0"
        is hessian from energies.
    stencil_size : {3, 5}
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size : float
        Displacement size [a0].

    Returns
    -------
    findifrec : dict
        Dictionary of finite difference data, specified below.
        The dictionary makes findifrec _extensible_. If you need a new field
        in the record, just add it.
        All fields should be present at all times, with two exceptions:
            1. Fields for computed quantities will not be available until
               after they are computed.
            2. Displacement specific overrides for globals will not be
               available unless the user specified the overrides.
               (Such overrides are not implemented at time of writing. An example
               is giving a displacement its own step dict.)

    step : dict
        A descriptor for the finite difference step.
        In future, this can be overriden by step fields for individual displacements.

        units : {'Bohr'}
            The units for the displacement. The code currently assumes "bohr," per MolSSI standards.
        size : float
            The step size for the displacement.

    stencil_size : {3, 5}
        Number of points to evaluate at for each displacement basis vector. Count
        includes the central reference point.

    displacement_space : {'CdSalc'}
        A string specifying the vector space in which displacements are performed.
        Currently, only CdSalc is supported.

    project_translations : bool
        Whether translations are to be projected out of the displacements.

    project_rotations : bool
        Whether rotations are to be projected out of the displacements.

    molecule : dict
        The reference molecule, in MolSSI schema. See
        https://molssi-qc-schema.readthedocs.io/en/latest/auto_topology.html

    displacements : dict
        A dictionary mapping labels specifying the displacement to data about
        the geometry. Labels are of the form "A: a, B: b" where A and B index the
        basis vector in displacement space and A < B, and a and b index the step
        magnitude. For instance, "0: 1, 1: -1" specifies displacing +1 in
        displacement vector 0 and -1 in displacement vector 1. "1: -1, 0: 1" is
        forbidden for breaking ordering. Generalizes to arbitrary numbers of
        simultaneous displacements in the obvious way.

        The possible geometry data is as follows:

        geometry: list of floats
            (3 * nat) The molecular geometry as a flat list in bohr. All coordinates
            are given for one atom before proceeding to the next atom.

        energy: int
            The last computed electronic energy at the geometry.

        gradient: list of floats
            (3 * nat) The last computed gradient of energy with respect to changes in
            geometry at the geometry, as a flat list. All coordinates are given for
            displacing one atom before proceeding to the next atom.

    reference : dict
         A geometry data dict, as described above, for the reference geometry.
    """

    msg_dict = {
        "1_0":
        "energies to determine gradients",
        "2_1":
        "gradients to determine vibrational frequencies and \n"
        "  normal modes. Resulting frequencies are only valid at stationary points",
        "2_0":
        "gradients to determine vibrational frequencies and \n"
        "  normal modes. Resulting frequencies are only valid at stationary points"
    }

    try:
        print_msg = msg_dict[mode]
    except KeyError:
        raise ValidationError("FINDIF: Mode {} not recognized.".format(mode))

    def init_string(data):
        return f"""  Using finite-differences of {print_msg}.
    Generating geometries for use with {data["stencil_size"]}-point formula.
    Displacement size will be {data["step_size"]:6.2e}.\n"""

    # Genuine support for qcdb molecules would be nice. But that requires qcdb CdSalc tech.
    # Until then, silently swap the qcdb molecule out for a psi4.core.molecule.
    if isinstance(mol, qcdb.Molecule):
        mol = core.Molecule.from_dict(mol.to_dict())

    data = _initialize_findif(mol, freq_irrep_only, mode, stencil_size, step_size, init_string, 1)

    # We can finally start generating displacements.
    ref_geom = np.array(mol.geometry())

    # Now we generate the metadata...
    findifrec = {
        "step": {
            "units": "bohr",
            "size": data["step_size"]
        },
        "stencil_size": data["stencil_size"],
        "displacement_space": "CdSALC",
        "project_translations": data["project_translations"],
        "project_rotations": data["project_rotations"],
        "molecule": mol.to_schema(dtype=1, units='Bohr'),
        "displacements": {},
        "reference": {}
    }

    def append_geoms(indices, steps):
        """Given a list of indices and a list of steps to displace each, append the corresponding geometry to the list."""

        # Next, to make this salc/magnitude composite.
        disp_geom, label = _displace_cart(findifrec['molecule']['molecule']['masses'], ref_geom, data["salc_list"],
                                          zip(indices, steps), data["step_size"])
        if data["print_lvl"] > 2:
            core.print_out("\nDisplacement '{}'\n{}\n".format(label, np.array_str(disp_geom, **array_format)))
        findifrec["displacements"][label] = {"geometry": disp_geom}

    for h in range(data["n_irrep"]):
        active_indices = data["salc_indices_pi"][h]

        for index in active_indices:
            # Displace along the diagonal.
            # Remember that the totally symmetric irrep has special displacements.
            for val in data["disps"]["sym_irr" if h == 0 else "asym_irr"]:
                append_geoms((index, ), val)

        # Hessian from energies? We have off-diagonal displacements to worry about.
        if mode == "2_0":
            # i indexes SALC indices of the current irrep.
            for i, index in enumerate(active_indices):
                for index2 in active_indices[:i]:
                    for val in data["disps"]["off"]:
                        append_geoms((index, index2), val)

    if data["print_lvl"] > 2:
        core.print_out("\nReference\n{}\n".format(np.array_str(ref_geom, **array_format)))
    findifrec["reference"]["geometry"] = ref_geom

    if data["print_lvl"] > 1:
        core.print_out("\n-------------------------------------------------------------\n")

    return findifrec


def assemble_gradient_from_energies(findifrec):
    """Compute the gradient by finite difference of energies.

    Parameters
    ----------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.

    Returns
    -------
    gradient : ndarray
        (nat, 3) Cartesian gradient [Eh/a0].
    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    def init_string(data):
        return f"""  Computing gradient from energies.
    Using {findifrec["stencil_size"]}-point formula.
    Energy without displacement: {findifrec["reference"]["energy"]:15.10f}
    Check energies below for precision!
    Forces are for mass-weighted, symmetry-adapted cartesians [a0]."""

    data = _initialize_findif(mol, -1, "1_0", findifrec['stencil_size'], findifrec['step']['size'], init_string)
    salc_indices = data["salc_indices_pi"][0]

    # Extract the energies, and turn then into an ndarray for easy manipulating
    # E(i, j) := Energy on displacing the ith SALC we care about in the jth step
    # Steps are ordered, for example, -2, -1, 1, 2
    max_disp = (findifrec["stencil_size"] - 1) // 2  # The numerator had better be divisible by two.
    e_per_salc = 2 * max_disp
    E = np.zeros((len(salc_indices), e_per_salc))

    for i, salc_index in enumerate(salc_indices):
        for j in range(1, max_disp + 1):
            E[i, max_disp - j] = findifrec["displacements"][f"{salc_index}: {-j}"]["energy"]
            E[i, max_disp + j - 1] = findifrec["displacements"][f"{salc_index}: {j}"]["energy"]

    # Perform the finite difference.
    if findifrec["stencil_size"] == 3:
        g_q = (E[:, 1] - E[:, 0]) / (2.0 * findifrec["step"]["size"])
    elif findifrec["stencil_size"] == 5:
        g_q = (E[:, 0] - 8.0 * E[:, 1] + 8.0 * E[:, 2] - E[:, 3]) / (12.0 * findifrec["step"]["size"])
    else:  # This error SHOULD have already been caught, but just in case...
        raise ValidationError("FINDIF: {} is an invalid number of points.".format(findifrec["stencil_size"]))
    g_q = np.asarray(g_q)

    if data["print_lvl"]:
        energy_string = ""
        for i in range(1, max_disp + 1):
            energy_string = f"Energy(-{i})        " + energy_string + f"Energy(+{i})        "
        core.print_out("\n     Coord      " + energy_string + "    Force\n")
        for salc in range(data["n_salc"]):
            print_str = "    {:5d}" + " {:17.10f}" * (e_per_salc) + " {force:17.10f}" + "\n"
            energies = E[salc]
            core.print_out(print_str.format(salc, force=g_q[salc], *energies))
        core.print_out("\n")

    # Transform the gradient from mass-weighted SALCs to non-mass-weighted Cartesians
    B = data["salc_list"].matrix()
    g_cart = np.dot(g_q, B)
    g_cart = g_cart.reshape(data["n_atom"], 3)
    massweighter = np.array([mol.mass(a) for a in range(data["n_atom"])])**(0.5)
    g_cart = (g_cart.T * massweighter).T

    if data["print_lvl"]:
        core.print_out("\n-------------------------------------------------------------\n")

    return g_cart


def _process_hessian_symmetry_block(H_block, B_block, massweighter, irrep, print_lvl):
    """Perform post-construction processing for a symmetry block of the Hessian.
       Statements need to be printed, and the Hessian must be made orthogonal.

    Parameters
    ---------
    H_block : ndarray
        A block of the Hessian for an irrep, in mass-weighted salcs.
        Dimensions # cdsalcs by # cdsalcs.
    B_block : ndarray
        A block of the B matrix for an irrep, which transforms CdSalcs to Cartesians.
        Dimensions # cdsalcs by # cartesians.
    massweighter : ndarray
        The mass associated with each atomic coordinate.
        Dimension # cartesians. Due to x, y, z, values appear in groups of three.
    irrep : str
        A string identifying the irrep H_block and B_block are of.
    print_lvl : int
        The level of printing information requested by the user.

    Returns
    -------
    H_block : ndarray
        H_block, but made into an orthogonal array.
    """

    # Symmetrize our Hessian block.
    # The symmetric structure is lost due to errors in the computation
    H_block = (H_block + H_block.T) / 2.0

    if print_lvl >= 3:
        core.print_out("\n    Force Constants for irrep {} in mass-weighted, ".format(irrep))
        core.print_out("symmetry-adapted cartesian coordinates.\n")
        core.print_out("\n{}\n".format(np.array_str(H_block, **array_format)))

    evals, evects = np.linalg.eigh(H_block)
    # Get our eigenvalues and eigenvectors in descending order.
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evects = evects[:, idx]

    normal_irr = np.dot((B_block * massweighter).T, evects)

    if print_lvl >= 2:
        core.print_out("\n    Normal coordinates (non-mass-weighted) for irrep {}:\n".format(irrep))
        core.print_out("\n{}\n".format(np.array_str(normal_irr, **array_format)))

    return H_block


def _process_hessian(H_blocks, B_blocks, massweighter, print_lvl):
    """Perform post-construction processing for the Hessian.
       Statements need to be printed, and the Hessian must be transformed.

    Parameters
    ----------
    H_blocks : list of ndarray
        A list of blocks of the Hessian per irrep, in mass-weighted salcs.
        Each is dimension # cdsalcs-in-irrep by # cdsalcs-in-irrep.
    B_blocks : list of ndarray
        A block of the B matrix per irrep, which transforms CdSalcs to Cartesians.
        Each is dimensions # cdsalcs-in-irrep by # cartesians.
    massweighter : ndarray
        The mass associated with each atomic coordinate.
        Dimension 3 * natom. Due to x, y, z, values appear in groups of three.
    print_lvl : int
        The level of printing information requested by the user.

    Returns
    -------
    Hx : ndarray
        The Hessian in non-mass weighted cartesians.
    """

    # We have the Hessian in each irrep! The final task is to perform coordinate transforms.
    H = p4util.block_diagonal_array(*H_blocks)
    B = np.vstack(B_blocks)

    if print_lvl >= 3:
        core.print_out("\n    Force constant matrix for all computed irreps in mass-weighted SALCS.\n")
        core.print_out("\n{}\n".format(np.array_str(H, **array_format)))

    # Transform the massweighted Hessian from the CdSalc basis to Cartesians.
    # The Hessian is the matrix not of a linear transformation, but of a (symmetric) bilinear form
    # As such, the change of basis is formula A' = Xt A X, no inverses!
    # More conceptually, it's A'_kl = A_ij X_ik X_jl; Each index transforms linearly.
    Hx = np.dot(np.dot(B.T, H), B)
    if print_lvl >= 3:
        core.print_out("\n    Force constants in mass-weighted Cartesian coordinates.\n")
        core.print_out("\n{}\n".format(np.array_str(Hx, **array_format)))

    # Un-massweight the Hessian.
    Hx = np.transpose(Hx / massweighter) / massweighter

    if print_lvl >= 3:
        core.print_out("\n    Force constants in Cartesian coordinates.\n")
        core.print_out("\n{}\n".format(np.array_str(Hx, **array_format)))

    if print_lvl:
        core.print_out("\n-------------------------------------------------------------\n")

    return Hx


def assemble_hessian_from_gradients(findifrec, freq_irrep_only):
    """Compute the Hessian by finite difference of gradients.

    Parameters
    ----------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.
    freq_irrep_only : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.

    Returns
    -------
    hessian : ndarray
        (3 * nat, 3 * nat) Cartesian Hessian [Eh/a0^2]
    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    displacements = findifrec["displacements"]

    def init_string(data):
        return ("  Computing second-derivative from gradients using projected, \n"
                "  symmetry-adapted, cartesian coordinates.\n\n"
                "  {:d} gradients passed in, including the reference geometry.\n".format(len(displacements) + 1))

    data = _initialize_findif(mol, freq_irrep_only, "2_1", findifrec['stencil_size'], findifrec['step']['size'],
                              init_string)

    # For non-totally symmetric CdSALCs, a symmetry operation can convert + and - displacements.
    # Good News: By taking advantage of that, we (potentially) ran less computations.
    # Bad News: We need to find the - displacements from the + computations now.
    # The next ~80 lines of code are dedicated to that task.
    if data["print_lvl"]:
        core.print_out("  Generating complete list of displacements from unique ones.\n\n")

    pg = mol.point_group()
    ct = pg.char_table()
    order = pg.order()

    # Determine what atoms map to what other atoms under the point group operations.
    # The py-side compute_atom_map will work whether mol is a Py-side or C-side object.
    atom_map = qcdb.compute_atom_map(mol)
    if data["print_lvl"] >= 3:
        core.print_out("    The atom map:\n")
        for atom, sym_image_list in enumerate(atom_map):
            core.print_out(f"     {atom + 1:d} : ")
            for image_atom in sym_image_list:
                core.print_out(f"{image_atom + 1:4d}")
            core.print_out("\n")
        core.print_out("\n")

    # A list of lists of gradients, per irrep
    gradients_pi = [[]]
    # Extract and print the symmetric gradients. These need no additional processing.
    max_disp = (findifrec["stencil_size"] - 1) // 2  # The numerator had better be divisible by two.
    for i in data["salc_indices_pi"][0]:
        for n in range(-max_disp, 0):
            grad_raw = displacements[f"{i}: {n}"]["gradient"]
            gradients_pi[0].append(np.reshape(grad_raw, (-1, 3)))
        for n in range(1, max_disp + 1):
            grad_raw = displacements[f"{i}: {n}"]["gradient"]
            gradients_pi[0].append(np.reshape(grad_raw, (-1, 3)))

    if data["print_lvl"] >= 3:
        core.print_out("    Symmetric gradients\n")
        for gradient in gradients_pi[0]:
            core.print_out("\n{}\n".format(np.array_str(gradient, **array_format)))

    # Asymmetric gradient. There's always SOME operation that transforms a positive
    # into a negative displacement.By doing extra things here, we can find the
    # gradients at the positive displacements.
    for h in range(1, data["n_irrep"]):

        # If there are no CdSALCs in this irrep, let's skip it.
        if not data["n_disp_pi"][h]:
            gradients_pi.append([])
            continue

        gamma = ct.gamma(h)
        if data["print_lvl"] >= 3:
            core.print_out(f"Characters for irrep {h}\n")
            for group_op in range(order):
                core.print_out(" {:5.1f}".format(gamma.character(group_op)))
            core.print_out("\n")

        # Find the group operation that converts + to - displacements.
        for group_op in range(order):
            if gamma.character(group_op) == -1:
                break
        else:
            raise ValidationError("A symmetric gradient passed for a non-symmetric one.")
        if data["print_lvl"]:
            core.print_out("    Operation {} takes plus displacements of irrep {} to minus ones.\n".format(
                group_op + 1, gamma.symbol()))

        sym_op = np.array(ct.symm_operation(group_op).matrix())
        gradients = []

        def recursive_gradients(i, n):
            """Populate gradients, with step -n, -n+1, ... -1, 1, ... n. Positive displacements are computed."""

            grad_raw = displacements[f"{i}: {-n}"]["gradient"]
            gradients.append(np.reshape(grad_raw, (-1, 3)))
            new_grad = np.zeros((data["n_atom"], 3))
            for atom, image in enumerate(atom_map):
                atom2 = image[group_op]
                new_grad[atom2] = np.einsum("xy,y->x", sym_op, gradients[-1][atom])
            if n > 1:
                recursive_gradients(i, n - 1)
            gradients.append(new_grad)

        for i in data["salc_indices_pi"][h]:
            recursive_gradients(i, max_disp)
        gradients_pi.append(gradients)

    # Massweight all gradients.
    # Remember, the atom currently corresponds to our 0 axis, hence these transpose tricks.
    massweighter = np.asarray([mol.mass(a) for a in range(data["n_atom"])])**(-0.5)
    gradients_pi = [[(grad.T * massweighter).T for grad in gradients] for gradients in gradients_pi]

    if data["print_lvl"] >= 3:
        core.print_out("    All mass-weighted gradients\n")
        for gradients in gradients_pi:
            for grad in gradients:
                core.print_out("\n{}\n".format(np.array_str(grad, **array_format)))

    # We have all our gradients generated now!
    # Next, time to get our Hessian.

    H_pi = []
    B_pi = []
    irrep_lbls = mol.irrep_labels()
    massweighter = np.repeat(massweighter, 3)

    for h in range(data["n_irrep"]):
        n_disp = data["n_disp_pi"][h]
        Nindices = len(data["salc_indices_pi"][h])
        gradients = gradients_pi[h]

        if not Nindices:
            continue

        # Flatten each gradient, and turn it into a COLUMN of the matrix.
        gradient_matrix = np.array([grad.flatten() for grad in gradients]).T
        # Transform disps from Cartesian to CdSalc coordinates.
        # For future convenience, we transpose.
        # Rows are gradients and columns are coordinates with respect to a particular CdSALC.
        B_pi.append(data["salc_list"].matrix_irrep(h))
        grads_adapted = np.dot(B_pi[-1], gradient_matrix).T

        if data["print_lvl"] >= 3:
            core.print_out("Gradients in B-matrix coordinates\n")
            for disp in range(n_disp):
                core.print_out(f" disp {disp}: ")
                for salc in grads_adapted[disp]:
                    core.print_out(f"{salc:15.10f}")
                core.print_out("\n")

        H_pi.append(np.empty([Nindices, Nindices]))

        if findifrec["stencil_size"] == 3:
            H_pi[-1] = (grads_adapted[1::2] - grads_adapted[::2]) / (2.0 * findifrec["step"]["size"])
        elif findifrec["stencil_size"] == 5:
            H_pi[-1] = (grads_adapted[::4] - 8 * grads_adapted[1::4] + 8 * grads_adapted[2::4] -
                        grads_adapted[3::4]) / (12.0 * findifrec["step"]["size"])

        H_pi[-1] = _process_hessian_symmetry_block(H_pi[-1], B_pi[-1], massweighter, irrep_lbls[h], data["print_lvl"])

    # All blocks of the Hessian are now constructed!
    return _process_hessian(H_pi, B_pi, massweighter, data["print_lvl"])


def assemble_hessian_from_energies(findifrec, freq_irrep_only):
    """Compute the Hessian by finite difference of energies.

    Parameters
    ----------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.
    freq_irrep_only : int
        The 0-indexed Cotton ordered irrep to get frequencies for. Choose -1 for all irreps.

    Returns
    -------
    hessian : ndarray
        (3 * nat, 3 * nat) Cartesian Hessian [Eh/a0^2].
    """

    # This *must* be a Psi molecule at present - CdSalcList generation panics otherwise
    mol = core.Molecule.from_schema(findifrec["molecule"], nonphysical=True, verbose=0)

    displacements = findifrec["displacements"]
    ref_energy = findifrec["reference"]["energy"]

    def init_string(data):
        max_label_len = str(max([len(label) for label in displacements]))
        out_str = ""
        for label, disp_data in displacements.items():
            out_str += ("    {:" + max_label_len + "s} : {:20.10f}\n").format(label, disp_data["energy"])
        return ("  Computing second-derivative from energies using projected, \n"
                "  symmetry-adapted, cartesian coordinates.\n\n"
                "  {:d} energies passed in, including the reference geometry.\n"
                "    Using {:d}-point formula.\n"
                "    Energy without displacement: {:15.10f}\n"
                "    Check energies below for precision!\n{}".format(
                    len(displacements) + 1, findifrec["stencil_size"], ref_energy, out_str))

    data = _initialize_findif(mol, freq_irrep_only, "2_0", findifrec['stencil_size'], findifrec['step']['size'],
                              init_string)

    massweighter = np.repeat([mol.mass(a) for a in range(data["n_atom"])], 3)**(-0.5)
    B_pi = []
    H_pi = []
    irrep_lbls = mol.irrep_labels()
    max_disp = (findifrec["stencil_size"] - 1) // 2
    e_per_diag = 2 * max_disp

    # Unlike in the gradient case, we have no symmetry transformations to worry about.
    # We get to the task directly: assembling the force constants in each irrep block.
    for h in range(data["n_irrep"]):
        salc_indices = data["salc_indices_pi"][h]
        if not salc_indices: continue

        n_salcs = len(salc_indices)
        E = np.zeros((len(salc_indices), e_per_diag))

        # Step One: Diagonals
        # For asymmetric irreps, the energy at a + disp is the same as at a - disp
        # Just reuse the - disp energy for the + disp energy

        for i, salc_index in enumerate(salc_indices):
            for j in range(1, max_disp + 1):
                E[i, max_disp - j] = displacements[f"{salc_index}: {-j}"]["energy"]
                k = -j if h else j  # Because of the +- displacement trick
                E[i, max_disp + j - 1] = displacements[f"{salc_index}: {k}"]["energy"]
        # Now determine all diagonal force constants for this irrep.
        if findifrec["stencil_size"] == 3:
            diag_fcs = E[:, 0] + E[:, 1]
            diag_fcs -= 2 * ref_energy
            diag_fcs /= (findifrec["step"]["size"]**2)
        elif findifrec["stencil_size"] == 5:
            diag_fcs = -E[:, 0] + 16 * E[:, 1] + 16 * E[:, 2] - E[:, 3]
            diag_fcs -= 30 * ref_energy
            diag_fcs /= (12 * findifrec["step"]["size"]**2)
        H_irr = np.diag(diag_fcs)

        # TODO: It's a bit ugly to use the salc indices to grab the off-diagonals but the indices
        # within the irrep to grab the diagonals. Is there a better way to do this?

        # Step Two: Off-diagonals
        # We need off-diagonal energies, diagonal energies, AND the reference energy
        # Grabbing off-diagonal energies is a pain, so once we know our SALCs...
        # ...define offdiag_en to do that for us.
        for i, salc in enumerate(salc_indices):
            for j, salc2 in enumerate(salc_indices[:i]):
                offdiag_en = lambda index: displacements["{l}: {}, {k}: {}".format(k=salc, l=salc2, *data["disps"]["off"][index])]["energy"]
                if findifrec["stencil_size"] == 3:
                    fc = (+offdiag_en(0) + offdiag_en(1) + 2 * ref_energy - E[i][0] - E[i][1] - E[j][0] - E[j][1]) / (
                        2 * findifrec["step"]["size"]**2)
                elif findifrec["stencil_size"] == 5:
                    fc = (-offdiag_en(0) - offdiag_en(1) + 9 * offdiag_en(2) - offdiag_en(3) - offdiag_en(4) +
                          9 * offdiag_en(5) - offdiag_en(6) - offdiag_en(7) + E[i][0] - 7 * E[i][1] - 7 * E[i][2] +
                          E[i][3] + E[j][0] - 7 * E[j][1] - 7 * E[j][2] + E[j][3] + 12 * ref_energy) / (
                              12 * findifrec["step"]["size"]**2)
                H_irr[i, j] = fc
                H_irr[j, i] = fc

        B_pi.append(data["salc_list"].matrix_irrep(h))
        H_pi.append(_process_hessian_symmetry_block(H_irr, B_pi[-1], massweighter, irrep_lbls[h], data["print_lvl"]))

    # All blocks of the Hessian are now constructed!
    return _process_hessian(H_pi, B_pi, massweighter, data["print_lvl"])


def gradient_from_energies_geometries(molecule, stencil_size=3, step_size=0.005):
    """
    Generate geometries for a gradient by finite difference of energies.

    Parameters
    ----------
    molecule : qcdb.molecule or :py:class:`~psi4.core.Molecule`
        The molecule to compute the gradient of.
    stencil_size : {3, 5}, optional
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size : float, optional
        Displacement size [a0].

    Returns
    -------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.

    Notes
    -----
    Only symmetric displacements are necessary, so user specification of
    symmetry is disabled.
    """
    return _geom_generator(molecule, -1, "1_0", stencil_size, step_size)


def hessian_from_gradients_geometries(molecule, irrep, stencil_size=3, step_size=0.005):
    """
    Generate geometries for a hessian by finite difference of gradients.

    Parameters
    ----------
    molecule : qcdb.molecule or :py:class:`~psi4.core.Molecule`
        The molecule to compute the frequencies of.
    irrep : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    stencil_size : {3, 5}, optional
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size : float, optional

    Returns
    -------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.
    """
    return _geom_generator(molecule, irrep, "2_1", stencil_size, step_size)


def hessian_from_energies_geometries(molecule, irrep, stencil_size=3, step_size=0.005):
    """
    Generate geometries for a hessian by finite difference of energies.

    Parameters
    ----------
    molecule : qcdb.molecule or :py:class:`~psi4.core.Molecule`
        The molecule to compute the frequencies of.
    irrep : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    stencil_size : {3, 5}, optional
        Number of points to evaluate for each displacement basis vector inclusive of central reference geometry.
    step_size : float, optional

    Returns
    -------
    findifrec : dict
        Dictionary of finite difference data, specified in _geom_generator docstring.
    """
    return _geom_generator(molecule, irrep, "2_0", stencil_size, step_size)


class FinDifComputer(BaseTask):

    molecule: Any
    driver: str
    metameta: Dict[str, Any] = {}
    task_list: Dict[str, SingleResult] = {}
    findifrec: Dict[str, Any] = {}

    @pydantic.validator('driver')
    def set_driver(cls, driver):
        egh = ['energy', 'gradient', 'hessian']
        if driver not in egh:
            raise ValidationError(f"""Wrapper is unhappy to be calling function ({driver}) not among {egh}.""")

        return driver

    @pydantic.validator('molecule')
    def set_molecule(cls, mol):
        mol.update_geometry()
        mol.fix_com(True)
        mol.fix_orientation(True)
        return mol

    def __init__(self, **data):
        BaseTask.__init__(self, **data)

        print('FINDIFREC CLASS INIT DATA')
        pp.pprint(data)
        if 'dft_functional' in data:
            data['keywords']['dft_functional'] = data.pop('dft_functional')
        data['keywords']['PARENT_SYMMETRY'] = self.molecule.point_group().full_name()

        self.metameta['mode'] = str(data['findif_mode'][0]) + '_' + str(data['findif_mode'][1])
        self.metameta['irrep'] = data.pop('irrep', -1)

        if self.metameta['mode'] == '1_0':
            self.metameta['proxy_driver'] = 'energy'
            self.findifrec = gradient_from_energies_geometries(self.molecule, data.pop('findif_stencil_size'),
                                                               data.pop('findif_step_size'))

        elif self.metameta['mode'] == '2_1':
            self.metameta['proxy_driver'] = 'gradient'
            self.findifrec = hessian_from_gradients_geometries(self.molecule, self.metameta['irrep'],
                                                               data.pop('findif_stencil_size'),
                                                               data.pop('findif_step_size'))

        elif self.metameta['mode'] == '2_0':
            self.metameta['proxy_driver'] = 'energy'
            self.findifrec = hessian_from_energies_geometries(self.molecule, self.metameta['irrep'],
                                                              data.pop('findif_stencil_size'),
                                                              data.pop('findif_step_size'))

        print('FINDIFREC CLASS META DATA')
        pp.pprint(self.metameta)
        print('FINDIFREC CLASS')
        pp.pprint(self.findifrec)

        ndisp = len(self.findifrec["displacements"]) + 1
        print(f""" {ndisp} displacements needed ...""", end='')

        # var_dict = core.variables()
        reftask = SingleResult(
            **{
                "molecule": self.molecule,
                "driver": self.metameta['proxy_driver'],
                "method": data["method"],
                "basis": data["basis"],
                "keywords": data["keywords"] or {},
            })
        self.task_list["reference"] = reftask

        parent_group = self.molecule.point_group()
        for label, displacement in self.findifrec["displacements"].items():
            clone = self.molecule.clone()
            clone.reinterpret_coordentry(False)
            #clone.fix_orientation(True)

            # Load in displacement into the active molecule
            clone.set_geometry(core.Matrix.from_array(displacement["geometry"]))

            # If the user insists on symmetry, weaken it if some is lost when displacing.
            # or 'fix_symmetry' in self.findifrec.molecule
            print('SYMM', clone.schoenflies_symbol())
            if self.molecule.symmetry_from_input():
                disp_group = clone.find_highest_point_group()
                new_bits = parent_group.bits() & disp_group.bits()
                new_symm_string = qcdb.PointGroup.bits_to_full_name(new_bits)
                clone.reset_point_group(new_symm_string)

            self.task_list[label] = SingleResult(
                **{
                    "molecule": clone,
                    "driver": self.metameta['proxy_driver'],
                    "method": data["method"],
                    "basis": data["basis"],
                    "keywords": data["keywords"] or {},
                })

#        for n, displacement in enumerate(findif_meta_dict["displacements"].values(), start=2):
#            _process_displacement(energy, lowername, molecule, displacement, n, ndisp, write_orbitals=False, **kwargs)

    def build_tasks(self, obj, **kwargs):
        # permanently a dummy function
        pass

    def plan(self):
        return [x.plan() for x in self.task_list.values()]

    def compute(self):
        with p4util.hold_options_state():
            # gof = core.get_output_file()
            # core.close_outfile()

            for x in self.task_list.values():
                x.compute()

            # core.set_output_file(gof, True)

    def _prepare_results(self):
        results_list = {k: v.get_results() for k, v in self.task_list.items()}

        #for i, x in self.task_list.items():
        #    print('\nTASK', i)
        #    pp.pprint(x)
        #for i, x in results_list.items():
        #    print('\nRESULT', i)
        #    pp.pprint(x)

        # load SingleResult results into findifrec[reference]
        reference = self.findifrec["reference"]
        task = results_list["reference"]
        response = task['return_result']

        if task['driver'] == 'energy':
            reference['energy'] = response

        elif task['driver'] == 'gradient':
            reference['gradient'] = plump_qcvar(response, 'gradient')
            reference['energy'] = task['psi4:qcvars']['CURRENT ENERGY']

        elif task['driver'] == 'hessian':
            reference['hessian'] = plump_qcvar(response, 'hessian')
            reference['energy'] = task['psi4:qcvars']['CURRENT ENERGY']
            if 'CURRENT GRADIENT' in task['psi4:qcvars']:
                reference['gradient'] = plump_qcvar(task['psi4:qcvars']['CURRENT GRADIENT'], 'gradient')

        # load SingleResult results into findifrec[displacements]
        for label, displacement in self.findifrec["displacements"].items():
            task = results_list[label]
            response = task['return_result']

            if task['driver'] == 'energy':
                displacement['energy'] = response

            elif task['driver'] == 'gradient':
                displacement['gradient'] = plump_qcvar(response, 'gradient')
                displacement['energy'] = task['psi4:qcvars']['CURRENT ENERGY']

            elif task['driver'] == 'hessian':
                displacement['hessian'] = plump_qcvar(response, 'hessian')
                displacement['energy'] = task['psi4:qcvars']['CURRENT ENERGY']
                if 'CURRENT GRADIENT' in task['psi4:qcvars']:
                    displacement['gradient'] = plump_qcvar(task['psi4:qcvars']['CURRENT GRADIENT'], 'gradient')

        # apply finite difference formulas and load derivatives into findifrec[reference]
        if self.metameta['mode'] == '1_0':
            G0 = assemble_gradient_from_energies(self.findifrec)
            self.findifrec["reference"][self.driver] = G0

        elif self.metameta['mode'] == '2_1':
            H0 = assemble_hessian_from_gradients(self.findifrec, self.metameta['irrep'])
            self.findifrec["reference"][self.driver] = H0

        elif self.metameta['mode'] == '2_0':
            try:
                G0 = assemble_gradient_from_energies(self.findifrec)
            except KeyError:
                core.print_out("Unable to construct reference gradient from Hessian displacements.")
                # TODO: this happens properly when a subset irrep doesn't
                #  have the right displacements for grad. For both this case
                #  and distributed computing are-we-there-yet? queries,
                #  should have a probe as to whether all the
                #  findif[displacement] labels are present and whether
                #  all the findif[displacement][energy-or-gradient] values
                #  are ready. Not sure what type of error/query is best,
                #  so deferring for now. Also, possibly need to check if
                #  step size matches before using values from one findifrec
                #  to construct another quantity.
            else:
                self.findifrec["reference"]["gradient"] = G0

            H0 = assemble_hessian_from_energies(self.findifrec, self.metameta['irrep'])
            self.findifrec["reference"][self.driver] = H0


#        if core.get_option('FINDIF', 'GRADIENT_WRITE'):
#            filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".grad"
#            qcdb.gradparse.to_string(np.asarray(wfn.gradient()), filename, dtype='GRD', mol=molecule, energy=wfn.energy())

        print('\nFINDIF_RESULTS POST-LOAD')
        pp.pprint(self.findifrec)

    def get_results(self):
        assembled_results = self._prepare_results()

        # load QCVariables
        qcvars = self.task_list['reference'].get_results()['psi4:qcvars']

        #qcvars['CURRENT REFERENCE ENERGY'] = self.grand_need[0]['d_energy']
        #qcvars['CURRENT CORRELATION ENERGY'] = assembled_results['energy'] - self.grand_need[0]['d_energy']
        qcvars['FINDIF NUMBER'] = len(self.task_list)
        qcvars['NUCLEAR REPULSION ENERGY'] = self.molecule.nuclear_repulsion_energy()
        qcvars['CURRENT ENERGY'] = self.findifrec['reference']['energy']

        G0 = self.findifrec['reference'].get('gradient')
        if G0 is not None:
            qcvars['CURRENT GRADIENT'] = G0

        H0 = self.findifrec['reference'].get('hessian')
        if H0 is not None:
            qcvars['CURRENT HESSIAN'] = H0

        findifjob = {
            'findif_record': copy.deepcopy(self.findifrec),
            'driver': self.driver,
            # 'keywords':
            'model': {
                # 'basis':
                # 'method':
            },
            'molecule': self.molecule.to_schema(dtype=1)['molecule'],
            # 'properties':
            'provenance': p4util.provenance_stamp(__name__),
            'psi4:qcvars': qcvars,
            'return_result': self.findifrec['reference'][self.driver],
            'schema_name': 'qc_schema_output',
            'schema_version': 1,
            # 'success': True,
        }

        findifjob = unnp(findifjob, flat=True)
        print('\nFINDIF QCSchema:')
        pp.pprint(findifjob)
        return findifjob

    def get_psi_results(self, return_wfn=False):

        findifjob = self.get_results()

        ret_ptype = plump_qcvar(findifjob['return_result'], shape_clue=findifjob['driver'], ret='psi4')
        wfn = _findif_schema_to_wfn(findifjob)

        if return_wfn:
            return (ret_ptype, wfn)
        else:
            return ret_ptype


def _findif_schema_to_wfn(findifjob):
    """Helper function to keep Wavefunction dependent on FinDif-flavored QCSchemus."""

    # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
    mol = core.Molecule.from_schema(findifjob)
    basis = core.BasisSet.build(mol, "ORBITAL", 'def2-svp')
    wfn = core.Wavefunction(mol, basis)

#    wfn.set_energy(findifjob['psi4:qcvars'].get('CURRENT ENERGY'))  # catches Wfn.energy_
    for qcv, val in findifjob['psi4:qcvars'].items():
        for obj in [core, wfn]:
            obj.set_variable(qcv, plump_qcvar(val, qcv))

#    flat_grad = findifjob['psi4:qcvars'].get('CURRENT GRADIENT')
#    if flat_grad is not None:
#        finalgradient = plump_qcvar(flat_grad, 'gradient', ret='psi4')
#        wfn.set_gradient(finalgradient)
#
#        if finalgradient.rows(0) < 20:
#            core.print_out('CURRENT GRADIENT')
#            finalgradient.print_out()
#
#    flat_hess = findifjob['psi4:qcvars'].get('CURRENT HESSIAN')
#    if flat_hess is not None:
#        finalhessian = plump_qcvar(flat_hess, 'hessian', ret='psi4')
#        wfn.set_hessian(finalhessian)
#
#        if finalhessian.rows(0) < 20:
#            core.print_out('CURRENT HESSIAN')
#            finalhessian.print_out()

    return wfn
