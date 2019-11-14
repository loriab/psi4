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

"""Module with functions that call upon those in modules
:py:mod:`proc`, :py:mod:`driver`, and :py:mod:`wrappers`.

Place in this file quickly defined procedures such as
   - aliases for complex methods
   - simple modifications to existing methods

"""
import os
import re
import math
import warnings


# Python procedures like these can be run directly from the input file or integrated
# with the energy(), etc. routines by means of lines like those at the end
# of this file.


def fake_file11(wfn, filename='fake_file11.dat', **kwargs):
    r"""Function to print a file *filename* of the old file11 format
    from molecule and gradient information in *wfn*.

    .. versionadded:: 0.6
       *wfn* parameter passed explicitly

    :returns: None

    :type filename: string
    :param filename: destination file name for file11 file

    :type wfn: :py:class:`~psi4.core.Wavefunction`
    :param wfn: set of molecule, gradient from which to generate file11

    :examples:

    >>> # [1] file11 for CISD calculation
    >>> G, wfn = gradient('cisd', return_wfn=True)
    >>> fake_file11(wfn, 'mycalc.11')

    """
    molecule = wfn.molecule()
    molecule.update_geometry()
    gradient = wfn.gradient()

    with open(filename, 'w') as handle:
        handle.write('%d\n' % (molecule.natom()))

        for at in range(molecule.natom()):
            handle.write('%6s %16.8f %16.8f %16.8f\n' % (molecule.symbol(
                at), molecule.x(at), molecule.y(at), molecule.z(at)))

        for at in range(molecule.natom()):
            handle.write('%6s %16.8f %16.8f %16.8f\n' % (
                '', gradient.get(at, 0), gradient.get(at, 1), gradient.get(at, 2)))


def sherrill_gold_standard(**kwargs):
    r"""Function to call the quantum chemical method known as 'Gold Standard'
    in the Sherrill group. Uses :py:func:`~driver_cbs.complete_basis_set` to evaluate
    the following expression. Two-point extrapolation of the correlation energy
    performed according to :py:func:`~driver_cbs.corl_xtpl_helgaker_2`.

    .. math:: E_{total}^{\text{Au\_std}} = E_{total,\; \text{SCF}}^{\text{aug-cc-pVQZ}} \; + E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[TQ]Z}} \; + \delta_{\text{MP2}}^{\text{CCSD(T)}}\big\vert_{\text{aug-cc-pVTZ}}

    >>> # [1] single-point energy by this composite method
    >>> energy('sherrill_gold_standard')

    >>> # [2] finite-difference geometry optimization
    >>> optimize('sherrill_gold_standard')

    >>> # [3] finite-difference geometry optimization, overwriting some pre-defined sherrill_gold_standard options
    >>> optimize('sherrill_gold_standard', corl_basis='cc-pV[DT]Z', delta_basis='3-21g')

    """
    scf = {
        'wfn': 'hf',
        'basis': kwargs.pop('scf_basis', 'aug-cc-pVQZ'),
        'scheme': kwargs.pop('scf_scheme', 'xtpl_highest_1')
    }
    corl = {
        'wfn': kwargs.pop('corl_wfn', 'mp2'),
        'basis': kwargs.pop('corl_basis', 'aug-cc-pV[TQ]Z'),
        'scheme': kwargs.pop('corl_scheme', 'corl_xtpl_helgaker_2')
    }
    delta = {
        'wfn': kwargs.pop('delta_wfn', 'ccsd(t)'),
        'wfn_lesser': kwargs.pop('delta_wfn_lesser', 'mp2'),
        'basis': kwargs.pop('delta_basis', 'aug-cc-pVTZ'),
        'scheme': kwargs.pop('delta_scheme', 'xtpl_highest_1')
    }

    return [scf, corl, delta]


def allen_focal_point(**kwargs):
    r"""Function to call Wes Allen-style Focal
    Point Analysis. JCP 127 014306.  Uses
    :py:func:`~driver_cbs.complete_basis_set` to evaluate the following
    expression. SCF employs a three-point extrapolation according
    to :py:func:`~driver_cbs.scf_xtpl_helgaker_3`. MP2, CCSD, and
    CCSD(T) employ two-point extrapolation performed according to
    :py:func:`~driver_cbs.corl_xtpl_helgaker_2`.  CCSDT and CCSDT(Q)
    are plain deltas. This wrapper requires :ref:`Kallay's MRCC code <sec:mrcc>`.

    .. math:: E_{total}^{\text{FPA}} = E_{total,\; \text{SCF}}^{\text{cc-pV[Q56]Z}} \; + E_{corl,\; \text{MP2}}^{\text{cc-pV[56]Z}} \; + \delta_{\text{MP2}}^{\text{CCSD}}\big\vert_{\text{cc-pV[56]Z}} \; + \delta_{\text{CCSD}}^{\text{CCSD(T)}}\big\vert_{\text{cc-pV[56]Z}} \; + \delta_{\text{CCSD(T)}}^{\text{CCSDT}}\big\vert_{\text{cc-pVTZ}} \; + \delta_{\text{CCSDT}}^{\text{CCSDT(Q)}}\big\vert_{\text{cc-pVDZ}}

    >>> # [1] single-point energy by this composite method
    >>> energy('allen_focal_point')

    >>> # [2] finite-difference geometry optimization embarrasingly parallel
    >>> optimize('allen_focal_point', mode='sow')

    """

    scf = {  # HF
        'wfn': 'hf',
        'basis': kwargs.pop('scf_basis', 'cc-pV[Q56]Z'),
        'scheme': kwargs.pop('scf_scheme', 'scf_xtpl_helgaker_3'),
    }
    corl = {  # MP2 - HF
        'wfn': kwargs.pop('corl_wfn', 'mp2'),
        'basis': kwargs.pop('corl_basis', 'cc-pV[56]Z'),
        'scheme': kwargs.pop('corl_scheme', 'corl_xtpl_helgaker_2'),
    }
    delta = {  # CCSD - MP2
        'wfn': kwargs.pop('delta_wfn', 'mrccsd'),
        'wfn_lesser': kwargs.pop('delta_wfn_lesser', 'mp2'),
        'basis': kwargs.pop('delta_basis', 'cc-pV[56]Z'),
        'scheme': kwargs.pop('delta_scheme', 'corl_xtpl_helgaker_2'),
    }
    delta2 = {  # CCSD(T) - CCSD
        'wfn': kwargs.pop('delta2_wfn', 'mrccsd(t)'),
        'wfn_lesser': kwargs.pop('delta2_wfn_lesser', 'mrccsd'),
        'basis': kwargs.pop('delta2_basis', 'cc-pV[56]Z'),
        'scheme': kwargs.pop('delta2_scheme', 'corl_xtpl_helgaker_2'),
    }
    delta3 = {  # CCSDT - CCSD(T)
        'wfn': kwargs.pop('delta3_wfn', 'mrccsdt'),
        'wfn_lesser': kwargs.pop('delta3_wfn_lesser', 'mrccsd(t)'),
        'basis': kwargs.pop('delta3_basis', 'cc-pVTZ'),
        'scheme': kwargs.pop('delta3_scheme', 'xtpl_highest_1'),
    }
    delta4 = {  # CCSDT(Q) - CCSDT
        'wfn': kwargs.pop('delta4_wfn', 'mrccsdt(q)'),
        'wfn_lesser': kwargs.pop('delta4_wfn_lesser', 'mrccsdt'),
        'basis': kwargs.pop('delta4_basis', 'cc-pVDZ'),
        'scheme': kwargs.pop('delta4_scheme', 'xtpl_highest_1'),
    }

    return [scf, corl, delta, delta2, delta3, delta4]
