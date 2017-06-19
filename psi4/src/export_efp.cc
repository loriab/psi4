/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/pybind11.h"

#include "psi4/libefp_solver/efp_solver.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libparallel/process.h"

using namespace psi;
using namespace psi::efp;


void construct_from_pydict(py::list pyefp) {

    size_t efp_nfrag = len(pyefp);
    std::vector <std::string> efp_fnames;
    enum efp_coord_type { XYZABC, POINTS, ROTMAT };
    efp_coord_type efp_ctype;
    std::vector <double> hint_coords;

// TODO       if (efp_nfrag > 0) {

    // Collect and initialize all fragment names and library paths
    for (size_t fr = 0; fr < efp_nfrag; ++fr) {
        py::dict frag_dict = pyefp[fr].cast<py::dict>();
        efp_fnames.push_back(frag_dict["fragment_file"].cast<std::string>());
    }
    Process::environment.get_efp()->add_fragments(efp_fnames);

    // Collect geometry hints
    for (size_t fr = 0; fr < efp_nfrag; ++fr) {
        py::dict frag_dict = pyefp[fr].cast<py::dict>();

        if (frag_dict["efp_type"].cast<std::string>() == "xyzabc")
            efp_ctype = XYZABC;
        else if (frag_dict["efp_type"].cast<std::string>() == "points")
            efp_ctype = POINTS;

        hint_coords = frag_dict["coordinates_hint"].cast<std::vector <double>>();
        Process::environment.get_efp()->set_frag_coordinates(fr, efp_ctype, &hint_coords[0]);
    }

    // Finalize efp fragment composition
    Process::environment.get_efp()->finalize_fragments();

    Process::environment.get_efp()->print_efp_geometry();
}


void export_efp(py::module& m) {
    py::class_<EFP, std::shared_ptr<EFP> >(m, "EFP", "Class interfacing with libefp")
        // because there is no default constructor for libefp
        .def(py::init<Options&>())
#ifdef USING_libefp
        .def("compute", &EFP::compute, "Computes libefp energies and, if active, torque")
        .def("set_qm_atoms", &EFP::set_qm_atoms, "Provides libefp with QM fragment information")
        .def("print_out", &EFP::print_out, "Prints options settings and EFP and QM geometries")
#endif
        .def_static("construct_from_pydict", construct_from_pydict, "docstring")
        .def("nfragments", &EFP::get_frag_count,
             "Returns the number of EFP fragments in the molecule");
}
