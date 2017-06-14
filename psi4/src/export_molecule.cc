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
#include "psi4/physconst.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/corrtab.h"


#include <string>

using namespace psi;


void export_molecule(py::module& m)
{
    py::enum_<Molecule::GeometryUnits>(m, "GeometryUnits", "docstring")
        .value("Angstrom", Molecule::Angstrom)
        .value("Bohr", Molecule::Bohr)
        .export_values();

    py::enum_<Molecule::FragmentType>(m, "FragmentType",
        "Fragment status of neglect completely, include normally, or include as ghost atoms")
        .value("Absent", Molecule::Absent)
        .value("Real", Molecule::Real)
        .value("Ghost", Molecule::Ghost)
        .export_values();


    typedef void (Molecule::*matrix_set_geometry)(const Matrix&);
    typedef Vector3 (Molecule::*nuclear_dipole1)(const Vector3&) const;
    typedef Vector3 (Molecule::*nuclear_dipole2)() const;

    py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule",
                                                    "Class to store the elements, coordinates, "
                                                    "fragmentation pattern, basis sets, charge, "
                                                    "multiplicity, etc. of a molecule.")
        .def("set_geometry", matrix_set_geometry(&Molecule::set_geometry),
             "Sets the geometry, given a (Natom X 3) matrix arg2 of coordinates (in Bohr)")
        .def("nuclear_dipole", nuclear_dipole1(&Molecule::nuclear_dipole),
             "Gets the nuclear contribution to the dipole, withe respect to a specified origin")
        .def("nuclear_dipole", nuclear_dipole2(&Molecule::nuclear_dipole),
             "Gets the nuclear contribution to the dipole, withe respect to the origin")
        .def("set_name", &Molecule::set_name, "Sets molecule name")
        .def("name", &Molecule::name, "Gets molecule name")
        .def("reinterpret_coordentry", &Molecule::set_reinterpret_coordentry,
             "Do reinterpret coordinate entries during update_geometry().")
        .def("fix_orientation", &Molecule::set_orientation_fixed,
             "Fix the orientation at its current frame")
        .def("fix_com", &Molecule::set_com_fixed,
             "Whether to fix the Cartesian position, or to translate to the C.O.M.")
        .def("add_atom", &Molecule::add_atom,
             "Adds to Molecule arg1 an atom with atomic number arg2, Cartesian coordinates in Bohr "
             "(arg3, arg4, arg5), atomic symbol arg6, mass arg7, charge arg8 (optional), and "
             "lineno arg9 (optional)")
        .def("natom", &Molecule::natom, "Number of real atoms")
        .def("nallatom", &Molecule::nallatom, "Number of real and dummy atoms")
        .def("multiplicity", &Molecule::multiplicity, "Gets the multiplicity (defined as 2Ms + 1)")
        .def("nfragments", &Molecule::nfragments, "Gets the number of fragments in the molecule")
        .def("set_nuclear_charge", &Molecule::set_nuclear_charge,
             "Set the nuclear charge of the given atom to the value provided.")
        .def("basis_on_atom", &Molecule::basis_on_atom, py::return_value_policy::copy,
             "Gets the label of the orbital basis set on a given atom.")
        .def("print_in_input_format", &Molecule::print_in_input_format,
             "Prints the molecule as Cartesian or ZMatrix entries, just as inputted.")
        .def("create_psi4_string_from_molecule", &Molecule::create_psi4_string_from_molecule,
             "Gets a string reexpressing in input format the current states of the molecule")
        .def("save_xyz_file", &Molecule::save_xyz_file, "Saves an XYZ file to arg2")
        .def("save_string_xyz_file", &Molecule::save_string_xyz_file, "Saves an XYZ file to arg2")
        .def("save_string_xyz", &Molecule::save_string_xyz,
             "Saves the string of an XYZ file to arg2")
        .def("Z", &Molecule::Z, py::return_value_policy::copy, "Nuclear charge of atom")
        .def("x", &Molecule::x, "x position of atom")
        .def("y", &Molecule::y, "y position of atom")
        .def("z", &Molecule::z, "z position of atom")
        .def("fZ", &Molecule::Z, py::return_value_policy::copy,
             "Nuclear charge of atom arg1 (0-indexed including dummies)")
        .def("fx", &Molecule::x, "x position of atom arg1 (0-indexed including dummies in Bohr)")
        .def("fy", &Molecule::y, "y position of atom arg1 (0-indexed including dummies in Bohr)")
        .def("fz", &Molecule::z, "z position of atom arg1 (0-indexed including dummies in Bohr)")
        .def("center_of_mass", &Molecule::center_of_mass,
             "Computes center of mass of molecule (does not translate molecule)")
        .def("translate", &Molecule::translate, "Translates molecule by arg2")
        .def("move_to_com", &Molecule::move_to_com, "Moves molecule to center of mass")
        .def("mass", &Molecule::mass, "Gets mass of atom arg2")
        .def("set_mass", &Molecule::set_mass, "Gets mass of atom arg2")
        .def("symbol", &Molecule::symbol,
             "Gets the cleaned up label of atom arg2 (C2 => C, H4 = H)")
        .def("label", &Molecule::label,
             "Gets the original label of the atom as given in the input file (C2, H4)")
        .def("charge", &Molecule::charge, "Gets charge of atom")
        .def("molecular_charge", &Molecule::molecular_charge, "Gets the molecular charge")
        .def("extract_subsets", &Molecule::py_extract_subsets_1,
             "Returns copy of arg1 with arg2 fragments Real and arg3 fragments Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_2,
             "Returns copy of arg1 with arg2 fragments Real and arg3 fragment Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_3,
             "Returns copy of arg1 with arg2 fragment Real and arg3 fragments Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_4,
             "Returns copy of arg1 with arg2 fragment Real and arg3 fragment Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_5,
             "Returns copy of arg1 with arg2 fragments Real")
        .def("extract_subsets", &Molecule::py_extract_subsets_6,
             "Returns copy of arg1 with arg2 fragment Real")
        .def("activate_all_fragments", &Molecule::activate_all_fragments,
             "Sets all fragments in the molecule to be active")
        .def("deactivate_all_fragments", &Molecule::deactivate_all_fragments,
             "Sets all fragments in the molecule to be inactive")
        .def("set_active_fragments", &Molecule::set_active_fragments,
             "Sets the specified list arg2 of fragments to be Real")
        .def("set_active_fragment", &Molecule::set_active_fragment,
             "Sets the specified fragment arg2 to be Real")
        .def("set_ghost_fragments", &Molecule::set_ghost_fragments,
             "Sets the specified list arg2 of fragments to be Ghost")
        .def("set_ghost_fragment", &Molecule::set_ghost_fragment,
             "Sets the specified fragment arg2 to be Ghost")
        .def("atom_at_position", &Molecule::atom_at_position1,
             "Tests to see if an atom is at the position arg2 with a given tolerance arg3")
        .def("print_out", &Molecule::print, "Prints the molecule in Cartesians in input units")
        .def("print_out_in_bohr", &Molecule::print_in_bohr,
             "Prints the molecule in Cartesians in Bohr")
        .def("print_out_in_angstrom", &Molecule::print_in_angstrom,
             "Prints the molecule in Cartesians in Angstroms")
        .def("print_cluster", &Molecule::print_cluster,
             "Prints the molecule in Cartesians in input units adding fragment separators")
        .def("rotational_constants", &Molecule::rotational_constants,
             "Prints the rotational constants of the molecule")
        .def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy,
             "Computes nuclear repulsion energy")
        .def("find_point_group", &Molecule::find_point_group,
             "Finds computational molecular point group, user can override this with the symmetry "
             "keyword")
        .def("reset_point_group", &Molecule::reset_point_group,
             "Overrides symmetry from outside the molecule string")
        .def("set_point_group", &Molecule::set_point_group,
             "Sets the molecular point group to the point group object arg2")
        .def("get_full_point_group", &Molecule::full_point_group,
             "Gets point group name such as C3v or S8")
        .def("point_group", &Molecule::point_group, "Returns the current point group object")
        .def("schoenflies_symbol", &Molecule::schoenflies_symbol, "Returns the Schoenflies symbol")
        .def("form_symmetry_information", &Molecule::form_symmetry_information,
             "Uses the point group object obtain by calling point_group()")
        .def("symmetrize", &Molecule::symmetrize_to_abelian_group,
             "Finds the highest point Abelian point group within the specified tolerance, and "
             "forces the geometry to have that symmetry.")
        .def_static(
             "create_molecule_from_string", &Molecule::create_molecule_from_string,
             "Returns a new Molecule with member data from the geometry string arg1 in psi4 format")
        .def("is_variable", &Molecule::is_variable,
             "Checks if variable arg2 is in the list, returns true if it is, and returns false if "
             "not")
        .def("set_variable", &Molecule::set_variable,
             "Assigns the value arg3 to the variable arg2 in the list of geometry variables, then "
             "calls update_geometry()")
        .def("get_variable", &Molecule::get_variable,
             "Checks if variable arg2 is in the list, sets it to val and returns true if it is, "
             "and returns false if not")
        .def("update_geometry", &Molecule::update_geometry,
             "Reevaluates the geometry with current variable values, orientation directives, etc. "
             "Must be called after initial Molecule definition by string.")
        .def("set_molecular_charge", &Molecule::set_molecular_charge, "Sets the molecular charge")
        .def("set_multiplicity", &Molecule::set_multiplicity,
             "Sets the multiplicity (defined as 2Ms + 1)")
        .def("set_basis_all_atoms", &Molecule::set_basis_all_atoms,
             "Sets basis set arg2 to all atoms")
        .def("set_basis_by_symbol", &Molecule::set_basis_by_symbol,
             "Sets basis set arg3 to all atoms with symbol (e.g., H) arg2")
        .def("set_basis_by_label", &Molecule::set_basis_by_label,
             "Sets basis set arg3 to all atoms with label (e.g., H4) arg2")
        .def("distance_matrix", &Molecule::distance_matrix, "Returns Matrix of interatom distances")
        .def("print_distances", &Molecule::print_distances,
             "Print the interatomic distance geometrical parameters")
        .def("print_bond_angles", &Molecule::print_bond_angles,
             "Print the bond angle geometrical parameters")
        .def("print_out_of_planes", &Molecule::print_out_of_planes,
             "Print the out-of-plane angle geometrical parameters")
        .def("irrep_labels",
             [](Molecule& mol) {
                 std::vector<std::string> ret;
                 char** labels = mol.irrep_labels();
                 int nirrep = mol.point_group()->char_table().nirrep();
                 for (size_t h = 0; h < nirrep; h++) {
                     std::string lh(labels[h]);
                     ret.push_back(lh);
                 }
                 return ret;
             })
        .def_property("units", py::cpp_function(&Molecule::units),
                      py::cpp_function(&Molecule::set_units),
                      "Units (Angstrom or Bohr) used to define the geometry")
        .def("clone", &Molecule::clone, "Returns a new Molecule identical to arg1")
        .def("geometry", &Molecule::geometry,
             "Gets the geometry as a (Natom X 3) matrix of coordinates (in Bohr)");

}
