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


py::dict to_dict(Molecule& mol) {

    py::dict mol_init;

    if ((mol.name() != "") && (mol.name() != "default"))
        mol_init[py::str("name")] = mol.name();
    if (mol.units() == Molecule::Angstrom)
        mol_init[py::str("units")] = "Angstrom";
    else if (mol.units() == Molecule::Bohr)
        mol_init[py::str("units")] = "Bohr";
    mol_init[py::str("input_units_to_au")] = mol.input_units_to_au();
    if (mol.charge_specified())
        mol_init[py::str("system_charge")] = mol.molecular_charge();
    if (mol.multiplicity_specified())
        mol_init[py::str("system_multiplicity")] = mol.multiplicity();
    if (mol.com_fixed())
        mol_init[py::str("fix_com")] = true;
    if (mol.orientation_fixed())
        mol_init[py::str("fix_orientation")] = true;
    if (mol.symmetry_from_input() != "")
        mol_init[py::str("fix_symmetry")] = mol.symmetry_from_input();

    // TODO zmat, geometry_variables
    py::list fragment_types_str;
    std::vector <Molecule::FragmentType> fragment_types_enum = mol.fragment_types();
    for (int ifr = 0; ifr < mol.nfragments(); ++ifr) {
        if (fragment_types_enum[ifr] == Molecule::Real)
            fragment_types_str.append("Real");
        else if (fragment_types_enum[ifr] == Molecule::Ghost)
            fragment_types_str.append("Ghost");
        else if (fragment_types_enum[ifr] == Molecule::Absent)
            fragment_types_str.append("Absent");
    }
    mol_init[py::str("fragments")] = mol.fragments();
    mol_init[py::str("fragment_types")] = fragment_types_str;
    mol_init[py::str("fragment_charges")] = mol.fragment_charges();
    mol_init[py::str("fragment_multiplicities")] = mol.fragment_multiplicities();

    py::list full_atoms;
    for (int iat = 0; iat < mol.nallatom(); ++iat) {
        py::dict at_init;
        //CoordEntry at_entry = mol.atom_entry(iat);
        at_init[py::str("Z")] = mol.fZ(iat);
        at_init[py::str("charge")] = mol.fcharge(iat);
        at_init[py::str("mass")] = mol.fmass(iat);
        at_init[py::str("symbol")] = mol.fsymbol(iat);
        at_init[py::str("label")] = mol.flabel(iat);
    //    atdict['ghosted'] = at.ghosted
        at_init[py::str("x")] = mol.fx(iat);
        at_init[py::str("y")] = mol.fy(iat);
        at_init[py::str("z")] = mol.fz(iat);
        // TODO xyz only for now
        full_atoms.append(at_init);
    }
    mol_init[py::str("full_atoms")] = full_atoms;

    return mol_init;
}

std::shared_ptr<Molecule> from_dict(py::dict pymol) {

    std::shared_ptr <Molecule> mol(new Molecule);

    if (pymol["units"].cast<std::string>() == "Angstrom")
        mol->set_units(Molecule::Angstrom);
    else if (pymol["units"].cast<std::string>() == "Bohr")
        mol->set_units(Molecule::Bohr);
    else
        printf("bad construct_from_pydict\n");
    mol->set_input_units_to_au(pymol["input_units_to_au"].cast<double>());

    for (auto item : pymol) {
        if (std::string(py::str(item.first)) == "name")
            mol->set_name(pymol["name"].cast<std::string>());
        if (std::string(py::str(item.first)) == "fix_com")
            mol->set_com_fixed(pymol["fix_com"].cast<bool>());
        if (std::string(py::str(item.first)) == "fix_orientation")
            mol->set_orientation_fixed(pymol["fix_orientation"].cast<bool>());
        if (std::string(py::str(item.first)) == "fix_symmetry")
            mol->reset_point_group(pymol["fix_symmetry"].cast<std::string>());
        if (std::string(py::str(item.first)) == "system_charge")
            mol->set_molecular_charge(pymol["system_charge"].cast<int>());
        if (std::string(py::str(item.first)) == "system_multiplicity")
            mol->set_multiplicity(pymol["system_multiplicity"].cast<int>());
        if (std::string(py::str(item.first)) == "full_atoms") {
            for (auto at : pymol["full_atoms"]) {
                mol->add_atom(at["Z"].cast<double>(),
                              at["x"].cast<double>(),
                              at["y"].cast<double>(),
                              at["z"].cast<double>(),
                              at["symbol"].cast<std::string>().c_str(),
                              at["mass"].cast<double>(),
                              at["charge"].cast<double>(),
                              at["label"].cast<std::string>().c_str());
            }
        }
    }

    std::vector<Molecule::FragmentType> fragment_types;
    for (auto item : pymol["fragment_types"]) {
        if (item.cast<std::string>() == "Real")
            fragment_types.push_back(Molecule::Real);
        else if (item.cast<std::string>() == "Ghost")
            fragment_types.push_back(Molecule::Ghost);
        else if (item.cast<std::string>() == "Absent")
            fragment_types.push_back(Molecule::Absent);
        else
            printf("bad2\n");
    }
    mol->set_has_zmatrix(false);  // TODO

    mol->set_fragment_pattern(pymol["fragments"].cast<std::vector <std::pair <int, int>>>(),
                              fragment_types,
                              pymol["fragment_charges"].cast<std::vector <int>>(),
                              pymol["fragment_multiplicities"].cast<std::vector <int>>());


//        moldict['units'] = self.PYunits
//        moldict['input_units_to_au'] = self.input_units_to_au
//        moldict['zmat'] = self.zmat
//
//        moldict['reinterpret_coordentries'] = self.PYreinterpret_coordentries
//        moldict['lock_frame'] = self.lock_frame
//
//        for at in self.full_atoms:
//            atdict['symbol'] = at.PYsymbol
//            atdict['label'] = at.PYlabel

//    mol->set_orientation_fixed(false);
//    mol->set_com_fixed(false);
//    mol->set_reinterpret_coordentry(false);

//    printf("ASAD1 %d\n", mol->natom());
//    mol->print_full();
//    mol->print_cluster();
//    mol->update_geometry();
//    printf("ASAD2 %d\n", mol->natom());
//    mol->print_full();
//    mol->print_cluster();
    return mol;
}


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
    typedef void (Molecule::*matrix_set_full_geometry)(const Matrix&);
    typedef Vector3 (Molecule::*nuclear_dipole1)(const Vector3&) const;
    typedef Vector3 (Molecule::*nuclear_dipole2)() const;

    py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule",
                                                    "Class to store the elements, coordinates, "
                                                    "fragmentation pattern, basis sets, charge, "
                                                    "multiplicity, etc. of a molecule.")
        .def("geometry", &Molecule::geometry,
             "Gets the geometry [Bohr] as a (Natom X 3) matrix of coordinates (excluding dummies)")
        .def("set_geometry", matrix_set_geometry(&Molecule::set_geometry),
             "Sets the geometry, given a (Natom X 3) matrix arg0 of coordinates (in Bohr) (excluding dummies")
        .def("full_geometry", &Molecule::full_geometry,
             "Gets the geometry [Bohr] as a (Natom X 3) matrix of coordinates (including dummies)")
        .def("set_full_geometry", matrix_set_full_geometry(&Molecule::set_full_geometry),
             "Sets the geometry, given a (Natom X 3) matrix arg0 of coordinates (in Bohr) (including dummies")
        .def("nuclear_dipole", nuclear_dipole1(&Molecule::nuclear_dipole),
             "Gets the nuclear contribution to the dipole, with respect to a specified origin arg0")
        .def("nuclear_dipole", nuclear_dipole2(&Molecule::nuclear_dipole),
             "Gets the nuclear contribution to the dipole, with respect to the origin")
        .def("name", &Molecule::name, "Gets molecule name")
        .def("set_name", &Molecule::set_name, "Sets molecule name")
        .def("fix_orientation", &Molecule::set_orientation_fixed,
             "Fix the orientation at its current frame")
        .def("orientation_fixed", &Molecule::orientation_fixed,
             "Get whether or not orientation is fixed")
        .def("fix_com", &Molecule::set_com_fixed,
             "Sets whether to fix the Cartesian position, or to translate to the C.O.M.")
        .def("com_fixed", &Molecule::com_fixed,
             "Gets whether or not center of mass is fixed")
        .def("add_atom", &Molecule::add_atom,
             "Adds to self Molecule an atom with atomic number arg0, Cartesian coordinates in Bohr "
             "(arg1, arg2, arg3), atomic symbol arg4, mass arg5, charge arg6 (optional), and "
             "lineno arg7 (optional)")
        .def("natom", &Molecule::natom, "Number of real atoms")
        .def("nallatom", &Molecule::nallatom, "Number of real and dummy atoms")
        .def("molecular_charge", &Molecule::molecular_charge, "Gets the molecular charge")
        .def("set_molecular_charge", &Molecule::set_molecular_charge, "Sets the molecular charge")
        .def_property("chg", &Molecule::molecular_charge,
                             &Molecule::set_molecular_charge,
                      "Molecular charge")
        .def("charge_specified", &Molecule::charge_specified,
             "Gets whether the charge was given by the user")
        .def("multiplicity", &Molecule::multiplicity, "Gets the multiplicity (defined as 2Ms + 1)")
        .def("set_multiplicity", &Molecule::set_multiplicity,
             "Sets the multiplicity (defined as 2Ms + 1)")
        .def_property("mult", &Molecule::multiplicity,
                              &Molecule::set_multiplicity,
                      "Molecular multiplicity (defined as 2Ms + 1)")
        .def("multiplicity_specified", &Molecule::multiplicity_specified,
             "Gets whether the multiplicity was given by the user")
        .def("nfragments", &Molecule::nfragments, "Gets the number of fragments in the molecule")
        .def("n_active_fragments", &Molecule::nactive_fragments, "Gets the number of active fragments in the molecule")
        .def("fragment_atom_pair", &Molecule::fragment_atom_pair,
             "Gets pair with index of first atom and last atom + 1 belonging to fragment arg0")
        .def("set_nuclear_charge", &Molecule::set_nuclear_charge,
             "Set the nuclear charge of the given atom arg0 to the value arg1 (primarily for ECP).")
        .def("basis_on_atom", &Molecule::basis_on_atom, py::return_value_policy::copy,
             "Gets the label of the orbital basis set on a given atom arg0")
        .def("print_in_input_format", &Molecule::print_in_input_format,
             "Prints the molecule as Cartesian or ZMatrix entries, just as inputted.")
        .def("create_psi4_string_from_molecule", &Molecule::create_psi4_string_from_molecule,
             "Gets a string re-expressing in input format the current state of the molecule."
             "Contains Cartesian geometry info, fragmentation, charges and multiplicities, "
             "and any frame restriction.")
        .def("save_xyz_file", &Molecule::save_xyz_file, "Saves an XYZ file to arg0")
        .def("save_string_xyz_file", &Molecule::save_string_xyz_file, "Returns string of an XYZ file")
        .def("save_string_xyz", &Molecule::save_string_xyz,
             "Returns string of an XYZ file that includes charge and multiplicity")
        .def("Z", &Molecule::Z, py::return_value_policy::copy,
             "Nuclear charge of atom arg0 (0-indexed without dummies)")
        .def("x", &Molecule::x, "x position [Bohr] of atom arg0 (0-indexed without dummies)")
        .def("y", &Molecule::y, "y position [Bohr] of atom arg0 (0-indexed without dummies)")
        .def("z", &Molecule::z, "z position [Bohr] of atom arg0 (0-indexed without dummies)")
        .def("fZ", &Molecule::Z, py::return_value_policy::copy,
             "Nuclear charge of atom arg0 (0-indexed including dummies)")
        .def("fx", &Molecule::x, "x position of atom arg0 (0-indexed including dummies in Bohr)")
        .def("fy", &Molecule::y, "y position of atom arg0 (0-indexed including dummies in Bohr)")
        .def("fz", &Molecule::z, "z position of atom arg0 (0-indexed including dummies in Bohr)")
        .def("center_of_mass", &Molecule::center_of_mass,
             "Returns center of mass of molecule (does not translate molecule)")
        .def("translate", &Molecule::translate, "Translates molecule by arg0")
        .def("move_to_com", &Molecule::move_to_com, "Moves molecule to center of mass")
        .def("mass", &Molecule::mass, "Gets mass of atom arg0")
        .def("set_mass", &Molecule::set_mass, "Gets mass of atom arg0 to value arg1 (good for "
             "isotopic substitutions)")
        .def("fmass", &Molecule::fmass, "Gets mass of atom arg0 (0-indexed including dummies)")
        .def("symbol", &Molecule::symbol,
             "Gets the cleaned up label of atom arg0 (C2 => C, H4 = H) (0-indexed without dummies)")
        .def("fsymbol", &Molecule::fsymbol,
             "Gets the cleaned up label of atom arg0 (C2 => C, H4 = H) (0-indexed including dummies)")
        .def("label", &Molecule::label,
             "Gets the original label of the atom arg0 as given in the input file (C2, H4)"
             "(0-indexed without dummies)")
        .def("flabel", &Molecule::flabel,
             "Gets the original label of the atom arg0 as given in the input file (C2, H4)"
             "(0-indexed including dummies)")
        .def("charge", &Molecule::charge, "Gets charge of atom arg0 (0-indexed without dummies)")
        .def("fcharge", &Molecule::fcharge, "Gets charge of atom arg0 (0-indexed including dummies)")
        .def("true_atomic_number", &Molecule::true_atomic_number, "Gets atomic number of "
             "atom arg0 from element (0-indexed without dummies)")
        .def("ftrue_atomic_number", &Molecule::ftrue_atomic_number, "Gets atomic number of "
             "atom arg0 from element (0-indexed including dummies)")
        .def("extract_subsets", &Molecule::py_extract_subsets_1,
             "Returns copy of self with arg0 fragments Real and arg1 fragments Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_2,
             "Returns copy of self with arg0 fragments Real and arg1 fragment Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_3,
             "Returns copy of self with arg0 fragment Real and arg1 fragments Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_4,
             "Returns copy of self with arg0 fragment Real and arg1 fragment Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_5,
             "Returns copy of self with arg0 fragments Real")
        .def("extract_subsets", &Molecule::py_extract_subsets_6,
             "Returns copy of self with arg0 fragment Real")
        .def("activate_all_fragments", &Molecule::activate_all_fragments,
             "Sets all fragments in the molecule to be active")
        .def("deactivate_all_fragments", &Molecule::deactivate_all_fragments,
             "Sets all fragments in the molecule to be inactive")
        .def("set_active_fragments", &Molecule::set_active_fragments,
             "Sets the specified list arg0 of fragments to be Real")
        .def("set_active_fragment", &Molecule::set_active_fragment,
             "Sets the specified fragment arg0 to be Real")
        .def("set_ghost_fragments", &Molecule::set_ghost_fragments,
             "Sets the specified list arg0 of fragments to be Ghost")
        .def("set_ghost_fragment", &Molecule::set_ghost_fragment,
             "Sets the specified fragment arg0 to be Ghost")
        .def("fragments", &Molecule::fragments,
              "Returns list of pairs of atom ranges defining each fragment from parent molecule"
              "(fragments[frag_ind] = <Afirst,Alast+1>)")
        .def("fragment_types", &Molecule::fragment_types,
              "A list describing how to handle each fragment {Real, Ghost, Absent}")
        .def("fragment_charges", &Molecule::fragment_charges,
              "Gets the charge of each fragment")
        .def("fragment_multiplicities", &Molecule::fragment_multiplicities,
              "Gets the multiplicity of each fragment")
        .def("atom_at_position", &Molecule::atom_at_position1,
             "Tests to see if an atom is at the position arg0 with a given tolerance arg1")
        .def("print_out", &Molecule::print, "Prints the molecule in Cartesians in input units")
        .def("print_out_in_bohr", &Molecule::print_in_bohr,
             "Prints the molecule in Cartesians in Bohr")
        .def("print_out_in_angstrom", &Molecule::print_in_angstrom,
             "Prints the molecule in Cartesians in Angstroms")
        .def("print_cluster", &Molecule::print_cluster,
             "Prints the molecule in Cartesians in input units adding fragment separators")
        .def("rotational_constants", &Molecule::rotational_constants,
             "Returns the rotational constants [cm^-1] of the molecule")
        .def("print_rotational_constants", &Molecule::print_rotational_constants,
             "Print the rotational constants")
        .def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy,
             "Computes nuclear repulsion energy")
        .def("nuclear_repulsion_energy_deriv1", &Molecule::nuclear_repulsion_energy_deriv1,
             "Computes nuclear repulsion energy first derivatives")
        .def("nuclear_repulsion_energy_deriv2", &Molecule::nuclear_repulsion_energy_deriv2,
             "Computes nuclear repulsion energy second derivatives")
        .def("find_point_group", &Molecule::find_point_group,
             "Finds computational molecular point group, user can override this with the symmetry "
             "keyword")
        .def("reset_point_group", &Molecule::reset_point_group,
             "Overrides symmetry by arg0 from outside the molecule string")
        .def("set_point_group", &Molecule::set_point_group,
             "Sets the molecular point group to the point group object arg0")
        .def("symmetry_specified", &Molecule::symmetry_from_input,
             "docstring")
        .def("get_full_point_group", &Molecule::full_point_group,
             "Gets point group name such as C3v or S8")
        .def("get_full_point_group_with_n", &Molecule::full_point_group_with_n,
             "Gets point group name such as Cnv or Sn")
        .def("full_pg_n", &Molecule::full_pg_n,
             "Gets n in Cnv, etc.; If there is no n (e.g. Td) it's the highest-order rotation axis")
        .def("point_group", &Molecule::point_group, "Returns the current point group object")
        .def("schoenflies_symbol", &Molecule::schoenflies_symbol, "Returns the Schoenflies symbol")
        .def("form_symmetry_information", &Molecule::form_symmetry_information,
             "Uses the point group object obtain by calling point_group()")
        .def("symmetrize", &Molecule::symmetrize_to_abelian_group,
             "Finds the highest point Abelian point group within the specified tolerance arg0, and "
             "forces the geometry to have that symmetry.")
        .def_static(
             "create_molecule_from_string", &Molecule::create_molecule_from_string,
             "Returns a new Molecule with member data from the geometry string arg0 in psi4 format")
        .def("is_variable", &Molecule::is_variable,
             "Checks if variable arg0 is in the structural variable list")
        .def("set_variable", &Molecule::set_variable,
             "Sets the value arg1 to the variable arg0 in the list of structure variables, then "
             "calls update_geometry()")
        .def("get_variable", &Molecule::get_variable,
             "Returns the value of variable arg0 in the structural variables list")
        .def("update_geometry", &Molecule::update_geometry,
             "Reevaluates the geometry with current variable values, orientation directives, etc. "
             "by clearing the atoms list and rebuilding it. Idempotent. Use liberally."
             "Must be called after initial Molecule definition by string.")
        .def_property("zmat", &Molecule::has_zmatrix,
                              &Molecule::set_has_zmatrix,
             "Whether this molecule contains at least one zmatrix entry")
        .def("set_basis_all_atoms", &Molecule::set_basis_all_atoms,
             "Sets basis set arg0 to all atoms")
        .def("set_basis_by_symbol", &Molecule::set_basis_by_symbol,
             "Sets basis set arg1 to all atoms with symbol (e.g., H) arg0")
        .def("set_basis_by_label", &Molecule::set_basis_by_label,
             "Sets basis set arg1 to all atoms with label (e.g., H4) arg0")
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
                      "Units (Angstrom or Bohr) used to define the geometry. Imposes Psi4 physical constants conversion for input_units_to_au.")
        .def("input_units_to_au", &Molecule::input_units_to_au,
             "Gets unit conversion to Bohr for geometry")
        .def("set_input_units_to_au", &Molecule::set_input_units_to_au,
             "Sets unit conversion to Bohr for geometry to *arg0*")
        .def("to_dict", to_dict, "Gets python dictionary representation of instance")
        .def_static("from_dict", from_dict,
                    "Returns a new Molecule constructed from entires in python dictionary arg0")
        .def("clone", &Molecule::clone, "Returns a new Molecule identical to arg0");
}


//    /// Number of frozen core for molecule given freezing state
//    int nfrozen_core(std::shared_ptr<BasisSet> ecpbasis, const std::string& depth = "");
//
//    /// Sets the frame so update_geometry can act correctly
//    void set_lock_frame(bool tf) { lock_frame_ = tf; }
//
//    /// Do we reinterpret coordentries during a call to update_geometry?
//    void set_reinterpret_coordentry(bool rc);
//
//    /**
//     * Rotates the molecule using rotation matrix R
//     */
//    void rotate(const Matrix& R);
//    void rotate_full(const Matrix& R);
//
//    /**
//     * Reinterpret the fragments for reals/ghosts and build the atom list
//     */
//    void reinterpret_coordentries();
//
//    /**
//     *  Reorient molecule to standard frame. See input/reorient.cc
//     *  If you want the molecule to be reoriented about the center of mass
//     *  make sure you call move_to_com() prior to calling reorient()
//     */
////    void reorient();
//
//    /// Compute inertia tensor.
//    Matrix* inertia_tensor() const;
//
//
//    /// Return the rotor type
//    RotorType rotor_type(double tol = FULL_PG_TOL) const;
//
//
