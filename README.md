# <img src="https://github.com/psi4/psi4media/blob/master/logos-psi4/psi4square.png" height=150>

| **Status** | [![Travis build](https://img.shields.io/travis/psi4/psi4/master.svg?logo=linux)](https://travis-ci.org/psi4/psi4) [![Appveyor build](https://img.shields.io/appveyor/ci/psi4/psi4.svg?logo=windows)](https://ci.appveyor.com/project/psi4/psi4) [![Codecov coverage](https://codecov.io/gh/psi4/psi4/branch/master/graph/badge.svg)](https://codecov.io/gh/psi4/psi4) [![LGTM analysis](https://img.shields.io/lgtm/grade/python/g/psi4/psi4.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/psi4/psi4/context:python) |
| :------ | :------- |
| **Latest Release** | [![Last release tag](https://img.shields.io/github/release/psi4/psi4.svg)](https://github.com/psi4/psi4/releases)  [![Commits since release](https://img.shields.io/github/commits-since/psi4/psi4/v1.2.svg)](https://github.com/psi4/psi4/releases/tag/v1.2) [![python](https://img.shields.io/badge/python-2.7%2C%203.5%2C%203.6-blue.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) |
| **Communication** | [![User site](https://img.shields.io/badge/home-Psi4-5077AB.svg)](http://www.psicode.org) [![docs latest](https://img.shields.io/badge/docs-latest-5077AB.svg?logo=read%20the%20docs)](http://psicode.org/psi4manual/master/index.html) [![chat on forum](https://img.shields.io/badge/chat-on_forum-808493.svg)](http://forum.psicode.org/) [![dev chat on slack](https://img.shields.io/badge/dev_chat-on_slack-808493.svg?logo=slack)](https://join.slack.com/t/psi4/shared_invite/enQtNDUyOTYzNTE0NjQ3LWExZDhkY2U4MTM1ZDZlNTBkNjMyMDcxZmFkN2NmYmZkMzliNzY2ZDc2OTBlYTk5ZTA2OGRkNWYxNzJmN2QyYWM) |
| **Foundation** | [![license](https://img.shields.io/github/license/psi4/psi4.svg)](https://opensource.org/licenses/LGPL-3.0) [![platforms](https://img.shields.io/badge/Platforms-Linux%2C%20MacOS%2C%20Windows%20WSL-orange.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) [![python](https://img.shields.io/badge/python-3.5+-blue.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) |
| **Installation** | [![obtain latest](https://img.shields.io/badge/obtain-latest-green.svg)](http://vergil.chemistry.gatech.edu/nu-psicode/install-v1.2.1.html) [![Conda](https://img.shields.io/conda/v/psi4/psi4.svg)](https://anaconda.org/psi4/psi4) [![Anaconda-Server Badge](https://anaconda.org/psi4/psi4/badges/latest_release_relative_date.svg)](https://anaconda.org/psi4/psi4) |

Psi4 is an open-source suite of *ab initio* quantum chemistry programs
designed for efficient, high-accuracy simulations of a variety of
molecular properties. We can routinely perform computations with more
than 2500 basis functions running serially or on multi-core machines.

With computationally demanding portions written in C++, exports
of many C++ classes into Python via Pybind11, and a flexible Python driver, Psi4
strives to be friendly to both users and developers.

* **Users' Website**  www.psicode.org

* **Downloading and Installing Psi4** http://psicode.org/psi4manual/master/build_faq.html (for the CMake adept, see [CMakeLists.txt](CMakeLists.txt)

* **Manual**  [http://bit.ly/psi4manual](http://psicode.org/psi4manual/master/index.html) (built nightly from master branch)

* **Tutorial** http://psicode.org/psi4manual/master/tutorial.html

* **Forum** http://forum.psicode.org

* **Communication & Support** http://psicode.org/psi4manual/master/introduction.html#technical-support

* **Github**  https://github.com/psi4/psi4 (authoritative repository)

* **Travis CI build status** [![Build Status](https://travis-ci.org/psi4/psi4.svg?branch=master)](https://travis-ci.org/psi4/psi4)

* **Anaconda**  https://anaconda.org/psi4 (binary available for Linux, Mac, and WSL Windows [![Binstar Badge](https://anaconda.org/psi4/psi4/badges/downloads.svg)](https://anaconda.org/psi4/psi4) ) [instructions](http://psicode.org/psi4manual/master/conda.html#how-to-install-a-psi4-binary-with-the-psi4conda-installer-download-site)

* **CodeCov** Test case coverage for quicktests. [![codecov](https://codecov.io/gh/psi4/psi4/branch/master/graph/badge.svg)](https://codecov.io/gh/psi4/psi4) (More like 70% -- we'll get this back soon.)

* **Interested Developers**  http://psicode.org/developers.php (welcome to fork psi4/psi4 and follow [GitHub contribution procedure](http://psicode.org/psi4manual/master/build_obtaining.html#faq-githubworkflow))

* **Sample Inputs**  http://www.psicode.org/psi4manual/master/testsuite.html (also in share/psi4/samples)

* **Download Tarball** https://github.com/psi4/psi4/releases 

* **Build Dashboard** https://testboard.org/cdash/index.php?project=Psi

* **YouTube Channel** https://www.youtube.com/psitutorials


License
=======

Psi4: an open-source quantum chemistry software package

Copyright (c) 2007-2018 The Psi4 Developers.

The copyrights for code used from other parties are included in
the corresponding files.

Psi4 is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, version 3.

Psi4 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with Psi4; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

The full text of the GNU Lesser General Public License (version 3) is included in the
COPYING.LESSER file of this repository, and can also be found
[here](https://www.gnu.org/licenses/lgpl.txt).


Citation
========

The journal article reference describing Psi4 is:

R. M. Parrish, L. A. Burns, D. G. A. Smith, A. C. Simmonett,
A. E. DePrince III, E. G. Hohenstein, U. Bozkaya, A. Yu. Sokolov,
R. Di Remigio, R. M. Richard, J. F. Gonthier, A. M. James,
H. R. McAlexander, A. Kumar, M. Saitow, X. Wang, B. P. Pritchard,
P. Verma, H. F. Schaefer III, K. Patkowski, R. A. King, E. F. Valeev,
F. A. Evangelista, J. M. Turney, T. D. Crawford, and C. D. Sherrill,
J. Chem. Theory Comput. 13(7) 3185--3197 (2017).
[doi: 10.1021/acs.jctc.7b00174](http://dx.doi.org/10.1021/acs.jctc.7b00174).

