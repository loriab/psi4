/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libmints/typedefs.h"
//#include <boost/shared_ptr.hpp>
//#include <libmints/mints.h>
//#include <libqt/qt.h>
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
//#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include <math.h>
#include "fitter.h"

namespace psi {

DFChargeFitter::DFChargeFitter() :
    print_(0), debug_(0)
{
printf("constructed\n");
}
DFChargeFitter::~DFChargeFitter()
{
printf("destroyed\n");
}
SharedVector DFChargeFitter::fit()
{
    int naux = auxiliary_->nbf();
    int nso  = primary_->nbf();

    SharedVector d (new Vector("d", naux));

    double* dp = d->pointer();
    double** Dp = D_->pointer();

    /* 3-index */ {

    auto factory = std::make_shared<IntegralFactory>(auxiliary_,BasisSet::zero_ao_basis_set(),
        primary_,primary_);
    printf("factory\n");
//    auto eri = std::make_shared<TwoBodyAOInt>(factory->eri());
    std::shared_ptr<TwoBodyAOInt> eri(factory->eri());
    const double* buffer = eri->buffer();
    printf("eri\n");

    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        int nq = auxiliary_->shell(Q).nfunction();
        int sq = auxiliary_->shell(Q).function_index();
        for (int M = 0; M < primary_->nshell(); M++) {
            int nm = primary_->shell(M).nfunction();
            int sm = primary_->shell(M).function_index();
            for (int N = 0; N < primary_->nshell(); N++) {
                int nn = primary_->shell(N).nfunction();
                int sn = primary_->shell(N).function_index();

                eri->compute_shell(Q,0,M,N);
//                const double* buffer = eri->buffer();

                for (int oq = 0, index = 0; oq < nq; oq++) {
                    for (int om = 0; om < nm; om++) {
                        for (int on = 0; on < nn; on++, index++) {
                            dp[sq + oq] += Dp[sm + om][sn + on] * buffer[oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }
    }
    /* End 3-index */ }
    printf("end 3-index\n");
    /* 2-index */ {

    auto J = std::make_shared<Matrix>("J", naux, naux);
    double** Jp = J->pointer();
    printf("2i A\n");

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto factory2 = std::make_shared<IntegralFactory>(auxiliary_, zero, auxiliary_, zero);
    printf("2i B\n");
//    auto eri2 = std::make_shared<TwoBodyAOInt>(factory2->eri());
    std::shared_ptr<TwoBodyAOInt> eri2(factory2->eri());
    printf("2i C\n");

    printf("aux_nshell=%d prim_nshell=%d\n", auxiliary_->nshell(), primary_->nshell());
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        int nq = auxiliary_->shell(Q).nfunction();
        int sq = auxiliary_->shell(Q).function_index();
//        printf("2I %d %d %d\n", Q, nq, sq);
        for (int P = 0; P < auxiliary_->nshell(); P++) {
            int np = auxiliary_->shell(P).nfunction();
            int sp = auxiliary_->shell(P).function_index();

//            printf("2I: Q=%d P=%d\n", Q, P);
            eri2->compute_shell(Q,0,P,0);
    const double* buffer2 = eri2->buffer();
            for (int oq = 0; oq < nq; oq++) {
                for (int op = 0; op < np; op++) {
//                    printf("2I: Q=%d P=%d oq=%d op=%d, val=%f\n", Q, P, oq, op, buffer2[0]);//oq * np + op]);
                    Jp[sq + oq][sp + op] = buffer2[oq * np + op];
                }
            }
        }
    }

    printf("end 2-index\n");
    int info;
    info = C_DPOTRF('L',naux,Jp[0],naux);
    if (info) throw PSIEXCEPTION("DFChargeFitter: C_DPOTRF Failed");
    info = C_DPOTRS('L',naux,1,Jp[0],naux,dp,naux);
    if (info) throw PSIEXCEPTION("DFChargeFitter: C_DPOTRS Failed");

    /* End 2-index */ }
    printf("end linalg\n");

    d_ = d;
    return d;
}


} // Namespace psi
