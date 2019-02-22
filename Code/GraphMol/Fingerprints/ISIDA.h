//
//  Copyright (C) 2007-2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file ISIDA.h

*/
#include <RDGeneral/export.h>
#ifndef __RD_ATOMPAIRS_H__
#define __RD_ATOMPAIRS_H__

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitVects.h>
#include <cstdint>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
namespace RDKit {

namespace ISIDA {
const std::string isidaVersion = "0.1";

//! returns the atom-pair fingerprint for a molecule
/*!
  The algorithm used is described here:
  R.E. Carhart, D.H. Smith, R. Venkataraghavan; "Atom Pairs as
  Molecular Features in Structure-Activity Studies: Definition
  and Applications" JCICS 25, 64-73 (1985).


  \param mol:   the molecule to be fingerprinted
  \param minLength:   minimum distance between atoms to be
                      considered in a pair. Default is 1 bond.
  \param maxLength:   maximum distance between atoms to be
                      considered in a pair.
                      Default is maxPathLen-1 bonds.
  \param fromAtoms:   if provided, only atom pairs that involve
                      the specified atoms will be included in the
                      fingerprint
  \param ignoreAtoms: if provided, any atom pairs that include
                      the specified atoms will not be included in the
                      fingerprint
  \param atomInvariants: a list of invariants to use for the atom hashes
                         note: only the first \c codeSize bits of each
                         invariant are used.
  \param includeChirality: if set, chirality will be used in the atom invariants
                           (note: this is ignored if atomInvariants are
  provided)
  \param use2D:       if set, the 2D (topological) distance matrix is used.
  \param confId:      the conformation to use if 3D distances are being used


  \return a pointer to the fingerprint. The client is
  responsible for calling delete on this.

*/
RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::int32_t> *getISIDAFingerprint(
    const ROMol &mol, unsigned int minLength, unsigned int maxLength,
    const std::vector<std::uint32_t> *fromAtoms = 0,
    const std::vector<std::uint32_t> *ignoreAtoms = 0,
    const std::vector<std::uint32_t> *atomInvariants = 0,
    bool includeChirality = false, bool use2D = true, int confId = -1);

}  // namespace ISIDA 
}  // namespace RDKit

#endif
