//
//  Copyright (C) 2007-2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/ISIDA.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <DataStructs/SparseIntVect.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace ISIDA {

template <typename T1, typename T2>
void updateElement(SparseIntVect<T1> &v, T2 elem) {
  v.setVal(elem, v.getVal(elem) + 1);
}

template <typename T1>
void updateElement(ExplicitBitVect &v, T1 elem) {
  v.setBit(elem % v.getNumBits());
}

template <typename T>
void setAtomPairBit(std::uint32_t i, std::uint32_t j,
                    std::uint32_t nAtoms,
                    const std::vector<std::uint32_t> &atomCodes,
                    const double *dm, T *bv, unsigned int minLength,
                    unsigned int maxLength, bool includeChirality) {
  unsigned int dist = static_cast<unsigned int>(floor(dm[i * nAtoms + j]));
  if (dist >= minLength && dist <= maxLength) {
    std::uint32_t bitId =
        getAtomPairCode(atomCodes[i], atomCodes[j], dist, includeChirality);
    updateElement(*bv, static_cast<std::uint32_t>(bitId));
  }
}


SparseIntVect<std::int32_t> *getISIDAFingerprint(
    const ROMol &mol, unsigned int minLength, unsigned int maxLength,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality,
    bool use2D, int confId) {
  PRECONDITION(minLength <= maxLength, "bad lengths provided");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;

  if (use2D) {
    dm = MolOps::getDistanceMat(*lmol);
  } else {
    dm = MolOps::get3DDistanceMat(*lmol, confId);
  }
  const unsigned int nAtoms = lmol->getNumAtoms();

  std::vector<std::uint32_t> atomCodes;
  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    } else {
      atomCodes.push_back((*atomInvariants)[(*atomItI)->getIdx()] %
                          ((1 << codeSize) - 1));
    }
  }

  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    unsigned int i = (*atomItI)->getIdx();
    if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(), i) !=
                           ignoreAtoms->end()) {
      continue;
    }
    if (!fromAtoms) {
      for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
           atomItJ != lmol->endAtoms(); ++atomItJ) {
        unsigned int j = (*atomItJ)->getIdx();
        if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                     j) != ignoreAtoms->end()) {
          continue;
        }
        setAtomPairBit(i, j, nAtoms, atomCodes, dm, res, minLength, maxLength,
                       includeChirality);
      }
    } else {
      BOOST_FOREACH (std::uint32_t j, *fromAtoms) {
        if (j != i) {
          if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                       j) != ignoreAtoms->end()) {
            continue;
          }
          setAtomPairBit(i, j, nAtoms, atomCodes, dm, res, minLength, maxLength,
                         includeChirality);
        }
      }
    }
  }
  return res;
}

SparseIntVect<std::int32_t> *getHashedAtomPairFingerprint(
    const ROMol &mol, unsigned int nBits, unsigned int minLength,
    unsigned int maxLength, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality,
    bool use2D, int confId) {
  PRECONDITION(minLength <= maxLength, "bad lengths provided");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }
  auto *res = new SparseIntVect<std::int32_t>(nBits);
  const double *dm;
  try {
    if (use2D) {
      dm = MolOps::getDistanceMat(*lmol);
    } else {
      dm = MolOps::get3DDistanceMat(*lmol, confId);
    }
  } catch (const ConformerException &) {
    delete res;
    throw;
  }

  const unsigned int nAtoms = lmol->getNumAtoms();

  std::vector<std::uint32_t> atomCodes;
  atomCodes.reserve(nAtoms);
  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    } else {
      atomCodes.push_back((*atomInvariants)[(*atomItI)->getIdx()]);
    }
  }

  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    unsigned int i = (*atomItI)->getIdx();
    if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(), i) !=
                           ignoreAtoms->end()) {
      continue;
    }
    if (!fromAtoms) {
      for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
           atomItJ != lmol->endAtoms(); ++atomItJ) {
        unsigned int j = (*atomItJ)->getIdx();
        if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                     j) != ignoreAtoms->end()) {
          continue;
        }
        unsigned int dist =
            static_cast<unsigned int>(floor(dm[i * nAtoms + j]));
        if (dist >= minLength && dist <= maxLength) {
          std::uint32_t bit = 0;
          gboost::hash_combine(bit, std::min(atomCodes[i], atomCodes[j]));
          gboost::hash_combine(bit, dist);
          gboost::hash_combine(bit, std::max(atomCodes[i], atomCodes[j]));
          updateElement(*res, static_cast<std::int32_t>(bit % nBits));
        }
      }
    } else {
      BOOST_FOREACH (std::uint32_t j, *fromAtoms) {
        if (j != i) {
          if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                       j) != ignoreAtoms->end()) {
            continue;
          }
          unsigned int dist =
              static_cast<unsigned int>(floor(dm[i * nAtoms + j]));
          if (dist >= minLength && dist <= maxLength) {
            std::uint32_t bit = 0;
            gboost::hash_combine(bit, std::min(atomCodes[i], atomCodes[j]));
            gboost::hash_combine(bit, dist);
            gboost::hash_combine(bit, std::max(atomCodes[i], atomCodes[j]));
            updateElement(*res, static_cast<std::int32_t>(bit % nBits));
          }
        }
      }
    }
  }
  return res;
}



namespace {
template <typename T>
void TorsionFpCalc(T *res, const ROMol &mol, unsigned int nBits,
                   unsigned int targetSize,
                   const std::vector<std::uint32_t> *fromAtoms,
                   const std::vector<std::uint32_t> *ignoreAtoms,
                   const std::vector<std::uint32_t> *atomInvariants,
                   bool includeChirality) {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }
  std::vector<std::uint32_t> atomCodes;
  atomCodes.reserve(lmol->getNumAtoms());
  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    } else {
      // need to add to the atomCode here because we subtract off up to 2 below
      // as part of the branch correction
      atomCodes.push_back(((*atomInvariants)[(*atomItI)->getIdx()] << 1) + 1);
    }
  }

  boost::dynamic_bitset<> *fromAtomsBV = nullptr;
  if (fromAtoms) {
    fromAtomsBV = new boost::dynamic_bitset<>(lmol->getNumAtoms());
    BOOST_FOREACH (std::uint32_t fAt, *fromAtoms) { fromAtomsBV->set(fAt); }
  }
  boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
  if (ignoreAtoms) {
    ignoreAtomsBV = new boost::dynamic_bitset<>(lmol->getNumAtoms());
    BOOST_FOREACH (std::uint32_t fAt, *ignoreAtoms) {
      ignoreAtomsBV->set(fAt);
    }
  }

  PATH_LIST paths = findAllPathsOfLengthN(*lmol, targetSize, false);
  for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
       ++pathIt) {
    bool keepIt = true;
    if (fromAtomsBV) {
      keepIt = false;
    }
    const PATH_TYPE &path = *pathIt;
    if (fromAtomsBV) {
      if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
          fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
        keepIt = true;
      }
    }
    if (keepIt && ignoreAtomsBV) {
      BOOST_FOREACH (int pElem, path) {
        if (ignoreAtomsBV->test(pElem)) {
          keepIt = false;
          break;
        }
      }
    }
    if (keepIt) {
      std::vector<std::uint32_t> pathCodes(targetSize);
      for (unsigned int i = 0; i < targetSize; ++i) {
        unsigned int code = atomCodes[path[i]] - 1;
        // subtract off the branching number:
        if (i > 0 && i < targetSize - 1) {
          --code;
        }
        pathCodes[i] = code;
      }
      size_t bit = getTopologicalTorsionHash(pathCodes);
      updateElement(*res, bit % nBits);
    }
  }
  delete fromAtomsBV;
  delete ignoreAtomsBV;
}
}  // namespace


}  // end of namespace ISIDA 
}  // end of namespace RDKit
