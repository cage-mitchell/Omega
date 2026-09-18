#ifndef OMEGA_AUX_PHYSICALMIXING_H
#define OMEGA_AUX_PHYSICALMIXING_H
#include "DataTypes.h"
#include "Field.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

/// Diagnostic physical mixing, per Burchard & Rennau (2008):
///   M_phy(l, i, k) = 2 * Kv(i, k) * (dC/dz)^2
/// computed per tracer, cell, and layer from the current tracer field and
/// the vertical diffusivity already computed by VertMix (Kv is passed in,
/// not owned here, to match the existing AuxVars pattern of receiving
/// dependencies as arguments -- see TracerAuxVars::computeVarsOnCells).
/// This is diagnostic only: it does not feed back into the prognostic
/// tracer tendencies.
class PhysicalMixingAuxVars {
 public:
   Array3DReal PhysicalMixingCell;

   PhysicalMixingAuxVars(const std::string &AuxStateSuffix,
                         const HorzMesh *Mesh, const VertCoord *VCoord,
                         const I4 NTracers);

   KOKKOS_FUNCTION void
   computeVarsOnCells(int L, int ICell, int KChunk,
                      const Array2DReal &VertDiff,
                      const Array3DReal &TrCell) const {
      const int KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));
      const int KBottom = MaxLayerCell(ICell);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K = KStart + KVec;
         Real DCDZ    = 0;
         // TODO(confirm): one-sided difference against the layer below,
         // matching the stencil TracerVertMixSetupOnCell uses for
         // VertDiff(ICell, K+1) in VertMix.h. Bottom-most layer left as
         // zero gradient (no layer below to difference against) -- revisit
         // once verified against MPAS-O's mpas_ocn_discrete_variance_decay.F
         // reference implementation.
         if (K < KBottom) {
            const Real DZ = GeomZMid(ICell, K) - GeomZMid(ICell, K + 1);
            DCDZ = (TrCell(L, ICell, K) - TrCell(L, ICell, K + 1)) / DZ;
         }
         PhysicalMixingCell(L, ICell, K) =
             2._Real * VertDiff(ICell, K) * DCDZ * DCDZ;
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array2DReal GeomZMid;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
