#include "PhysicalMixingAuxVars.h"
#include "Field.h"
#include <limits>

namespace OMEGA {

PhysicalMixingAuxVars::PhysicalMixingAuxVars(const std::string &AuxStateSuffix,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord,
                                             const I4 NTracers)
    : PhysicalMixingCell("PhysicalMixingCell" + AuxStateSuffix, NTracers,
                         Mesh->NCellsSize, VCoord->NVertLayers),
      GeomZMid(VCoord->GeomZMid), // TODO(confirm): verify GeomZMid actually
                                  // lives on VertCoord and is populated by
                                  // the time this runs (computeGeomZHeight
                                  // is called in computeMomVertAux, which
                                  // computeAll invokes before this would --
                                  // confirm ordering once wired in).
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

void PhysicalMixingAuxVars::registerFields(const std::string &AuxGroupName,
                                           const std::string &MeshName) const {
   int NDims = 3;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix = (MeshName == "Default") ? "" : MeshName;

   DimNames[0] = "NTracers";
   DimNames[1] = "NCells" + DimSuffix;
   DimNames[2] = "NVertLayers";

   auto PhysicalMixingCellField = Field::create(
       PhysicalMixingCell.label(),
       "physical (turbulent) tracer variance decay, 2*Kv*(dC/dz)^2, "
       "diagnostic only -- Burchard and Rennau (2008)",
       "",                                // units -- TODO confirm convention
       "",                                // CF standard name
       std::numeric_limits<Real>::min(),  // min valid value
       std::numeric_limits<Real>::max(),  // max valid value
       NDims,                             // number of dimensions
       DimNames                           // dimension names
   );

   // NOTE: registered into AuxGroupName as passed in, which the caller
   // (AuxiliaryState.cpp) should set to a dedicated "DVDDiagnostics" group
   // rather than the general AuxiliaryState group -- see Carolyn's note
   // about keeping all DVD-related fields together.
   FieldGroup::addFieldToGroup(PhysicalMixingCell.label(), AuxGroupName);

   PhysicalMixingCellField->attachData<Array3DReal>(PhysicalMixingCell);
}

void PhysicalMixingAuxVars::unregisterFields() const {
   Field::destroy(PhysicalMixingCell.label());
}

} // namespace OMEGA
