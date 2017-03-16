
#include "DataFormats/Common/interface/Wrapper.h"
#include "Math/Cartesian3D.h" 
#include "DataFormats/Math/interface/Vector3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "DataFormats/PCC/interface/PCC.h"
#include "DataFormats/Math/interface/Vector.h" 

namespace DataFormats_PCC {
  struct dictionary {
    reco::PCC b;
    edm::Wrapper<reco::PCC> b_w;
  };      
}
