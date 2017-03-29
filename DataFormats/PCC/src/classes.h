#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/PCC/interface/PCC.h"

namespace DataFormats_PCC {
  struct dictionary {
    reco::PCC b;
    edm::Wrapper<reco::PCC> b_w;
  };      
}
