#include "NDLArTMSUniqueMatchRecoFiller.h"

namespace cafmaker
{
  NDLArTMSUniqueMatchRecoFiller::NDLArTMSUniqueMatchRecoFiller()
    : IRecoBranchFiller("LArTMSMatcher")
  {
    // nothing to do
    SetConfigured(true);
  }

  void
  NDLArTMSUniqueMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const
  {
    void
    
  }

  // todo: this is a placeholder
  std::deque<Trigger> NDLArTMSUniqueMatchRecoFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    return std::deque<Trigger>();
  }

}