/// Match ND-LAr and TMS tracks.
///
/// \author  Q. Weyrich <qweyrich@yorku.ca>
/// \date    July 2025

#ifndef ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H

#include "IRecoBranchFiller.h"
#include "MLNDLArRecoBranchFiller.h"
#include "TMSRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{
  class NDLArTMSUniqueMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArTMSUniqueMatchRecoFiller();

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::Matcher; }


    private:
      void MatchTracks(caf::StandardRecord &sr) const;

      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;
        
      std::vector<caf::SRNDTrackAssn> matchVectorSPINETMS;

      std::vector<caf::SRNDTrackAssn> matchVectorPandoraTMS;
  };
}

#endif //ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H