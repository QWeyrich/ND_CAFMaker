#include "NDLArTMSUniqueMatchRecoFiller.h"
#include <cmath>

namespace cafmaker
{
  NDLArTMSUniqueMatchRecoFiller::NDLArTMSUniqueMatchRecoFiller(double sigmaX, double sigmaY, bool singleAngle, double sigmaTh, double sigmaThX, double sigmaThY, bool useTime, double sigmaT, double fCut)
    : IRecoBranchFiller("LArTMSMatcher")
  {
    sigma_x = sigmaX;
    sigma_y = sigmaY;
    single_angle = singleAngle;
    sigma_angle = sigmaTh;
    sigma_angle_x = sigmaThX;
    sigma_angle_y = sigmaThY;
    use_time = useTime;
    sigma_t = sigmaT;
    fcut = fCut;
    // nothing to do
    SetConfigured(true);
  }

  std::vector<double> NDLArTMSUniqueMatchRecoFiller::Project_track(caf::SRTrack track, bool forward) const
  {
    if (forward) { // projects a LAr track forward to TMS
      double x = track.end.x;
      double y = track.end.y;
      double z = track.end.z;

      double dir_x = track.enddir.x;
      double dir_y = track.enddir.y;
      double dir_z = track.enddir.z;

      double proj_z = TMSLims.tms_z_lim1 - z;
      double proj_x = dir_x*proj_z/dir_z + x;
      double proj_y = dir_y*proj_z/dir_z + y;
    }
    else { // projects a TMS track backward to LAr
      double x = track.start.x;
      double y = track.start.y;
      double z = track.start.z;

      double dir_x = track.dir.x;
      double dir_y = track.dir.y;
      double dir_z = track.dir.z;

      double proj_z = z - LArLims.lar_z_lim2;
      double proj_x = -dir_x*proj_z/dir_z + x;
      double proj_y = -dir_y*proj_z/dir_z + y;
    }
    std::vector<double> proj_point = {proj_x, proj_y, proj_z};
    return proj_point;
  }

  std::vector<double> NDLArTMSUniqueMatchRecoFiller::Angle_between_tracks(const caf::SRTrack tms_track, const caf::SRTrack lar_track) const
  {
      double tms_dir_x = tms_track.dir.x;
      double tms_dir_y = tms_track.dir.y;
      double tms_dir_z = tms_track.dir.z;

      double lar_dir_x = lar_track.dir.x;
      double lar_dir_y = lar_track.dir.y;
      double lar_dir_z = lar_track.dir.z;

      double xz_dot_prod = tms_dir_x*lar_dir_x + tms_dir_z*lar_dir_z;
      if (xz_dot_prod != 0) {
        double xz_dot_prod = xz_dot_prod/(sqrt(pow(tms_dir_x,2)+pow(tms_dir_z,2))*sqrt(pow(lar_dir_x,2)+pow(lar_dir_z,2)));
      }
      double yz_dot_prod = tms_dir_y*lar_dir_y + tms_dir_z*lar_dir_z;
      if (yz_dot_prod != 0) {
        double yz_dot_prod = yz_dot_prod/(sqrt(pow(tms_dir_y,2)+pow(tms_dir_z,2))*sqrt(pow(lar_dir_y,2)+pow(lar_dir_z,2)));
      }
      double dot_prod = tms_dir_x*lar_dir_x + tms_dir_y*lar_dir_y + tms_dir_z*lar_dir_z;
      double angle_x = 180.0/std::numbers::pi * acos(dot_prod_x);
      double angle_y = 180.0/std::numbers::pi * acos(dot_prod_y);
      double angle_overall = 180.0/std::numbers::pi * acos(dot_prod);
      std::vector<double> angles = {angle_x,angle_y,angle_overall};
      return angles;
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_TMS_track(const caf::SRTrack tms_track, const double tms_z_cutoff) const
  {
    double x_start = tms_track.start.x;
    double y_start = tms_track.start.y;
    double z_start = tms_track.start.z;

    double x_end = tms_track.end.x;
    double y_end = tms_track.end.y;
    double z_end = tms_track.end.z;

    double dir_x = tms_track.dir.x;
    double dir_y = tms_track.dir.y;
    double dir_z = tms_track.dir.z;

    if ((x_start > TMSLims.tms_x_lim1)&&(x_start < TMSLims.tms_x_lim2) &&
        (y_start > TMSLims.tms_y_lim1)&&(y_start < TMSLims.tms_y_lim2) &&
        (z_start > TMSLims.tms_z_lim1)&&(z_start < TMSLims.tms_z_lim1 + tms_z_cutoff) && // checks track begins within fiducial volume and close enough to front
      
        (Project_track(tms_track,false)[0] > LArLims.lar_x_lim1)&&(Project_track(tms_track,false)[0] < LArLims.lar_x_lim2) &&
        (Project_track(tms_track,false)[1] > LArLims.lar_y_lim1)&&(Project_track(tms_track,false)[1] < LArLims.lar_y_lim2)) // checks that direction would have allowed it to originate from LAr
          {return true} 
    else {return false}
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_LAr_track(const caf::SRTrack lar_track, const double lar_z_cutoff) const
  {
    double x_start = lar_track.start.x;
    double y_start = lar_track.start.y;
    double z_start = lar_track.start.z;

    double x_end = lar_track.end.x;
    double y_end = lar_track.end.y;
    double z_end = lar_track.end.z;

    double dir_x = lar_track.enddir.x;
    double dir_y = lar_track.enddir.y;
    double dir_z = lar_track.enddir.z;

    if ((x_start > LArLims.lar_x_lim1)&&(x_start < LArLims.lar_x_lim2) &&
        (y_start > LArLims.lar_y_lim1)&&(y_start < LArLims.lar_y_lim2) &&
        (z_start > LArLims.lar_z_lim1)&&(z_start < LArLims.lar_z_lim2) && // checks track begins within fiducial volume

        (x_end > LArLims.lar_x_lim1)&&(x_end < LArLims.lar_x_lim2) &&
        (y_end > LArLims.lar_y_lim1)&&(y_end < LArLims.lar_y_lim2) &&
        (z_end > LArLims.lar_z_lim2 - lar_z_cutoff)&&(z_end < LArLims.lar_z_lim2) && // checks track ends close enough to back of LAr
      
        (Project_track(lar_track,true)[0] > TMSLims.tms_x_lim1)&&(Project_track(lar_track,true)[0] < TMSLims.tms_x_lim2) &&
        (Project_track(lar_track,true)[1] > TMSLims.tms_y_lim1)&&(Project_track(lar_track,true)[1] < TMSLims.tms_y_lim2)) // checks that direction would allow it to hit TMS
        {return true} 
    else {return false}
  }

  bool NDLArTMSUniqueMatchRecoFiller::Track_match_sorter(const caf::SRNDTrackAssn trackMatch1, const caf::SRNDTrackAssn trackMatch2) const {
    double fScore1 = trackMatch1.matchScore;
    double fScore2 = trackMatch2.matchScore;

    if (fScore1 < fScore2) {return true}
    else {return false}
  }

  void
  NDLArTMSUniqueMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const
  {
    std::vector<caf::SRNDTrackAssn> possiblePandoraMatches; // vector will store all possible matched tracks between Pandora and TMS
    std::vector<caf::SRNDTrackAssn> possibleSPINEMatches; // vector will store all possible matched tracks between SPINE and TMS

    double tms_z_cutoff = 10;
    double lar_z_cutoff = 10; // tracks must overlap last/first 10 cm of the detectors

    for (unsigned int ixn_tms = 0; ixn_tms < sr.nd.tms.nixn; ixn_tms++)
    {
      caf::SRTMSInt tms_int = sr.nd.tms.ixn[ixn_tms];
      unsigned int n_tms_tracks = tms_int.ntracks;

      for (unsigned int itms = 0; itms < n_tms_tracks; itms++)
      {
        caf::SRTrack tms_trk = tms_int.tracks[itms];

        if (!Consider_TMS_track(tms_trk,tms_z_cutoff)) {
          continue; //skips the tms track if it isn't suitable according to the function
        }

        for (unsigned int ixn_pan = 0; ixn_pan < sr.nd.lar.pandora.npandora; ixn_pan++)
        {
          caf::SRNDLArInt pan_int = sr.nd.lar.pandora.ixn[ixn_pan];
          unsigned int n_pan_tracks = pan_int.ntracks;

          for (unsigned int ipan = 0; ipan < n_pan_tracks; ipan++)
          {
            caf::SRTrack pan_trk = pan_int.tracks[ipan];

            if (!Consider_LAr_track(pan_trk,lar_z_cutoff)) {
              continue; //skips the lar track if it isn't suitable according to the function
            }

            std::vector<double> proj_vec = Project_track(pan_trk,true);

            double delta_x = tms_trk.start.x - proj_vec[0];
            double delta_y = tms_trk.start.y - proj_vec[1];

            std::vector<double> angles = Angle_between_tracks(tms_trk,pan_trk);

            if (single_angle) {
              double angle = angles.end();
              double fScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle/sigma_angle,2);
            }
            else {
              double angle_x = angles[0];
              double angle_y = angles[1];
              double fScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle_x/sigma_angle_x,2)+ pow(angle_y/sigma_angle_y,2);
            }
            if (use_time) {
              double delta_t = 0; // placeholder until time is added
              fScore += pow(delta_t/sigma_t,2)
            }
        
            caf::SRTMSID tmsid;
            tmsid.ixn = ixn_tms;
            tmsid.idx = itms;
            caf::SRNDLArID panid;
            panid.ixn = ixn_pan;
            panid.idx = ipan;

            caf::SRNDTrackAssn potential_match;
            potential_match.tmsid = tmsid;
            potential_match.larid = panid;
            potential_match.matchScore = fScore;
            potential_match.transdispl = sqrt(pow(delta_x,2)+pow(delta_y,2));
            potential_match.angdispl = cos(std::numbers::pi/180.0 * angles.end())

            possiblePandoraMatches.push_back(potential_match)
          }
        }

        for (unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.dlp.ndlp; ixn_dlp++)
        {
          caf::SRNDLArInt dlp_int = sr.nd.lar.dlp.ixn[ixn_dlp];
          unsigned int n_dlp_tracks = dlp_int.ntracks;

          for (unsigned int idlp = 0; idlp < n_dlp_tracks; idlp++)
          {
            caf::SRTrack dlp_trk = dlp_int.tracks[idlp];

            if (!Consider_LAr_track(dlp_trk,lar_z_cutoff)) {
              continue; //skips the lar track if it isn't suitable according to the function
            }

            std::vector<double> proj_vec = Project_track(dlp_trk,true);

            double delta_x = tms_trk.start.x - proj_vec[0];
            double delta_y = tms_trk.start.y - proj_vec[1];

            std::vector<double> angles = Angle_between_tracks(tms_trk,dlp_trk);

            if (single_angle) {
              double angle = angles.end();
              double fScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle/sigma_angle,2);
            }
            else {
              double angle_x = angles[0];
              double angle_y = angles[1];
              double fScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle_x/sigma_angle_x,2)+ pow(angle_y/sigma_angle_y,2);
            }
            if (use_time) {
              double delta_t = 0; // placeholder until time is added
              fScore += pow(delta_t/sigma_t,2)
            }
        
            caf::SRTMSID tmsid;
            tmsid.ixn = ixn_tms;
            tmsid.idx = itms;
            caf::SRNDLArID dlpid;
            dlpid.ixn = ixn_dlp;
            dlpid.idx = idlp;

            caf::SRNDTrackAssn potential_match;
            potential_match.tmsid = tmsid;
            potential_match.larid = dlpid;
            potential_match.matchScore = fScore;
            potential_match.transdispl = sqrt(pow(delta_x,2)+pow(delta_y,2));
            potential_match.angdispl = cos(std::numbers::pi/180.0 * angles.end());

            possibleSPINEMatches.push_back(potential_match);
          }
        }
      }
    }

    std::sort(possiblePandoraMatches.begin(),possiblePandoraMatches.end(),Track_match_sorter);

    std::vector<caf::SRNDLArID> matched_pan; 
    std::vector<caf::SRTMSID> matched_tms; // stores LAr (Pandora) and TMS indices that have already been matched

    for (unsigned int match_idx = 0; match_idx < possiblePandoraMatches.size(), match_idx++) {
      caf::SRNDTrackAssn track_match = possiblePandoraMatches[match_idx];
      double score = track_match.matchScore;
      if (score > f_cut) {break;}
      caf::SRNDLArID panid = track_match.larid;
      if (std::count(matched_pan.begin(),matched_pan.end(),panid) == 0) {
        caf::SRTMSID tmsid = track_match.tmsid;
        if (std::count(matched_tms.begin(),matched_tms.end(),tmsid) == 0) {
          matched_tms.push_back(tmsid);
          matched_pan.push_back(panid);
          sr.nd.trkmatch.push_back(track_match);
          sr.nd.ntrkmatch += 1;
        }
      }
    }

    std::sort(possibleSPINEMatches.begin(),possibleSPINEMatches.end(),Track_match_sorter);

    std::vector<caf::SRNDLArID> matched_dlp; 
    std::vector<caf::SRTMSID> matched_tms; // stores LAr (SPINE) and TMS indices that have already been matched

    for (unsigned int match_idx = 0; match_idx < possibleSPINEMatches.size(), match_idx++) {
      caf::SRNDTrackAssn track_match = possibleSPINEMatches[match_idx];
      double score = track_match.matchScore;
      if (score > f_cut) {break;}
      caf::SRNDLArID dlpid = track_match.dlpid;
      if (std::count(matched_dlp.begin(),matched_dlp.end(),dlpid) == 0) {
        caf::SRTMSID tmsid = track_match.tmsid;
        if (std::count(matched_tms.begin(),matched_tms.end(),tmsid) == 0) {
          matched_tms.push_back(tmsid);
          matched_dlp.push_back(dlpid);
          sr.nd.trkmatch.push_back(track_match);
          sr.nd.ntrkmatch += 1;
        }
      }
    }
  }

  // todo: this is a placeholder
  std::deque<Trigger> NDLArTMSUniqueMatchRecoFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    return std::deque<Trigger>();
  }

}