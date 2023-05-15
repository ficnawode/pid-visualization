void momentum() {
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  const int NEvents = treeIn->GetEntries();

  // declare branches and hook them up to the session
  auto *rec_header = new AnalysisTree::EventHeader();
  treeIn->SetBranchAddress("RecEventHeader.", &rec_header);
  auto *sim_tracks = new AnalysisTree::Particles();
  treeIn->SetBranchAddress("SimParticles.", &sim_tracks);
  auto *vtx_tracks = new AnalysisTree::TrackDetector();
  treeIn->SetBranchAddress("VtxTracks.", &vtx_tracks);
  auto *trd_tracks = new AnalysisTree::TrackDetector();
  treeIn->SetBranchAddress("TrdTracks.", &trd_tracks);
  auto *tof_hits = new AnalysisTree::HitDetector();
  treeIn->SetBranchAddress("TofHits.", &tof_hits);
  auto *vtx_tof_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2TofHits.", &vtx_tof_matching);
  auto *tof_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("TofHits2SimParticles.", &tof_sim_matching);
  auto *vtx_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2SimParticles.", &vtx_sim_matching);

  TFile *fileOut = TFile::Open("momentum.root", "recreate");

  // declare fields to be accessed and get their id
  AnalysisTree::Configuration *treeConfig = treeIn->GetConfiguration();
  // SIMULATED
  const int sp = treeConfig->GetBranchConfig("SimParticles").GetFieldId("p");
  const int sx = treeConfig->GetBranchConfig("SimParticles").GetFieldId("x");
  const int sy = treeConfig->GetBranchConfig("SimParticles").GetFieldId("y");
  const int sz = treeConfig->GetBranchConfig("SimParticles").GetFieldId("z");
  const int smother_id =
      treeConfig->GetBranchConfig("SimParticles").GetFieldId("mother_id");

  // STS+MVD
  const int rp = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("p");
  const int rpT = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("pT");
  const int rEta = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("eta");

  // TOF
  const int mass2 = treeConfig->GetBranchConfig("TofHits").GetFieldId("mass2");
  const int qp_tof =
      treeConfig->GetBranchConfig("TofHits").GetFieldId("qp_tof");
  const int x_tof = treeConfig->GetBranchConfig("TofHits").GetFieldId("x");
  const int y_tof = treeConfig->GetBranchConfig("TofHits").GetFieldId("y");

  // TRD
  const int trd_p = treeConfig->GetBranchConfig("TrdTracks").GetFieldId("p");
  const int trd_eloss0 =
      treeConfig->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_0");
  const int trd_eloss1 =
      treeConfig->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_1");
  const int trd_eloss2 =
      treeConfig->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_2");
  const int trd_eloss3 =
      treeConfig->GetBranchConfig("TrdTracks").GetFieldId("energy_loss_3");

  // declare histograms
  TH2F hc_qp_mass2("hc_qp_mass2",
                   "correlation qp_tof mass2; sign(q)*p (GeV/c);mass^2 (GeV)^2",
                   500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_protons(
      "hc_qp_mass2 protons",
      "correlation qp_tof mass2 protons; sign(q)*p (GeV/c);mass^2 (GeV)^2", 500,
      -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_pion_plus(
      "hc_qp_mass2 pi+",
      "correlation qp_tof mass2 pi +; sign(q)*p (GeV/c);mass^2 (GeV)^2", 500,
      -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_pion_minus(
      "hc_qp_mass2 pi-",
      "correlation qp_tof mass2 pi-; sign(q)*p (GeV/c);mass^2 (GeV)^2", 500,
      -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_kaon_plus(
      "hc_qp_mass2 K+",
      "correlation qp_tof mass2 K+; sign(q)*p (GeV/c);mass^2 (GeV)^2", 500, -16,
      16, 500, -5, 10);
  TH2F hc_qp_mass2_kaon_minus(
      "hc_qp_mass2 K-",
      "correlation qp_tof mass2 K-; sign(q)*p (GeV/c);mass^2 (GeV)^2", 500, -16,
      16, 500, -5, 10);
  TH2F hc_qp_mass2_others("hc_qp_mass2 others",
                          "correlation qp_tof mass2 other particles; "
                          "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                          500, -16, 16, 500, -5, 10);

  TH2F p_mass2("p_mass2", "Matching STS p TOF mass2 ;p ;mass^2", 400, -2, 20,
               400, -2, 10);
  TH2F p_mass2_protons("p_mass2 protons",
                       "Matching STS p TOF mass2 protons;p ;mass^2", 400, -2,
                       20, 400, -2, 10);
  TH2F p_mass2_pion_plus("p_mass2 pi+",
                         "Matching STS p TOF mass2 pi+;p ;mass^2", 400, -2, 20,
                         400, -2, 10);
  TH2F p_mass2_pion_minus("p_mass2 pi-",
                          "Matching STS p TOF mass2 pi-;p ;mass^2", 400, -2, 20,
                          400, -2, 10);
  TH2F p_mass2_kaon_plus("p_mass2 K+", "Matching STS p TOF mass2 K+;p ;mass^2",
                         400, -2, 20, 400, -2, 10);
  TH2F p_mass2_kaon_minus("p_mass2 K-", "Matching STS p TOF mass2 K-;p ;mass^2",
                          400, -2, 20, 400, -2, 10);
  TH2F p_mass2_others("p_mass2 others",
                      "Matching STS p TOF mass2 other particles;p ;mass^2", 400,
                      -2, 20, 400, -2, 10);

  TH2F pT_eta("pT_eta", "STS pT_eta;pT ;eta", 400, -1, 14, 400, -1, 8);
  TH2F pT_eta_protons("pT_eta protons", "STS pT_eta protons;pT ;eta", 400, -1,
                      14, 400, -1, 8);
  TH2F pT_eta_kaon_plus("pT_eta K+", "STS pT_eta K+;pT ;eta", 400, -1, 14, 400,
                        -1, 8);
  TH2F pT_eta_kaon_minus("pT_eta K-", "STS pT_eta K-;pT ;eta", 400, -1, 14, 400,
                         -1, 8);
  TH2F pT_eta_pion_plus("pT_eta pi+", "STS pT_eta pi+;pT ;eta", 400, -1, 14,
                        400, -1, 8);
  TH2F pT_eta_pion_minus("pT_eta pi-", "STS pT_eta pi-;pT ;eta", 400, -1, 14,
                         400, -1, 8);
  TH2F pT_eta_others("pT_eta others", "STS pT_eta others;pT ;eta", 400, -1, 14,
                     400, -1, 8);

  TH2F eloss0_p("eloss0_p", " TRD deltaE_0 p ;p ;eloss0", 400, -2, 20, 400, -4,
                100);
  TH2F eloss1_p("eloss1_p", " TRD deltaE_1 p ;p ;eloss1", 400, -2, 20, 400, -4,
                100);
  TH2F eloss2_p("eloss2_p", " TRD deltaE_2 p ;p ;eloss2", 400, -2, 20, 400, -4,
                100);
  TH2F eloss3_p("eloss3_p", " TRD deltaE_3 p ;p ;eloss3", 400, -2, 20, 400, -4,
                100);
  TH2F tof_x_y("TOF x y", "TOF x y; x; y", 400, -600, 600, 400, -400, 400);

  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

    // tof hits
    for (const auto &tof_hit : *tof_hits) {
      float tof_mass2 = tof_hit.GetField<float>(mass2);
      const float tof_qp_tof = tof_hit.GetField<float>(qp_tof);
      const float tof_x = tof_hit.GetField<float>(x_tof);
      const float tof_y = tof_hit.GetField<float>(y_tof);

      hc_qp_mass2.Fill(tof_qp_tof, tof_mass2);
      tof_x_y.Fill(tof_x, tof_y);

      const int tof_id = vtx_tof_matching->GetMatch(tof_hit.GetId());
      if (tof_id > 0) {
        const float vtx_p = vtx_tracks->GetChannel(tof_id).GetField<float>(rp);
        p_mass2.Fill(vtx_p, tof_mass2);
      }

      const int sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      if (sim_id > 0) {
        int particle_id = sim_tracks->GetChannel(sim_id).GetPid();
        switch (particle_id) {
        case 2212:
          hc_qp_mass2_protons.Fill(tof_qp_tof, tof_mass2);
          break;
        case 321:
          hc_qp_mass2_kaon_plus.Fill(tof_qp_tof, tof_mass2);
          break;
        case -321:
          hc_qp_mass2_kaon_minus.Fill(tof_qp_tof, tof_mass2);
          break;
        case 211:
          hc_qp_mass2_pion_plus.Fill(tof_qp_tof, tof_mass2);
          break;
        case -211:
          hc_qp_mass2_pion_minus.Fill(tof_qp_tof, tof_mass2);
          break;
        default:
          hc_qp_mass2_others.Fill(tof_qp_tof, tof_mass2);
        }
      }

      if (sim_id > 0 && tof_id > 0) {
        const float vtx_p = vtx_tracks->GetChannel(tof_id).GetField<float>(rp);
        int particle_id = sim_tracks->GetChannel(sim_id).GetPid();
        switch (particle_id) {
        case 2212:
          p_mass2_protons.Fill(vtx_p, tof_mass2);
          break;
        case 321:
          p_mass2_kaon_plus.Fill(vtx_p, tof_mass2);
          break;
        case -321:
          p_mass2_kaon_minus.Fill(vtx_p, tof_mass2);
          break;
        case 211:
          p_mass2_pion_plus.Fill(vtx_p, tof_mass2);
          break;
        case -211:
          p_mass2_pion_minus.Fill(vtx_p, tof_mass2);
          break;
        default:
          p_mass2_others.Fill(vtx_p, tof_mass2);
        }
      }
    }
    for (const auto &track : *vtx_tracks) {
      const float _vtx_eta = track.GetField<float>(rEta);
      const float _vtx_pT = track.GetField<float>(rpT);
      pT_eta.Fill(_vtx_eta, _vtx_pT);

      const int sim_id = vtx_sim_matching->GetMatch(track.GetId());
      if (sim_id > 0) {
        int particle_id = sim_tracks->GetChannel(sim_id).GetPid();
        switch (particle_id) {
        case 2212:
          pT_eta_protons.Fill(_vtx_eta, _vtx_pT);
          break;
        case 321:
          pT_eta_kaon_plus.Fill(_vtx_eta, _vtx_pT);
          break;
        case -321:
          pT_eta_kaon_minus.Fill(_vtx_eta, _vtx_pT);
          break;
        case 211:
          pT_eta_pion_plus.Fill(_vtx_eta, _vtx_pT);
          break;
        case -211:
          pT_eta_pion_minus.Fill(_vtx_eta, _vtx_pT);
          break;
        default:
          pT_eta_others.Fill(_vtx_eta, _vtx_pT);
        }
      }
    }

    for (const auto &track : *trd_tracks) {
      const float _trd_p = track.GetField<float>(trd_p);
      const float _trd_eloss0 = track.GetField<float>(trd_eloss0);
      const float _trd_eloss1 = track.GetField<float>(trd_eloss1);
      const float _trd_eloss2 = track.GetField<float>(trd_eloss2);
      const float _trd_eloss3 = track.GetField<float>(trd_eloss3);

      eloss0_p.Fill(_trd_p, _trd_eloss0);
      eloss1_p.Fill(_trd_p, _trd_eloss1);
      eloss2_p.Fill(_trd_p, _trd_eloss2);
      eloss3_p.Fill(_trd_p, _trd_eloss3);
    }
  }

  // write to histograms
  hc_qp_mass2.Write();
  hc_qp_mass2_protons.Write();
  hc_qp_mass2_pion_plus.Write();
  hc_qp_mass2_pion_minus.Write();
  hc_qp_mass2_kaon_plus.Write();
  hc_qp_mass2_kaon_minus.Write();
  hc_qp_mass2_others.Write();

  p_mass2.Write();
  p_mass2_protons.Write();
  p_mass2_pion_plus.Write();
  p_mass2_pion_minus.Write();
  p_mass2_kaon_plus.Write();
  p_mass2_kaon_minus.Write();
  p_mass2_others.Write();

  pT_eta.Write();
  pT_eta_protons.Write();
  pT_eta_kaon_plus.Write();
  pT_eta_kaon_minus.Write();
  pT_eta_pion_plus.Write();
  pT_eta_pion_minus.Write();
  pT_eta_others.Write();

  eloss0_p.Write();
  eloss1_p.Write();
  eloss2_p.Write();
  eloss3_p.Write();

  tof_x_y.Write();

  fileOut->Close();
}