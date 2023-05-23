
void p_mass2() {
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  TFile *fileOut = TFile::Open("out/p_mass2.root", "recreate");
  const int NEvents = treeIn->GetEntries();

  // declare branches and hook them up to the session
  auto *sim_tracks = new AnalysisTree::Particles();
  treeIn->SetBranchAddress("SimParticles.", &sim_tracks);
  auto *vtx_tracks = new AnalysisTree::TrackDetector();
  treeIn->SetBranchAddress("VtxTracks.", &vtx_tracks);
  auto *tof_hits = new AnalysisTree::HitDetector();
  treeIn->SetBranchAddress("TofHits.", &tof_hits);
  auto *vtx_tof_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2TofHits.", &vtx_tof_matching);
  auto *tof_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("TofHits2SimParticles.", &tof_sim_matching);
  auto *vtx_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2SimParticles.", &vtx_sim_matching);

  // declare fields to be accessed and get their id
  AnalysisTree::Configuration *treeConfig = treeIn->GetConfiguration();

  // STS+MVD
  const int rp = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("p");

  // TOF
  const int mass2 = treeConfig->GetBranchConfig("TofHits").GetFieldId("mass2");

  // declare histograms
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

  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

    // tof hits
    for (const auto &tof_hit : *tof_hits) {
      float tof_mass2 = tof_hit.GetField<float>(mass2);

      const int tof_id = vtx_tof_matching->GetMatch(tof_hit.GetId());
      if (tof_id < 0)
        continue;

      const float vtx_p = vtx_tracks->GetChannel(tof_id).GetField<float>(rp);
      p_mass2.Fill(vtx_p, tof_mass2);

      const int sim_id_tof = tof_sim_matching->GetMatch(tof_hit.GetId());
      const int sim_id_vtx = vtx_sim_matching->GetMatch(tof_hit.GetId());
      // TODO!!!!
      if (sim_id_tof != sim_id_vtx)
        continue;

      if (sim_id < 0)
        continue;

      int particle_id = sim_tracks->GetChannel(sim_id).GetPid();
      switch (particle_id) {
      case 2212://protons
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

  // write to histograms
  p_mass2.Write();
  p_mass2_protons.Write();
  p_mass2_pion_plus.Write();
  p_mass2_pion_minus.Write();
  p_mass2_kaon_plus.Write();
  p_mass2_kaon_minus.Write();
  p_mass2_others.Write();

  fileOut->Close();
}
