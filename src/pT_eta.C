
void pT_eta() {
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  const int NEvents = treeIn->GetEntries();
  TFile *fileOut = TFile::Open("out/pT_eta.root", "recreate");

  auto *sim_tracks = new AnalysisTree::Particles();
  treeIn->SetBranchAddress("SimParticles.", &sim_tracks);
  auto *vtx_tracks = new AnalysisTree::TrackDetector();
  treeIn->SetBranchAddress("VtxTracks.", &vtx_tracks);
  auto *tof_hits = new AnalysisTree::HitDetector();
  auto *vtx_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2SimParticles.", &vtx_sim_matching);

  // declare fields to be accessed and get their id
  AnalysisTree::Configuration *treeConfig = treeIn->GetConfiguration();

  // STS+MVD
  const int rp = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("p");
  const int rpT = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("pT");
  const int rEta = treeConfig->GetBranchConfig("VtxTracks").GetFieldId("eta");

  // declare histograms
  TH2F pT_eta("pT_eta", "STS pT_eta;pT ;#eta", 400, -1, 9, 400, -1, 8);
  TH2F pT_eta_protons("pT_eta protons", "STS pT_eta protons;pT ;#eta", 400, -1,
                      9, 400, -1, 8);
  TH2F pT_eta_kaon_plus("pT_eta K+", "STS pT_eta K+;pT ;#eta", 400, -1, 9, 400,
                        -1, 8);
  TH2F pT_eta_kaon_minus("pT_eta K-", "STS pT_eta K-;pT ;#eta", 400, -1, 9, 400,
                         -1, 8);
  TH2F pT_eta_pion_plus("pT_eta pi+", "STS pT_eta pi+;pT ;#eta", 400, -1, 9,
                        400, -1, 8);
  TH2F pT_eta_pion_minus("pT_eta pi-", "STS pT_eta pi-;pT ;#eta", 400, -1, 9,
                         400, -1, 8);
  TH2F pT_eta_others("pT_eta others", "STS pT_eta others;pT ;#eta", 400, -1, 9,
                     400, -1, 8);

  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

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
  }

  // write to histograms
  pT_eta.Write();
  pT_eta_protons.Write();
  pT_eta_kaon_plus.Write();
  pT_eta_kaon_minus.Write();
  pT_eta_pion_plus.Write();
  pT_eta_pion_minus.Write();
  pT_eta_others.Write();

  fileOut->Close();
}