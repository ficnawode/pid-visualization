
void qp_beta() {
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  TFile *fileOut = TFile::Open("out/qp_beta.root", "recreate");
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
  const int qp_tof =
      treeConfig->GetBranchConfig("TofHits").GetFieldId("qp_tof");

  // declare histograms
  TH2F qp_beta("qp_beta", "correlation qp_tof #beta; sign(q)*p (GeV/c);#beta ",
               500, -16, 16, 500, 0, 3);
  TH2F qp_beta_protons("qp_beta protons sim pid",
                       "correlation qp_tof #beta protons sim pid; "
                       "sign(q)*p (GeV/c);#beta ",
                       500, -16, 16, 500, 0, 3);
  TH2F qp_beta_pion_plus(
      "qp_beta pi+ sim pid",
      "correlation qp_tof #beta pi + sim pid; sign(q)*p (GeV/c);#beta ", 500,
      -16, 16, 500, 0, 3);
  TH2F qp_beta_pion_minus(
      "qp_beta pi- sim pid",
      "correlation qp_tof #beta pi- sim pid; sign(q)*p (GeV/c);#beta ", 500,
      -16, 16, 500, 0, 3);
  TH2F qp_beta_kaon_plus(
      "qp_beta K+ sim pid",
      "correlation qp_tof #beta K+ sim pid; sign(q)*p (GeV/c);#beta ", 500, -16,
      16, 500, 0, 3);
  TH2F qp_beta_kaon_minus(
      "qp_beta K- sim pid",
      "correlation qp_tof #beta K- sim pid; sign(q)*p (GeV/c);#beta ", 500, -16,
      16, 500, 0, 3);
  TH2F qp_beta_others("qp_beta others sim pid",
                      "correlation qp_tof #beta other particles sim pid; "
                      "sign(q)*p (GeV/c);#beta ",
                      500, -16, 16, 500, 0, 3);

  // fill histograms
  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

    for (const auto &tof_hit : *tof_hits) {
      const float tof_mass2 = tof_hit.GetField<float>(mass2);
      const float tof_qp_tof = tof_hit.GetField<float>(qp_tof);

      const int tof_id = vtx_tof_matching->GetMatch(tof_hit.GetId());
      if (tof_id < 0)
        continue;

      const float vtx_p = vtx_tracks->GetChannel(tof_id).GetField<float>(rp);
      const float vtx_p2 = vtx_p * vtx_p;
      const float beta = TMath::Sqrt(1 / (1 + (tof_mass2 / vtx_p2)));

      qp_beta.Fill(tof_qp_tof, beta);

      const int sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      if (sim_id > 0) {
        int sim_pdg = sim_tracks->GetChannel(sim_id).GetPid();
        switch (sim_pdg) {
        case 2212:
          qp_beta_protons.Fill(tof_qp_tof, beta);
          break;
        case 321:
          qp_beta_kaon_plus.Fill(tof_qp_tof, beta);
          break;
        case -321:
          qp_beta_kaon_minus.Fill(tof_qp_tof, beta);
          break;
        case 211:
          qp_beta_pion_plus.Fill(tof_qp_tof, beta);
          break;
        case -211:
          qp_beta_pion_minus.Fill(tof_qp_tof, beta);
          break;
        default:
          qp_beta_others.Fill(tof_qp_tof, beta);
        }
      }
    }
  }

  // write to histograms
  qp_beta.Write();
  qp_beta_protons.Write();
  qp_beta_kaon_plus.Write();
  qp_beta_kaon_minus.Write();
  qp_beta_pion_plus.Write();
  qp_beta_pion_minus.Write();
  qp_beta_others.Write();

  fileOut->Close();
}