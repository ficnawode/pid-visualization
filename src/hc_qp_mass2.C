
void hc_qp_mass2() {
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  TFile *fileOut = TFile::Open("out/hc_qp_mass2.root", "recreate");
  const int NEvents = treeIn->GetEntries();

  // declare branches and hook them up to the session
  auto *sim_tracks = new AnalysisTree::Particles();
  treeIn->SetBranchAddress("SimParticles.", &sim_tracks);
  auto *tof_hits = new AnalysisTree::HitDetector();
  treeIn->SetBranchAddress("TofHits.", &tof_hits);
  auto *tof_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("TofHits2SimParticles.", &tof_sim_matching);

  // declare fields to be accessed and get their id
  AnalysisTree::Configuration *treeConfig = treeIn->GetConfiguration();
  // TOF
  const int mass2 = treeConfig->GetBranchConfig("TofHits").GetFieldId("mass2");
  const int qp_tof =
      treeConfig->GetBranchConfig("TofHits").GetFieldId("qp_tof");

  // declare histograms
  TH2F hc_qp_mass2("hc_qp_mass2",
                   "correlation qp_tof mass2; sign(q)*p (GeV/c);mass^2 (GeV)^2",
                   500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_protons("hc_qp_mass2 protons sim pid",
                           "correlation qp_tof mass2 protons sim pid; "
                           "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                           500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_pion_plus(
      "hc_qp_mass2 pi+ sim pid",
      "correlation qp_tof mass2 pi + sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_pion_minus(
      "hc_qp_mass2 pi- sim pid",
      "correlation qp_tof mass2 pi- sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_kaon_plus(
      "hc_qp_mass2 K+ sim pid",
      "correlation qp_tof mass2 K+ sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_kaon_minus(
      "hc_qp_mass2 K- sim pid",
      "correlation qp_tof mass2 K- sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      500, -16, 16, 500, -5, 10);
  TH2F hc_qp_mass2_others("hc_qp_mass2 others sim pid",
                          "correlation qp_tof mass2 other particles sim pid; "
                          "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                          500, -16, 16, 500, -5, 10);

  // fill histograms
  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

    for (const auto &tof_hit : *tof_hits) {
      float tof_mass2 = tof_hit.GetField<float>(mass2);
      const float tof_qp_tof = tof_hit.GetField<float>(qp_tof);

      hc_qp_mass2.Fill(tof_qp_tof, tof_mass2);

      const int sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      if (sim_id > 0) {
        int sim_pdg = sim_tracks->GetChannel(sim_id).GetPid();
        switch (sim_pdg) {
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
    }
  }

  // write to histograms
  hc_qp_mass2.Write();
  hc_qp_mass2_protons.Write();
  hc_qp_mass2_kaon_plus.Write();
  hc_qp_mass2_kaon_minus.Write();
  hc_qp_mass2_pion_plus.Write();
  hc_qp_mass2_pion_minus.Write();
  hc_qp_mass2_others.Write();

  fileOut->Close();
}