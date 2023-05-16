
void reco_sim_comp_tracks() {

  // TString fileName = "1.analysistree.root";
  // TFile* fileIn = TFile::Open(fileName, "read");
  // TTree* treeIn = fileIn->Get<TTree>("rTree");
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  TFile *fileOut1 = TFile::Open("out/reco_sim_comp_tracks.root", "recreate");

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
  auto *sim_vtx_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2SimParticles.", &sim_vtx_matching);

  // AnalysisTree::Configuration* config =
  // fileIn->Get<AnalysisTree::Configuration>("Configuration");
  // treeIn->GetConfiguration() instead of congif
  const int sp = treeIn->GetConfiguration()
                     ->GetBranchConfig("SimParticles")
                     .GetFieldId("p"); // simulated
  const int spT = treeIn->GetConfiguration()
                      ->GetBranchConfig("SimParticles")
                      .GetFieldId("pT");
  const int seta = treeIn->GetConfiguration()
                       ->GetBranchConfig("SimParticles")
                       .GetFieldId("eta");
  const int sphi = treeIn->GetConfiguration()
                       ->GetBranchConfig("SimParticles")
                       .GetFieldId("phi");
  const int srap = treeIn->GetConfiguration()
                       ->GetBranchConfig("SimParticles")
                       .GetFieldId("rapidity");
  const int sx = treeIn->GetConfiguration()
                     ->GetBranchConfig("SimParticles")
                     .GetFieldId("x");
  const int sy = treeIn->GetConfiguration()
                     ->GetBranchConfig("SimParticles")
                     .GetFieldId("y");
  const int sz = treeIn->GetConfiguration()
                     ->GetBranchConfig("SimParticles")
                     .GetFieldId("z");
  const int smother_id = treeIn->GetConfiguration()
                             ->GetBranchConfig("SimParticles")
                             .GetFieldId("mother_id");
  const int rp = treeIn->GetConfiguration()
                     ->GetBranchConfig("VtxTracks")
                     .GetFieldId("p"); // reconstructed
  const int rpT =
      treeIn->GetConfiguration()->GetBranchConfig("VtxTracks").GetFieldId("pT");
  const int reta = treeIn->GetConfiguration()
                       ->GetBranchConfig("VtxTracks")
                       .GetFieldId("eta");
  const int rphi = treeIn->GetConfiguration()
                       ->GetBranchConfig("VtxTracks")
                       .GetFieldId("phi");
  const int mass2 =
      treeIn->GetConfiguration()->GetBranchConfig("TofHits").GetFieldId(
          "mass2"); // tofhits
  const int qp_tof =
      treeIn->GetConfiguration()->GetBranchConfig("TofHits").GetFieldId(
          "qp_tof");

  int N = 0;

  const int Nevents = treeIn->GetEntries();

  // zadania
  TH1F hspx("hspx", "simulated px;px [GeV/c];dN/dpx", 1000, -.4,
            .4); // simulatedulated
  TH1F hspy("hspy", "simulated py;py [GeV/c];dN/dpy", 1000, -.4, .4);
  TH1F hspz("hspz", "simulated pz;pz [GeV/c];dN/dpz", 1000, -.2, .4);
  TH1F hsphi("hsphi", "simulated #phi; #phi [#circ];dN/d#phi; ", 500, -3.5,
             3.5);
  TH1F hsp("hsp", "simulated p;p [GeV/c];dN/dp", 1000, -.1, .6);
  TH1F hspt("hspt", "simulated pT;pT [GeV/c];dN/dpT", 1000, -.1, .6);
  TH1F hseta("hseta", "simulated #eta; #eta; dN/d#eta", 1000, -4, 10);
  TH1F hsrap("hsrap", "simulated rapidity; rapidity; counts", 1000, -3, 6);
  TH1F hsx("hsx", "simulated x; x; dN/dx", 100, -3, 3);
  TH1F hsy("hsy", "simulated y; y; dN/dy", 100, -3, 3);
  TH1F hsz("hsz", "simulated z; z; dN/dz", 2100, -100, 1500);
  TH1F hrpx("hrpx", "reconstructed px;px [GeV/c];dN/dpx", 1000, -2, 2);
  TH1F hrpy("hrpy", "reconstructed py;py [GeV/c];dN/dpy", 1000, -2, 2);
  TH1F hrpz("hrpz", "reconstructed pz;pz [GeV/c];dN/dpz", 1000, -.1, 10);
  TH1F hrphi("hrphi", "reconstructed #phi; #phi [#circ];dN/d#phi", 250, -3.5,
             3.5);
  TH1F hrp("hrp", "reconstructed p;p [GeV/c];dN/dp", 1000, -.1, 10);
  TH1F hrpt("hrpt", "reconstructed pT;pT [GeV/c];dN/dpT", 1000, -.1, 2.5);
  TH1F hreta("hreta", "reconstructed #eta; #eta; dN/d#eta", 1000, 0, 6);
  TH2F hcpx("hcpx", "correlation px; reconstructed px; simulated px", 1000, -3,
            3, 1000, -3, 3); // correlationrelation
  TH2F hcpy("hcpy", "correlation py; reconstructed py; simulated py", 1000, -3,
            3, 1000, -3, 3);
  TH2F hcpz("hcpz", "correlation pz; reconstructed pz; simulated pz", 1000, -1,
            10, 1000, -1, 10);
  TH2F hcphi("hcphi", "correlation #phi; reconstructed #phi; simulated #phi",
             1000, -3.5, 3.5, 1000, -3.5, 3.5);
  TH2F hcp("hcp", "correlation p; reconstructed p; simulated p", 1000, 0, 6,
           1000, 0, 6);
  TH2F hcpt("hcpt", "correlation pT; reconstructed pT; simulated pT", 1000, -1,
            2.5, 1000, -1, 2.5);
  TH2F hceta("hceta", "correlation #eta; reconstructed #eta; simulated #eta",
             1000, 0, 6, 1000, 0, 6);
  TH2F hcrpxpy("hcpxpy", "correlation reconstructed px py; px; py", 1000, -3, 3,
               1000, -3, 3);

  for (int i = 0; i < Nevents; i++) {
    treeIn->GetEntry(i);

    // simulated
    for (const auto &sim_track : *sim_tracks) {
      const int sim_mother_id = sim_track.GetField<int>(smother_id);
      if (sim_mother_id == -1)
        continue; // excluding bullet particles
      const float sim_px = sim_track.GetPx();
      const float sim_py = sim_track.GetPy();
      const float sim_pz = sim_track.GetPz();
      const float sim_p = sim_track.GetField<float>(sp);
      const float sim_pt = sim_track.GetField<float>(spT);
      const float sim_eta = sim_track.GetField<float>(seta);
      const float sim_phi = sim_track.GetField<float>(sphi);
      const float sim_rap = sim_track.GetField<float>(srap);
      const float sim_x = sim_track.GetField<float>(sx);
      const float sim_y = sim_track.GetField<float>(sy);
      const float sim_z = sim_track.GetField<float>(sz);
      hspx.Fill(sim_px);
      hspy.Fill(sim_py);
      hspz.Fill(sim_pz);
      hsp.Fill(sim_p);
      hspt.Fill(sim_pt);
      hseta.Fill(sim_eta);
      hsphi.Fill(sim_phi);
      hsrap.Fill(sim_rap);
      hsx.Fill(sim_x);
      hsy.Fill(sim_y);
      hsz.Fill(sim_z);
    }
    // reconstructed and correlations
    for (const auto &vtx_track : *vtx_tracks) {
      const float rec_px = vtx_track.GetPx(); // reconstructed
      const float rec_py = vtx_track.GetPy();
      const float rec_pz = vtx_track.GetPz();
      const float rec_p = vtx_track.GetField<float>(rp);
      const float rec_pt = vtx_track.GetField<float>(rpT);
      const float rec_eta = vtx_track.GetField<float>(reta);
      const float rec_phi = vtx_track.GetField<float>(rphi);
      hrpx.Fill(rec_px);
      hrpy.Fill(rec_py);
      hrpz.Fill(rec_pz);
      hrp.Fill(rec_p);
      hrpt.Fill(rec_pt);
      hreta.Fill(rec_eta);
      hrphi.Fill(rec_phi);
      const int sim_id = sim_vtx_matching->GetMatch(
          vtx_track.GetId()); // simulated, match with vtx
      if (sim_id < 0)
        continue;
      const float sim_px = sim_tracks->GetChannel(sim_id).GetPx();
      const float sim_py = sim_tracks->GetChannel(sim_id).GetPy();
      const float sim_pz = sim_tracks->GetChannel(sim_id).GetPz();
      const float sim_p = sim_tracks->GetChannel(sim_id).GetField<float>(sp);
      const float sim_pt = sim_tracks->GetChannel(sim_id).GetField<float>(spT);
      const float sim_eta =
          sim_tracks->GetChannel(sim_id).GetField<float>(seta);
      const float sim_phi =
          sim_tracks->GetChannel(sim_id).GetField<float>(sphi);
      const float sim_rap =
          sim_tracks->GetChannel(sim_id).GetField<float>(srap);
      const float sim_x = sim_tracks->GetChannel(sim_id).GetField<float>(sx);
      const float sim_y = sim_tracks->GetChannel(sim_id).GetField<float>(sy);
      const float sim_z = sim_tracks->GetChannel(sim_id).GetField<float>(sz);
      // correlations
      hcpx.Fill(rec_px, sim_px);
      hcpy.Fill(rec_py, sim_py);
      hcpz.Fill(rec_pz, sim_pz);
      hcp.Fill(rec_p, sim_p);
      hcpt.Fill(rec_pt, sim_pt);
      hceta.Fill(rec_eta, sim_eta);
      hcphi.Fill(rec_phi, sim_phi);
      hcrpxpy.Fill(rec_px, rec_py);
    }
  }

  hspx.Write();
  hspy.Write();
  hspz.Write();
  hsp.Write();
  hspt.Write();
  hseta.Write();
  hsphi.Write();
  hsrap.Write();
  hsx.Write();
  hsy.Write();
  hsz.Write();
  hrpx.Write();
  hrpy.Write();
  hrpz.Write();
  hrp.Write();
  hrpt.Write();
  hreta.Write();
  hrphi.Write();
  hcpx.GetXaxis()->SetTitle("reconstructed");
  hcpx.GetYaxis()->SetTitle("simulated");
  hcpx.Write();
  hcpy.GetXaxis()->SetTitle("reconstructed");
  hcpy.GetYaxis()->SetTitle("simulated");
  hcpy.Write();
  hcpz.GetXaxis()->SetTitle("reconstructed");
  hcpz.GetYaxis()->SetTitle("simulated");
  hcpz.Write();
  hcp.GetXaxis()->SetTitle("reconstructed");
  hcp.GetYaxis()->SetTitle("simulated");
  hcp.Write();
  hcp.GetXaxis()->SetTitle("reconstructed");
  hcp.GetYaxis()->SetTitle("simulated");
  hcpt.Write();
  hceta.GetXaxis()->SetTitle("reconstructed");
  hceta.GetYaxis()->SetTitle("simulated");
  hceta.Write();
  hcphi.GetXaxis()->SetTitle("reconstructed");
  hcphi.GetYaxis()->SetTitle("simulated");
  hcphi.Write();
  hcrpxpy.GetXaxis()->SetTitle("reconstructed");
  hcrpxpy.GetYaxis()->SetTitle("simulated");
  hcrpxpy.Write();

  fileOut1->Close();
}