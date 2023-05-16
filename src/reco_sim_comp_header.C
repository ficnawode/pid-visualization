void reco_sim_comp_header() {

  // TString fileName = "1.analysistree.root";
  // TFile* fileIn = TFile::Open(fileName, "read");
  // TTree* treeIn = fileIn->Get<TTree>("rTree");
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));

  auto *eve_header = new AnalysisTree::EventHeader();
  treeIn->SetBranchAddress("SimEventHeader.", &eve_header);
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

  AnalysisTree::Configuration *config =
      treeIn
          ->GetConfiguration(); // fileIn->Get<AnalysisTree::Configuration>("Configuration");
  const int sx = config->GetBranchConfig("SimEventHeader")
                     .GetFieldId("vtx_x"); // simulated
  const int sy = config->GetBranchConfig("SimEventHeader").GetFieldId("vtx_y");
  const int sz = config->GetBranchConfig("SimEventHeader").GetFieldId("vtx_z");
  const int spsi =
      config->GetBranchConfig("SimEventHeader").GetFieldId("psi_RP");
  const int rx = config->GetBranchConfig("RecEventHeader")
                     .GetFieldId("vtx_x"); // reconstructed
  const int ry = config->GetBranchConfig("RecEventHeader").GetFieldId("vtx_y");
  const int rz = config->GetBranchConfig("RecEventHeader").GetFieldId("vtx_z");
  const int rchi2 =
      config->GetBranchConfig("RecEventHeader").GetFieldId("vtx_chi2");
  const int rEpsd =
      config->GetBranchConfig("RecEventHeader").GetFieldId("Epsd");
  const int sb = config->GetBranchConfig("SimEventHeader").GetFieldId("b");
  const int rm = config->GetBranchConfig("RecEventHeader").GetFieldId("M");

  int N = 0;

  const int Nevents = treeIn->GetEntries();

  // zadania
  TFile *fileOut1 = TFile::Open("out/reco_sim_comp_header.root", "recreate");
  TH1F hsx("hsx", "simulated x; x [cm]; dN/dx", 100, -.6,
           .6); // simulatedulated
  TH1F hsy("hsy", "simulated y; y [cm]; dN/dy", 100, -.6, .6);
  TH1F hsz("hsz", "simulated z; z [cm]; dN/dz", 1000, -.08, .08);
  TH1F hspsi("hspsi", "simulated #psi; #psi [rad]; dN/dx", 200, -.1, 6.5);
  TH1F hrx("hrx", "reconstructed x; x [cm]; dN/dx", 100, -.6,
           .6); // reconstructedonstructed
  TH1F hry("hry", "reconstructed y; x [cm]; dN/dy", 100, -.6, .6);
  TH1F hrz("hrz", "reconstructed z; x [cm]; dN/dz", 100, -.08, .08);
  TH1F hrchi2("hrchi2", "reconstructed #chi^{2}/NDF; #chi^{2}/NDF; counts", 100,
              0, 3);
  TH1F hEpsd("hEpsd", "reconstructed Epsd; Epsd [GeV]; dN/dEpsd", 100, -1, 60);
  TH1F hm("hm", "reconstructed multiplicity; multiplicty [STS]; counts", 100, 0,
          600);
  TH1F hb("hb", "reconstructed impact parameter; b [fm]; dN/db", 100, 0, 16);
  TH2F hcx("hcx", "correlation px; reconstructed x; simulated x", 1000, -.6, .6,
           1000, -.6, .6); // correlation
  TH2F hcy("hcy", "correlation py; reconstructed y; simulated y", 1000, -.6, .6,
           1000, -.6, .6);
  TH2F hcz("hcz", "correlation pz; reconstructed z; simulated z", 1000, -.08,
           .08, 1000, -.08, .08);
  TH2F hcbm("hcbm", "correlation (b, multiplicty); b [fm]; multiplicty [STS]",
            100, 0, 16, 600, 0, 600);
  TH1F hdx("hdx", "difference x_{rec} - x_{sim};#Delta x [cm]; dN/d#Delta x",
           100, -.01, .01); // simulated
  TH1F hdy("hdy", "diiference y_{rec} - y_{sim};#Delta y [cm]; dN/d#Delta y",
           100, -.01, .01);
  TH1F hdz("hdz", "difference z_{rec} -  z_{sim};#Delta z [cm]; dN/d#Delta z",
           1000, -.2, .2);

  for (int i = 0; i < Nevents; i++) {
    treeIn->GetEntry(i);
    const float sim_x = eve_header->GetField<float>(sx); // simulated
    hsx.Fill(sim_x);
    const float sim_y = eve_header->GetField<float>(sy);
    hsy.Fill(sim_y);
    const float sim_z = eve_header->GetField<float>(sz);
    hsz.Fill(sim_z);
    const float sim_psi = eve_header->GetField<float>(spsi);
    hspsi.Fill(sim_psi);
    const float rec_x = rec_header->GetField<float>(rx); // reconstructed
    const float rec_y = rec_header->GetField<float>(ry);
    const float rec_z = rec_header->GetField<float>(rz);
    hrx.Fill(rec_x);
    hry.Fill(rec_y);
    hrz.Fill(rec_z);
    hrchi2.Fill(rec_header->GetField<float>(rchi2));
    const float rec_epsd = rec_header->GetField<float>(rEpsd);
    hEpsd.Fill(rec_epsd);
    const int rec_m = rec_header->GetField<int>(rm);
    hm.Fill(rec_m);
    const float sim_b = eve_header->GetField<float>(sb);
    hb.Fill(sim_b);
    hcbm.Fill(sim_b, rec_m);
    hcx.Fill(sim_x, rec_x); // correlations
    hcy.Fill(sim_y, rec_y);
    hcz.Fill(sim_z, rec_z);
    hdx.Fill(sim_x - rec_x); // difference
    hdy.Fill(sim_y - rec_y);
    hdz.Fill(sim_z - rec_z);
  }
  gStyle->SetOptStat(0);
  hsx.Write();
  hsy.Write();
  hsz.Write();
  hspsi.Write();
  hrx.Write();
  hry.Write();
  hrz.Write();
  hrchi2.Write();
  hEpsd.Write();
  hm.Write();
  hb.Write();
  hcx.Write();
  hcy.Write();
  hcz.Write();
  hcbm.Write();
  hdx.Write();
  hdy.Write();
  hdz.Write();

  fileOut1->Close();
}