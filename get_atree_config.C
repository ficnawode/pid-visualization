
void get_atree_config(TString filename = "in/1001.analysistree.root") {
  auto file = new TFile(filename);
  cout << file << endl;
  auto config = (AnalysisTree::Configuration *)file->Get("Configuration");
  config->Print();
}