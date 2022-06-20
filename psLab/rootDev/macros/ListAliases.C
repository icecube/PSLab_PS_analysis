
void ListAliases(TTree* tree) {
  int n = tree->GetListOfAliases()->GetEntries();
  for (int i=0; i<n; ++i) {
    cout << i << ": \"";
    cout << tree->GetListOfAliases()->At(i)->GetName() << "\" , \"";
    cout << tree->GetListOfAliases()->At(i)->GetTitle() << "\"\n";
  }
}
