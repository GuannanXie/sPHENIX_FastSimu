void runBtoD(int npart = 20000, TString output = "D0.toyMc.root", TString particleName = "D0", Int_t mWriteType = 1, bool isCombinB = false)
{
  gSystem->Load("toyMcBtoD_C");
  toyMcBtoD(npart,output,particleName,mWriteType,isCombinB);
}
