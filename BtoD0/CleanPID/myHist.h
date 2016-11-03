#ifndef MYHIST_H_INCLUDED
#define MYHIST_H_INCLUDED
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
//hist
//signal test
TH1D* hD0;
TH1D* hD0Wg;
TH1F* hpt;
TH1F* hptWg;
TH1F* hppt;
TH1F* hpptWg;

//significance
TH2F* h2massSig;
TH2F* h2massBg;
TH2F* h2massSigDoubleCount;

//efficiency
TH1F* hD0pT_noCut;
TH1F* hD0pT_Cut;

//D0 Dca
TH2F* h2D0DcaSig;
TH2F* h2D0DcaBg;
TH2F* h2D0DcaSigDoubleCount;

//tree
TNtuple* nt_sig;
TNtuple* nt_bg;

//cuts
namespace anaCuts
{
    int   const nPtBins = 5;
    float const PtEdge[nPtBins+1] = {0., 1., 2., 3., 5., 15.};
    
    float const pt = 0.6;
    float const eta = 1.0;
    
    float const rapidity = 1.0;
    
    float const massMin = 1.828;
    float const massMax = 1.892;
    
    float const rightHalfLowEdge = -1.54696;
    float const rightHalfHighEdge = 1.59464;
    
    //float const dcaV0ToPv[nPtBins] = {60, 60, 60, 60, 60};
    //float const decayLength[nPtBins] = {80, 80, 80, 80, 80};
    //float const cosTheta[nPtBins] = {0.95, 0.95, 0.95, 0.95, 0.95};
    //float const dcaDaughters[nPtBins] = {100, 100, 100, 100, 100}; //0.0050;
    //float const kDca[nPtBins] = {60, 60, 60, 60, 60};//0.008, // minimum
    //float const pDca[nPtBins] = {60, 60, 60, 60, 60};//0.008
    float const dcaV0ToPv[nPtBins] = {61, 49, 38, 38, 40};
    float const decayLength[nPtBins] = {145, 181, 212, 247, 259};
    float const cosTheta[nPtBins] = {0.95, 0.95, 0.95, 0.95, 0.95};
    float const dcaDaughters[nPtBins] = {84, 66, 57, 50, 60}; //0.0050;
    float const kDca[nPtBins] = {103, 91, 95, 79, 58};//0.008, // minimum
    float const pDca[nPtBins] = {110, 111, 86, 81, 62};//0.008
    
    // all variables with a phys prefix are for final physics plots
    //int const physNCentralities = 4;
    //int const physCentralityEdges[physNCentralities+1] = {0,3,5,6,9}; // 40-80, 20-40, 10-20, 0-10
    //TString const physCentralityName[physNCentralities] = {"40-80%","20-40%","10-20%","0-10%"};
    int const physNCentralities = 9;
    int const physCentralityEdges[physNCentralities+1] = {0,0,1,2,3,4,5,6,7,8}; // 70-80, 60-70, 50-60, 40-50, 30-40, 20-30, 10-20, 5-10, 0-5
    TString const physCentralityName[physNCentralities] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
    int   const physNPtBins = 10;
    double const physPtEdge[physNPtBins+1] = {0., 0.5, 1.0, 1.50, 2.0, 2.50, 3.0, 3.50, 4.0, 5.0, 8.};
    //int   const physNPtBins = 7;
    //double const physPtEdge[physNPtBins+1] = {0., 0.7, 1.1, 1.6, 2.2, 3.0, 5.0, 8.0};
}

//define hist
void DefineHist() {
    TH1::SetDefaultSumw2();
    hD0 = new TH1D("hD0","hD0",2,0,2);
    hD0->SetDirectory(0);
    hD0Wg = new TH1D("hD0Wg","hD0Wg",2,0,2);
    hD0Wg->SetDirectory(0);
    hpt = new TH1F("hpt","hpt",200,0,20);
    hpt->SetDirectory(0);
    hptWg = new TH1F("hptWg","hptWg",200,0,20);
    hptWg->SetDirectory(0);
    hppt = new TH1F("hppt","hppt",200,0,20);
    hppt->SetDirectory(0);
    hpptWg = new TH1F("hpptWg","hpptWg",200,0,20);
    hpptWg->SetDirectory(0);
    h2massSig = new TH2F("h2massSig","h2massSig;D0 pT;mass(kpi)",200,0,20,50,1.6,2.1);
    h2massSig->SetDirectory(0);
    h2massBg = new TH2F("h2massBg","h2massBg;D0 pT;mass(kpi)",200,0,20,50,1.6,2.1);
    h2massBg->SetDirectory(0);
    h2massSigDoubleCount = new TH2F("h2massSigDoubleCount","h2massSigDoubleCount;D0 pT;mass(kpi)",200,0,20,50,1.6,2.1);
    h2massSigDoubleCount->SetDirectory(0);
    hD0pT_noCut = new TH1F("hD0pT_noCut","hD0pT_noCut",200,0,20);
    hD0pT_noCut->SetDirectory(0);
    hD0pT_Cut = new TH1F("hD0pT_Cut","hD0pT_Cut",200,0,20);
    hD0pT_Cut->SetDirectory(0);
    h2D0DcaSig = new TH2F("h2D0DcaSig","h2D0DcaSig;D0 pT;D0 Dca",200,0,20,200,0,0.1);
    h2D0DcaSig->SetDirectory(0);
    h2D0DcaBg = new TH2F("h2D0DcaBg","h2D0DcaBg;D0 pT;D0 Dca",200,0,20,200,0,0.1);
    h2D0DcaBg->SetDirectory(0);
    h2D0DcaSigDoubleCount = new TH2F("h2D0DcaSigDoubleCount","h2D0DcaSigDoubleCount;D0 pT;D0 Dca",200,0,20,200,0,0.1);
    h2D0DcaSigDoubleCount->SetDirectory(0);
}
//write hist
void WriteHist(TFile* fout) {
    fout->cd();
    hD0->Write();
    hD0Wg->Write();
    hpt->Write();
    hptWg->Write();
    hppt->Write();
    hpptWg->Write();
    h2massSig->Write();
    h2massBg->Write();
    h2massSigDoubleCount->Write();
    hD0pT_noCut->Write();
    hD0pT_Cut->Write();
    h2D0DcaSig->Write();
    h2D0DcaBg->Write();
    h2D0DcaSigDoubleCount->Write();
}
void DeleteHist() {
    if(hD0) { delete hD0; hD0 = NULL; }
    if(hD0Wg) { delete hD0Wg; hD0Wg = NULL; }
    if(hpt) { delete hpt; hpt = NULL; }
    if(hptWg) { delete hptWg; hptWg = NULL; }
    if(hppt) { delete hppt; hppt = NULL; }
    if(hpptWg) { delete hpptWg; hpptWg = NULL; }
    if(h2massSig) { delete h2massSig; h2massSig = NULL; }
    if(h2massBg) { delete h2massBg; h2massBg = NULL; }
    if(h2massSigDoubleCount) { delete h2massSigDoubleCount; h2massSigDoubleCount = NULL; }
    if(hD0pT_noCut) { delete hD0pT_noCut; hD0pT_noCut = NULL; }
    if(hD0pT_Cut) { delete hD0pT_Cut; hD0pT_Cut = NULL; }
    if(h2D0DcaSig) { delete h2D0DcaSig; h2D0DcaSig = NULL; }
    if(h2D0DcaBg) { delete h2D0DcaBg; h2D0DcaBg = NULL; }
    if(h2D0DcaSigDoubleCount) { delete h2D0DcaSigDoubleCount; h2D0DcaSigDoubleCount = NULL; }
}
void DefineTree() {
    int BufSize = (int)pow(2., 16.);
    nt_sig = new TNtuple("Signal", "", "pid:m:pt:eta:y:phi:pTWg:doubleCountWg:" // MC D0
                     "rM:rPt:rEta:rY:rPhi:" // Rc D0
                     "dca12:decayLength:dcaD0ToPv:dcaXY:dcaZ:cosTheta:" // Rc pair
                     "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                     "kRM:kRPt:kREta:kRY:kRPhi:kRDca:kRSDca:kRDcaXY:kRDcaZ:kTpc:kTof:" // Rc Kaon
                     "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion
                     "pRM:pRPt:pREta:pRY:pRPhi:pRDca:pRSDca:pRDcaXY:pRDcaZ:pTpc:pTof:" // Rc Pion
                     "kHft:pHft:primaryPt:primaryY:primaryPhi", BufSize); //add primary pt, y, phi
    nt_sig->SetAutoSave(-500000); // autosave every 0.5 Mbytes
    
    nt_bg = new TNtuple("BackGround", "", "kPt:kpTWg:kEta:kphi:kTpc:kTof:"
                        "piPt:pipTWg:piEta:piphi:piTpc:piTof:"
                        "rD0M:rD0Pt:rD0Y:cosTheta:decayLength:dcaToPv:kDca:piDca:dcaDaughters:kHft:piHft", BufSize); //RC LC
    nt_bg->SetAutoSave(-500000); // autosave every 0.5 Mbytes
}

#endif // VARIABLE_H_INCLUDED
