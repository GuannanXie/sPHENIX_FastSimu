/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0/B+/B0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Xiaolong Chen (xlchen@lbl.gov)
 *            Guannan Xie (guannanxie@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TStopwatch.h"
// #include "TSystem.h"
// #include "TMemStat.h"
#include "myHist.h"

using namespace std;

void setDecayChannels(int const mMode);
bool GetD0FromBDecay(int& kf, TLorentzVector& b, TVector3& v0_D0, TClonesArray& daughters);
void decayAndFill(int const kf, TLorentzVector* b, TVector3 const v0_D0, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void fillDouble(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00, TVector3 kPos, TVector3 pPos);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
//bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
//bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
//bool reconstructD0(int const centrality, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
//int getPhiIndexDca(double);
int getD0PtIndex(float const pt);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
//int getPhiIndexHftRatio(double);

TPythia6Decayer* pydecay;
TFile* result;
float primaryPt;
float primaryY;
float primaryPhi;

const Int_t nParticles = 3;  //0--pi, 1--k, 2--p
const Int_t nCentHftRatio = 1;
const Float_t PI = TMath::Pi();

// HFT ratio binning
const Int_t nEtasHftRatio = 1;
const Int_t nVzsHftRatio = 1;
const Int_t nPtBinsHftRatio = 17;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] = { -1.0, 1.0 };
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = { -6.0e4,  6.0e4};
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] = { 0.3 ,  0.4 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 2.6 , 3.0 , 4. , 20. };

// DCA binning
int const nVzsDca = 1;
float const VzEdgeDca[nVzsDca + 1] = { -6.e4, 6.e4};

int const nEtasDca = 1;
float const EtaEdgeDca[nEtasDca + 1] = { -1.0, 1.0};

const Int_t nPtBinsDca = 18;
const Double_t ptEdgeDca[nPtBinsDca + 1] =
{ 0.3 ,  0.4 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 2.6 , 3.0 , 4. , 6., 20. };

//input 1--hft match ratio
const float mHftMatchRatio = 1.0; //merged to tpc eff

//input 2--dcaXY~dcaZ 2D smearing
TH2F* h2Dca[nParticles][nPtBinsDca];    //0--pi, 1--k, 2--p

//input 3--momentum resolution
TF1* fPionMomResolution = NULL;
TF1* fKaonMomResolution = NULL; //currentlly use the same momentum resolution as PION
TF1* fProtonMomResolution = NULL; //not add proton first

//input 4--back ground particle number from hijing
TH3I* h3PiKPNumber = NULL;
//pT spectra from paper
TH1F* hpiPtWg = NULL;
TH1F* hkPtWg = NULL;
TH1F* hpPtWg = NULL;

//input 5--tpc tracking efficieny
TF1* fTpcK = NULL;
TF1* fTpcPi = NULL;
TF1* fTpcP = NULL;

//input 6--tof matching efficieny
TF1* fTofK = NULL;
TF1* fTofPi = NULL;
TF1* fTofP = NULL;

//input 7--D0/B pT weight
TH1F* hD0PtWgFonll = NULL;
TH1F* hBPtWgFonll = NULL;

TString outFileName;// = "D0.toyMc.root";
std::pair<float, float> const momentumRange(0, 20);

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.2; //for D0 from B decay may be larger
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCentHftRatio] = {4.};

int mMode = 0;  //1--D0,  2--B0,  3--Bpm
bool mUseHist = true;
bool mUseTree = false;
//const float M_D_0 = 1.86484;
const float M_B_0 = 5.27958;
const float M_B_PLUS = 5.27926;
//============== main  program ==================
void toyMcBtoD(int npart = 10, TString output = "D0.toyMc.root", TString particleName = "D0", Int_t mWriteType=1, bool isCombinB = false)
{
    // TMemStat mem;
    // mem.Enable();
    TStopwatch*   stopWatch = new TStopwatch();
    stopWatch->Start();
    gRandom->SetSeed();
    
    //In this part, if you want, add the random number to choose specific D0, B0, and B+ ratio
    //I plan to combine B0 and B+
    //
    
    if(particleName.CompareTo("D0")==0)      mMode = 1;
    else if(particleName.CompareTo("B0")==0) mMode = 2;
    else if(particleName.CompareTo("Bpm")==0) mMode = 3;
    else {
        cout << ">>>>>>>Please input the right particle name: D0, B0, Bpm ......<<<<<<<<" << endl;
        exit(1);
    }
    
    //loat input files
    outFileName = output;
    if(!output.Contains(".root")) outFileName += ".root";
    bookObjects();
    
    if(mWriteType==2) {mUseHist=false; mUseTree=true;}
    else if(mWriteType==3) {mUseHist=true; mUseTree=true;}
    else {mUseHist=true; mUseTree=false;}  //default
    if(mUseHist) DefineHist();
    if(mUseTree) DefineTree();
    
    bool mSignal = true;
    bool mBackGround = false;
    if(mSignal) {
        pydecay = TPythia6Decayer::Instance();
        pydecay->Init();
        TPythia6::Instance()->SetMRPY(1, 88158204); //random seed number
    }
    for (int ipart = 0; ipart < npart; ipart++)
    {
        if (npart>1000 && !(ipart % (npart/100)))
            cout << "____________ ipart = " << ipart / static_cast<float>(npart) << " ________________" << endl;
        
        //background reconstruction--primary k/pi from hijing
        if(mBackGround) {
            vector<float> vecPt_pi;
            vector<float> vecEta_pi;
            vector<float> vecPhi_pi;
            vector<TVector3> vecPos_pi;
            vector<int> vecType_pi;
            vector<float> vecPt_k;
            vector<float> vecEta_k;
            vector<float> vecPhi_k;
            vector<TVector3> vecPos_k;
            vector<int> vecType_k;
            vector<float> vecPt_piMisPid;
            vector<float> vecEta_piMisPid;
            vector<float> vecPhi_piMisPid;
            vector<TVector3> vecPos_piMisPid;
            vector<float> vecPt_kMisPid;
            vector<float> vecEta_kMisPid;
            vector<float> vecPhi_kMisPid;
            vector<TVector3> vecPos_kMisPid;
            vector<float> vecPt_pMisPid;
            vector<float> vecEta_pMisPid;
            vector<float> vecPhi_pMisPid;
            vector<TVector3> vecPos_pMisPid;
            Double_t nPion_tmp, nKaon_tmp, nProton_tmp;
            h3PiKPNumber->GetRandom3(nPion_tmp,nKaon_tmp,nProton_tmp);
            int nPion = floor(nPion_tmp);
            int nKaon = floor(nKaon_tmp);
            int nProton = floor(nProton_tmp);
            // k smear
            TVector3 kPos(0,0,0);
            TLorentzVector rckFourMom;
            TVector3 rckPos;
            TLorentzVector kFourMom;
            for(int ik=0; ik<nKaon; ik++) {
                float pt = gRandom->Uniform(0.6,20);
                float eta = gRandom->Uniform(-1,1);
                float phi = gRandom->Uniform(-PI,PI);
                kFourMom.SetPtEtaPhiM(pt, eta , phi, M_KAON_MINUS);
                rckFourMom = smearMom(kFourMom,fKaonMomResolution);
                rckPos = smearPosData(1, 0, 7, rckFourMom, kPos); //1--kaon 0--pi
                vecPt_k.push_back(rckFourMom.Perp());
                vecEta_k.push_back(rckFourMom.Eta());
                vecPhi_k.push_back(rckFourMom.Phi());
                vecPos_k.push_back(rckPos);
                vecType_k.push_back(1);  //0--pi, 1--k, 2--p
                if(pt>1.6)  // mis pid
                {
                    vecPt_kMisPid.push_back(rckFourMom.Perp());
                    vecEta_kMisPid.push_back(rckFourMom.Eta());
                    vecPhi_kMisPid.push_back(rckFourMom.Phi());
                    vecPos_kMisPid.push_back(rckPos);
                }
            }
            // pi smear
            TVector3 piPos(0,0,0);
            TLorentzVector rcpiFourMom;
            TVector3 rcpiPos;
            TLorentzVector piFourMom;
            for(int ipi=0; ipi<nPion; ipi++) {
                float pt = gRandom->Uniform(0.6,20);
                float eta = gRandom->Uniform(-1,1);
                float phi = gRandom->Uniform(-PI,PI);
                piFourMom.SetPtEtaPhiM(pt,eta,phi,M_PION_PLUS);
                rcpiFourMom = smearMom(piFourMom,fPionMomResolution);
                rcpiPos = smearPosData(0, 0, 7, rcpiFourMom, piPos); //1--k, 0--pi
                vecPt_pi.push_back(rcpiFourMom.Perp());
                vecEta_pi.push_back(rcpiFourMom.Eta());
                vecPhi_pi.push_back(rcpiFourMom.Phi());
                vecPos_pi.push_back(rcpiPos);
                vecType_pi.push_back(0);   //0--pi, 1--k, 2--p
                if(pt>1.6)  // mis pid
                {
                    vecPt_piMisPid.push_back(rcpiFourMom.Perp());
                    vecEta_piMisPid.push_back(rcpiFourMom.Eta());
                    vecPhi_piMisPid.push_back(rcpiFourMom.Phi());
                    vecPos_piMisPid.push_back(rcpiPos);
                }
            }
            // p smear
            TVector3 pPos(0,0,0);
            TLorentzVector rcpFourMom;
            TVector3 rcpPos;
            TLorentzVector pFourMom;
            for(int ip=0; ip<nProton; ip++) {
                float pt = gRandom->Uniform(0.6,20);
                float eta = gRandom->Uniform(-1,1);
                float phi = gRandom->Uniform(-PI,PI);
                pFourMom.SetPtEtaPhiM(pt,eta,phi,M_PROTON);
                rcpFourMom = smearMom(pFourMom,fProtonMomResolution);
                rcpPos = smearPosData(2, 0, 7, rcpFourMom, pPos); //1--k, 0--pi, 2--p
                if(pt>1.6)  // mis pid
                {
                    vecPt_pMisPid.push_back(rcpFourMom.Perp());
                    vecEta_pMisPid.push_back(rcpFourMom.Eta());
                    vecPhi_pMisPid.push_back(rcpFourMom.Phi());
                    vecPos_pMisPid.push_back(rcpPos);
                }
            }
            // pi mixed into k
            vecPt_k.insert(vecPt_k.end(),vecPt_piMisPid.begin(),vecPt_piMisPid.end());
            vecEta_k.insert(vecEta_k.end(),vecEta_piMisPid.begin(),vecEta_piMisPid.end());
            vecPhi_k.insert(vecPhi_k.end(),vecPhi_piMisPid.begin(),vecPhi_piMisPid.end());
            vecPos_k.insert(vecPos_k.end(),vecPos_piMisPid.begin(),vecPos_piMisPid.end());
            vecType_k.insert(vecType_k.end(),vecPos_piMisPid.size(),0);
            // p mixed into k
            vecPt_k.insert(vecPt_k.end(),vecPt_pMisPid.begin(),vecPt_pMisPid.end());
            vecEta_k.insert(vecEta_k.end(),vecEta_pMisPid.begin(),vecEta_pMisPid.end());
            vecPhi_k.insert(vecPhi_k.end(),vecPhi_pMisPid.begin(),vecPhi_pMisPid.end());
            vecPos_k.insert(vecPos_k.end(),vecPos_pMisPid.begin(),vecPos_pMisPid.end());
            vecType_k.insert(vecType_k.end(),vecPos_pMisPid.size(),2);
            // k mixed into pi
            vecPt_pi.insert(vecPt_pi.end(),vecPt_kMisPid.begin(),vecPt_kMisPid.end());
            vecEta_pi.insert(vecEta_pi.end(),vecEta_kMisPid.begin(),vecEta_kMisPid.end());
            vecPhi_pi.insert(vecPhi_pi.end(),vecPhi_kMisPid.begin(),vecPhi_kMisPid.end());
            vecPos_pi.insert(vecPos_pi.end(),vecPos_kMisPid.begin(),vecPos_kMisPid.end());
            vecType_pi.insert(vecType_pi.end(),vecPos_kMisPid.size(),1);
            // p mixed into pi
            vecPt_pi.insert(vecPt_pi.end(),vecPt_pMisPid.begin(),vecPt_pMisPid.end());
            vecEta_pi.insert(vecEta_pi.end(),vecEta_pMisPid.begin(),vecEta_pMisPid.end());
            vecPhi_pi.insert(vecPhi_pi.end(),vecPhi_pMisPid.begin(),vecPhi_pMisPid.end());
            vecPos_pi.insert(vecPos_pi.end(),vecPos_pMisPid.begin(),vecPos_pMisPid.end());
            vecType_pi.insert(vecType_pi.end(),vecPos_pMisPid.size(),2);
            // clear
            vecPt_kMisPid.clear(); vecEta_kMisPid.clear(); vecPhi_kMisPid.clear(); vecPos_kMisPid.clear();
            vecPt_piMisPid.clear(); vecEta_piMisPid.clear(); vecPhi_piMisPid.clear(); vecPos_piMisPid.clear();
            vecPt_pMisPid.clear(); vecEta_pMisPid.clear(); vecPhi_pMisPid.clear(); vecPos_pMisPid.clear();
            //cout << "nKaons = " << vecPt_k.size() << endl;
            //cout << "nPions = " << vecPt_pi.size() << endl;
            //reconstruction
            for(int ik=0; ik<vecPt_k.size(); ik++) {
                TVector3 vertex(0,0,0);
                float kpt = vecPt_k[ik];
                float keta = vecEta_k[ik];
                float kphi = vecPhi_k[ik];
                if(kpt<0.6) continue;
                //cout << kpt << "\t" << keta << "\t" << kphi << endl;
                rckFourMom.SetPtEtaPhiM(kpt, keta , kphi, M_KAON_MINUS);
                rckPos = vecPos_k[ik];
                float kRDca = dcaSigned(rckFourMom.Vect(), rckPos, vertex);
                if(fabs(kRDca)<50) continue;
                // get weight--- 0--pi, 1--k, 2--p
                float kPtWg = 1.;
                float ktof = 1.;
                float ktpc = 1.;
                if(vecType_k[ik]==0) {
                    kPtWg = hpiPtWg->GetBinContent(hpiPtWg->FindBin(kpt));
                    ktof = (kpt<=1.6 ? fTofPi->Eval(kpt) : 1.0);
                    ktpc = fTpcPi->Eval(kpt);
                }
                else if(vecType_k[ik]==1) {
                    kPtWg = hkPtWg->GetBinContent(hkPtWg->FindBin(kpt));
                    ktof = (kpt<=1.6 ? fTofK->Eval(kpt) : 1.0);
                    ktpc = fTpcK->Eval(kpt);
                }
                else if(vecType_k[ik]==2) {
                    kPtWg = hpPtWg->GetBinContent(hpPtWg->FindBin(kpt));
                    ktof = (kpt<=1.6 ? fTofP->Eval(kpt) : 1.0);
                    ktpc = fTpcP->Eval(kpt);
                }
                else {
                    continue;
                }
                for(int ipi=0; ipi<vecPt_pi.size(); ipi++) {
                    float pipt = vecPt_pi[ipi];
                    float pieta = vecEta_pi[ipi];
                    float piphi = vecPhi_pi[ipi];
                    if(pipt<0.6) continue;
                    if(fabs(kpt-pipt) +fabs(keta-pieta) + fabs(kphi-piphi) < 2e-6) continue;
                    
                    rcpiFourMom.SetPtEtaPhiM(pipt,pieta,piphi,M_PION_PLUS);
                    rcpiPos = vecPos_pi[ipi];
                    float piRDca = dcaSigned(rcpiFourMom.Vect(), rcpiPos, vertex);
                    if(fabs(piRDca)<50) continue;
                    // get weight--- 0--pi, 1--k, 2--p
                    float piPtWg = 1.;
                    float pitof = 1.;
                    float pitpc = 1.;
                    if(vecType_pi[ipi]==0) {
                        piPtWg = hpiPtWg->GetBinContent(hpiPtWg->FindBin(pipt));
                        pitof = (pipt<=1.6 ? fTofPi->Eval(pipt) : 1.0);
                        pitpc = fTpcPi->Eval(pipt);
                    }
                    else if(vecType_pi[ipi]==1) {
                        piPtWg = hkPtWg->GetBinContent(hkPtWg->FindBin(pipt));
                        pitof = (pipt<=1.6 ? fTofK->Eval(pipt) : 1.0);
                        pitpc = fTpcK->Eval(pipt);
                    }
                    else if(vecType_pi[ipi]==2) {
                        piPtWg = hpPtWg->GetBinContent(hpPtWg->FindBin(pipt));
                        pitof = (pipt<=1.6 ? fTofP->Eval(pipt) : 1.0);
                        pitpc = fTpcP->Eval(pipt);
                    }
                    else {
                        continue;
                    }
                    
                    TLorentzVector const rMom =  rckFourMom + rcpiFourMom;
                    if (rMom.M() > 2.1 || rMom.M() < 1.6) continue;
                    
                    TVector3 v0;
                    float const dca12 = dca1To2(rckFourMom.Vect(), rckPos, rcpiFourMom.Vect(), rcpiPos, v0);
                    if (dca12 > 100.) continue;
                    
                    float const decayLength = (v0 - vertex).Mag();
                    // if (decayLength < 40) return;
                    float const dcaToPv = dca(rMom.Vect(), v0, vertex);
                    // if (dcaToPv > 100.) return;//this cut was removed for b->D decay
                    float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
                    //if (cosTheta < 0.9) continue; //default low limit cut is 0.97
                    
                    // save
                    if(mUseTree) {
                        float arr[25];
                        int iArr = 0;
                        arr[iArr++] = kpt;
                        arr[iArr++] = kPtWg;
                        arr[iArr++] = keta;
                        arr[iArr++] = kphi;
                        arr[iArr++] = ktpc;
                        arr[iArr++] = ktof;
                        
                        arr[iArr++] = pipt;
                        arr[iArr++] = piPtWg;
                        arr[iArr++] = pieta;
                        arr[iArr++] = piphi;
                        arr[iArr++] = pitpc;
                        arr[iArr++] = pitof;
                        
                        arr[iArr++] = rMom.M();
                        arr[iArr++] = rMom.Perp();
                        arr[iArr++] = rMom.Rapidity();
                        arr[iArr++] = cosTheta;
                        arr[iArr++] = decayLength;
                        arr[iArr++] = dcaToPv;
                        arr[iArr++] = kRDca;
                        arr[iArr++] = piRDca;
                        arr[iArr++] = dca12;
                        
                        arr[iArr++] = mHftMatchRatio;
                        arr[iArr++] = mHftMatchRatio;
                        nt_bg->Fill(arr);
                        if (ipart % 1000 == 1) nt_bg->AutoSave("SaveSelf");
                    }
                    if(mUseHist) {
                        float weight = mHftMatchRatio * mHftMatchRatio * ktof * pitof * ktpc * pitpc * kPtWg * piPtWg;
                        int ptIndex = getD0PtIndex(rMom.Perp());
                        if(ptIndex < 0) continue;
                        bool isD0 = ( fabs(piRDca) > anaCuts::pDca[ptIndex]
                                     && fabs(kRDca) > anaCuts::kDca[ptIndex]
                                     && dca12 < anaCuts::dcaDaughters[ptIndex]
                                     && decayLength > anaCuts::decayLength[ptIndex]
                                     && cosTheta > anaCuts::cosTheta[ptIndex]
                                     //&& dcaToPv < anaCuts::dcaV0ToPv[ptIndex]
                                     && kpt > 0.6 && pipt > 0.6
                                     && fabs(rMom.Rapidity()) < 1.0
                                     );
                        bool isD0_prompt = isD0 && cosTheta>0.95 && dcaToPv < anaCuts::dcaV0ToPv[ptIndex];
                        if(isD0) {
                            h2massBg->Fill(rMom.Perp(),rMom.M(),weight);
                            if(rMom.M()>1.845 && rMom.M()<1.885) h2D0DcaBg->Fill(rMom.Perp(),dcaToPv/1.e4,weight); //about 2Sigma
                        }
                        if(isD0_prompt) {
                            h2massBg_D0->Fill(rMom.Perp(),rMom.M(),weight);
							hptTotal->Fill(rMom.Perp(),weight);
							if(kpt>1.6 && pipt>1.6) hpt1->Fill(rMom.Perp(),weight);
							else if(kpt<1.6 && pipt<1.6) hpt3->Fill(rMom.Perp(),weight);
							else hpt2->Fill(rMom.Perp(),weight);
                        }
                    }
                }//end pion
            }//end kaon
            vecPt_k.clear(); vecEta_k.clear(); vecPhi_k.clear(); vecPos_k.clear(); vecType_k.clear();
            vecPt_pi.clear(); vecEta_pi.clear(); vecPhi_pi.clear(); vecPos_pi.clear(); vecType_pi.clear();
        }//end background
        
        //signal reconstruction
        if(mSignal) {
            TClonesArray ptl("TParticle", 100);
            TLorentzVector* b_d = new TLorentzVector;
            for (int iD0 = 0; iD0 < 1; iD0++)
            {
                if(isCombinB) {
                    float FR_B0 = 0.4;
                    float FR_Bplus = 0.4;
                    float FR_Sum = FR_B0 + FR_Bplus;
                    float FR_Rdm = FR_Sum*gRandom->Rndm();
                    if(FR_Rdm<FR_B0) mMode = 2;
                    else if(FR_Rdm<FR_B0+FR_Bplus) mMode = 3;
                    else continue;  //this doesn't work, because I close some decay branch
                }
                
                setDecayChannels(mMode);
                //see particle mass: http://www.johnmarcampbell.com/lambdaMixingDoc/da/dae/phys__constants_8h_source.html
                if(mMode==1) getKinematics(*b_d, M_D_0);
                else if(mMode==2) getKinematics(*b_d, M_B_0);
                else if(mMode==3) getKinematics(*b_d, M_B_PLUS);
                else {
                    cout << ">>>>>>>>>>>No exact input particle<<<<<<<<<<" << endl;
                    exit(1);
                }
                
                int KF;
                TVector3 v0_D0(0.,0.,0.);
                if(mMode==1) {
                    KF = 421;
                    if(gRandom->Rndm()<0.5) KF = -KF;
                    decayAndFill(KF, b_d, v0_D0, 1.0, ptl);
                }
                else {
                    if(mMode==2) KF = 511;
                    if(mMode==3) KF = 521;
                    bool ifGetD0 = GetD0FromBDecay(KF, *b_d, v0_D0, ptl);  //B decay to D0
                    if(!ifGetD0) continue;
                    decayAndFill(KF, b_d, v0_D0, 1.0, ptl);
                }
            }//end nD0
            if(mUseTree && (ipart % 6000 == 1)) nt_sig->AutoSave("SaveSelf");
            delete b_d;
        }//end signal
    }//end npart
    
    write();
    if(mUseHist) DeleteHist();
    if(mUseTree) {
        //if(nt_sig) { delete nt_sig; nt_sig = NULL; }
        //if(nt_bg) { delete nt_bg; nt_bg = NULL; }
    }
    
    // mem.Show();
    stopWatch->Stop();
    stopWatch->Print();
}

void setDecayChannels(int const mMode)
{
    // Set some stable particle
    TPythia6::Instance()->SetMDCY(122,1,0); //Set D+/D- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(128,1,0); //Set Ds+/Ds- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(102,1,0); //Set pi0 to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(15,1,0); //Set tau+/tau- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(107,1,0); //Set rho+/rho- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(281,1,0); //Set a_1+/a_1- to be stable particle, not decay
    
    // Set intermediate particle channels -->D0
    // ---not set this, in order to keep the B->D0 branch ratio the same as Pythia
    
    // set D0 decay channel, constrained to D0->k-pi+, D0bar->k+pi-
    for (int idc = 747; idc < 807 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
    TPythia6::Instance()->SetMDME(763, 1, 1);
    
    // set B0 decay channels
    if(mMode==2) {
        for (int idc = 863; idc < 898 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
        for (int idc = 864; idc < 868 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
        for (int idc = 870; idc < 874 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
        TPythia6::Instance()->SetMDME(876, 1, 1);
        for (int idc = 880; idc < 882 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
        for (int idc = 885; idc < 886 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
    }
    
    // set B+ decay channels
    if(mMode==3) {
        for (int idc = 908; idc < 943 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
        for (int idc = 908; idc < 931 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
    }
}

bool GetD0FromBDecay(int& kf, TLorentzVector& b, TVector3& v0_D0, TClonesArray& daughters)
{
    if(gRandom->Rndm()<0.5) kf = -kf;
    pydecay->Decay(kf, &b);
    pydecay->ImportParticles(&daughters);
    
    TLorentzVector D0Mom;
    int nTrk = daughters.GetEntriesFast();
    //cout << "nTrk = " << nTrk << ", mMode = " << mMode << endl;
    bool ifGetD0 = false;
    for (int iTrk = 0; iTrk < nTrk; ++iTrk)
    {
        TParticle* ptl0_pre = (TParticle*)daughters.At(iTrk);
        //cout << "itrk = " << ptl0_pre->GetPdgCode() << endl;
        switch (ptl0_pre->GetPdgCode())
        {
            case 421:
                ptl0_pre->Momentum(D0Mom);
                v0_D0.SetXYZ(ptl0_pre->Vx() * 1000., ptl0_pre->Vy() * 1000., ptl0_pre->Vz() * 1000.); // converted to μm
                kf = 421;
                ifGetD0 = true;
                break;
            case -421:
                ptl0_pre->Momentum(D0Mom);
                v0_D0.SetXYZ(ptl0_pre->Vx() * 1000., ptl0_pre->Vy() * 1000., ptl0_pre->Vz() * 1000.); // converted to μm
                kf = -421;
                ifGetD0 = true;
                break;
            default:
                break;
        }
    }
    daughters.Clear();
    b = D0Mom;
    return ifGetD0;
}

void decayAndFill(int const kf, TLorentzVector* b, TVector3 const v0_D0, double const weight, TClonesArray& daughters)
{
    pydecay->Decay(kf, b);
    pydecay->ImportParticles(&daughters);
    
    TLorentzVector kMom;
    TLorentzVector pMom;
    TVector3 v00;
    
    int nTrk = daughters.GetEntriesFast();
    for (int iTrk = 0; iTrk < nTrk; ++iTrk)
    {
        TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
        
        switch (abs(ptl0->GetPdgCode()))
        {
            case 321:
                ptl0->Momentum(kMom);
                // v00.SetXYZ(0,0,0);
                v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
                break;
            case 211:
                ptl0->Momentum(pMom);
                break;
            default:
                break;
        }
    }
    daughters.Clear();
    
    TVector3 v0 = v00 + v0_D0;
    fill(kf, b, weight, kMom, pMom, v0);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00)
{
    int const centrality = floor(nCentHftRatio * gRandom->Rndm());
    
    TVector3 const vertex(0,0,0);
    // smear primary vertex
    // float const sigmaVertex = sigmaVertexCent[cent];
    // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));
    
    v00 += vertex;
    
    // smear momentum
    TLorentzVector kRMom;
    TLorentzVector pRMom;
    kRMom = smearMom(kMom, fKaonMomResolution);
    pRMom = smearMom(pMom, fPionMomResolution);
    
    // smear position
    TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, v00);
    TVector3 const pRPos = smearPosData(0, vertex.z(), centrality, pRMom, v00);
    
    if(pRMom.Perp()>1.6 && pRMom.Perp()>1.6) {
        TLorentzVector kMom_exc;
        TLorentzVector pMom_exc;
        kMom_exc.SetVectMag(kRMom.Vect(),pRMom.M());
        pMom_exc.SetVectMag(pRMom.Vect(),kRMom.M());
        fillDouble(kf, b, weight, kMom_exc, pMom_exc, v00, kRPos, pRPos);
    }
    
    // reconstruct
    TLorentzVector const rMom = kRMom + pRMom;
    float const kDca = dca(kMom.Vect(), v00, vertex);
    float const pDca = dca(pMom.Vect(), v00, vertex);
    float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
    float const kRSDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
    float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
    float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);
    float const pRDca = dca(pRMom.Vect(), pRPos, vertex);
    float const pRSDca = dcaSigned(pRMom.Vect(), pRPos, vertex);
    float const pRDcaXY = dcaXY(pRMom.Vect(), pRPos, vertex);
    float const pRDcaZ = dcaZ(pRMom.Vect(), pRPos, vertex);
    
    TVector3 v0;
    float const dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v0);
    float const decayLength = (v0 - vertex).Mag();
    float const dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
    float const D0DcaXY = dcaXY(rMom.Vect(), v0, vertex);
    float const D0DcaZ = dcaZ(rMom.Vect(), v0, vertex);
    float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
    float const angle12 = kRMom.Vect().Angle(pRMom.Vect());
    
    float kRPt = kRMom.Perp();
    float pRPt = pRMom.Perp();
    float kREta = kRMom.PseudoRapidity();
    float pREta = pRMom.PseudoRapidity();
    float D0RPt = rMom.Perp();
    float ktpc = (kRPt<0.2 ? 0 : fTpcK->Eval(kRPt));
    float pitpc = (pRPt<0.2 ? 0 : fTpcPi->Eval(pRPt));
    float ktof = (kRPt<1.6 ? fTofK->Eval(kRPt) : 1.0);
    float pitof = (pRPt<1.6 ? fTofK->Eval(pRPt) : 1.0);
    //float mDoubleCountWg = 1.0;
    float ptWg = 1;
    if(mMode == 1) ptWg = 1.2/1. * hD0PtWgFonll->GetBinContent(hD0PtWgFonll->FindBin(primaryPt)); //1.2 is due to |eta|<1.2,  6 is every event loop 6 D0/B
    else ptWg = 1.2/1. * hBPtWgFonll->GetBinContent(hD0PtWgFonll->FindBin(primaryPt));
    
    int const charge = kf > 0 ? 1 : -1;
    // save
    if(mUseTree) {
        float arr[110];
        int iArr = 0;
        arr[iArr++] = kf;
        arr[iArr++] = b->M();
        arr[iArr++] = b->Perp();
        arr[iArr++] = b->PseudoRapidity();
        arr[iArr++] = b->Rapidity();
        arr[iArr++] = b->Phi();
        arr[iArr++] = ptWg;
        arr[iArr++] = 1.0; //1.0 -- real signal, 0--double count
        
        arr[iArr++] = rMom.M();
        arr[iArr++] = rMom.Perp();
        arr[iArr++] = rMom.PseudoRapidity();
        arr[iArr++] = rMom.Rapidity();
        arr[iArr++] = rMom.Phi();
        
        arr[iArr++] = dca12;
        arr[iArr++] = decayLength;
        arr[iArr++] = dcaD0ToPv;
        arr[iArr++] = D0DcaXY;
        arr[iArr++] = D0DcaZ;
        arr[iArr++] = cosTheta;
        
        arr[iArr++] = kMom.M();
        arr[iArr++] = kMom.Perp();
        arr[iArr++] = kMom.PseudoRapidity();
        arr[iArr++] = kMom.Rapidity();
        arr[iArr++] = kMom.Phi();
        arr[iArr++] = kDca;
        
        arr[iArr++] = kRMom.M();
        arr[iArr++] = kRMom.Perp();
        arr[iArr++] = kRMom.PseudoRapidity();
        arr[iArr++] = kRMom.Rapidity();
        arr[iArr++] = kRMom.Phi();
        arr[iArr++] = kRDca;
        arr[iArr++] = kRSDca;
        arr[iArr++] = kRDcaXY;
        arr[iArr++] = kRDcaZ;
        arr[iArr++] = ktpc;
        arr[iArr++] = ktof;
        
        arr[iArr++] = pMom.M();
        arr[iArr++] = pMom.Perp();
        arr[iArr++] = pMom.PseudoRapidity();
        arr[iArr++] = pMom.Rapidity();
        arr[iArr++] = pMom.Phi();
        arr[iArr++] = pDca;
        
        arr[iArr++] = pRMom.M();
        arr[iArr++] = pRMom.Perp();
        arr[iArr++] = pRMom.PseudoRapidity();
        arr[iArr++] = pRMom.Rapidity();
        arr[iArr++] = pRMom.Phi();
        arr[iArr++] = pRDca;
        arr[iArr++] = pRSDca;
        arr[iArr++] = pRDcaXY;
        arr[iArr++] = pRDcaZ;
        arr[iArr++] = pitpc;
        arr[iArr++] = pitof;
        
        arr[iArr++] = mHftMatchRatio;
        arr[iArr++] = mHftMatchRatio;
        
        arr[iArr++] = primaryPt;
        arr[iArr++] = primaryY;
        arr[iArr++] = primaryPhi;
        
        nt_sig->Fill(arr);
        //if (ipart % 1000 == 1) nt_sig->AutoSave("SaveSelf");
    }
    if(mUseHist) {
        hpt->Fill(b->Perp());
        hptWg->Fill(b->Perp(),ptWg);
        hppt->Fill(primaryPt);
        hpptWg->Fill(primaryPt,ptWg);
        if(kf==421) {
            hD0->Fill("D^{0}",1);
            hD0Wg->Fill("D^{0}",ptWg);
        }
        else {
            hD0->Fill("#bar{D^{0}}",1);
            hD0Wg->Fill("#bar{D^{0}}",ptWg);
        }
        
        if(fabs(b->Rapidity())>=1.0) return;
        hD0pT_noCut->Fill(b->Perp(),ptWg);
        
        //cut
        float weight = mHftMatchRatio * mHftMatchRatio * ktof * pitof * ktpc * pitpc * ptWg;
        int ptIndex = getD0PtIndex(rMom.Perp());
        if(ptIndex < 0) return;
        bool isD0 = ( pRDca > anaCuts::pDca[ptIndex]
                     && kRDca > anaCuts::kDca[ptIndex]
                     && dca12 < anaCuts::dcaDaughters[ptIndex]
                     && decayLength > anaCuts::decayLength[ptIndex]
                     && cosTheta > anaCuts::cosTheta[ptIndex]
                     //&& dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex]
                     && kRPt > 0.6 && pRPt > 0.6
                     && fabs(kREta) < 1 && fabs(pREta) < 1
                     );
        bool isD0_prompt = isD0 && cosTheta>0.95 && dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex];
        if(isD0) {
            hD0pT_Cut->Fill(D0RPt,weight);
            hD0pT_Cut_doubleCount->Fill(D0RPt,weight);
            h2massSig->Fill(D0RPt,rMom.M(),weight);
            h2massSigDoubleCount->Fill(D0RPt,rMom.M(),weight);
            if(rMom.M()>1.845 && rMom.M()<1.885) { //about 2Sigma
                h2D0DcaSig->Fill(D0RPt,dcaD0ToPv/1.e4,weight);
                h2D0DcaSigDoubleCount->Fill(D0RPt,dcaD0ToPv/1.e4,weight);
            }
        }
        if(isD0_prompt) {
            hPromptD0pT_Cut->Fill(D0RPt,weight);
            h2massSig_D0->Fill(D0RPt,rMom.M(),weight);
        }
    }
}

void fillDouble(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00, TVector3 kPos, TVector3 pPos)
{
    int const centrality = floor(nCentHftRatio * gRandom->Rndm());
    
    TVector3 const vertex(0,0,0);
    // smear primary vertex
    // float const sigmaVertex = sigmaVertexCent[cent];
    // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));
    
    v00 += vertex;
    
    // smear momentum
    TLorentzVector kRMom;
    TLorentzVector pRMom;
    kRMom = kMom;//smearMom(kMom, fKaonMomResolution);
    pRMom = pMom;//smearMom(pMom, fPionMomResolution);
    
    // smear position
    TVector3 const kRPos = kPos;//smearPosData(1, vertex.z(), centrality, kRMom, v00);
    TVector3 const pRPos = pPos;//smearPosData(0, vertex.z(), centrality, pRMom, v00);
    
    // reconstruct
    TLorentzVector const rMom = kRMom + pRMom;
    float const kDca = dca(kMom.Vect(), v00, vertex);
    float const pDca = dca(pMom.Vect(), v00, vertex);
    float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
    float const kRSDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
    float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
    float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);
    float const pRDca = dca(pRMom.Vect(), pRPos, vertex);
    float const pRSDca = dcaSigned(pRMom.Vect(), pRPos, vertex);
    float const pRDcaXY = dcaXY(pRMom.Vect(), pRPos, vertex);
    float const pRDcaZ = dcaZ(pRMom.Vect(), pRPos, vertex);
    
    TVector3 v0;
    float const dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v0);
    float const decayLength = (v0 - vertex).Mag();
    float const dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
    float const D0DcaXY = dcaXY(rMom.Vect(), v0, vertex);
    float const D0DcaZ = dcaZ(rMom.Vect(), v0, vertex);
    float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
    float const angle12 = kRMom.Vect().Angle(pRMom.Vect());
    
    float kRPt = kRMom.Perp();
    float pRPt = pRMom.Perp();
    float kREta = kRMom.PseudoRapidity();
    float pREta = pRMom.PseudoRapidity();
    float D0RPt = rMom.Perp();
    float ktpc = (kRPt<0.2 ? 0 : fTpcK->Eval(kRPt));
    float pitpc = (pRPt<0.2 ? 0 : fTpcPi->Eval(pRPt));
    float ktof = (kRPt<1.6 ? fTofK->Eval(kRPt) : 1.0);
    float pitof = (pRPt<1.6 ? fTofK->Eval(pRPt) : 1.0);
    //float mDoubleCountWg = 1.0;
    float ptWg = 1;
    if(mMode == 1) ptWg = 1.2/1. * hD0PtWgFonll->GetBinContent(hD0PtWgFonll->FindBin(primaryPt)); //1.2 is due to |eta|<1.2,  6 is every event loop 6 D0/B
    else ptWg = 1.2/1. * hBPtWgFonll->GetBinContent(hD0PtWgFonll->FindBin(primaryPt));
    
    int const charge = kf > 0 ? 1 : -1;
    // save
    if(mUseTree) {
        float arr[110];
        int iArr = 0;
        arr[iArr++] = kf;
        arr[iArr++] = b->M();
        arr[iArr++] = b->Perp();
        arr[iArr++] = b->PseudoRapidity();
        arr[iArr++] = b->Rapidity();
        arr[iArr++] = b->Phi();
        arr[iArr++] = ptWg;
        arr[iArr++] = 0; //1.0 -- real signal, 0--double count
        
        arr[iArr++] = rMom.M();
        arr[iArr++] = rMom.Perp();
        arr[iArr++] = rMom.PseudoRapidity();
        arr[iArr++] = rMom.Rapidity();
        arr[iArr++] = rMom.Phi();
        
        arr[iArr++] = dca12;
        arr[iArr++] = decayLength;
        arr[iArr++] = dcaD0ToPv;
        arr[iArr++] = D0DcaXY;
        arr[iArr++] = D0DcaZ;
        arr[iArr++] = cosTheta;
        
        arr[iArr++] = kMom.M();
        arr[iArr++] = kMom.Perp();
        arr[iArr++] = kMom.PseudoRapidity();
        arr[iArr++] = kMom.Rapidity();
        arr[iArr++] = kMom.Phi();
        arr[iArr++] = kDca;
        
        arr[iArr++] = kRMom.M();
        arr[iArr++] = kRMom.Perp();
        arr[iArr++] = kRMom.PseudoRapidity();
        arr[iArr++] = kRMom.Rapidity();
        arr[iArr++] = kRMom.Phi();
        arr[iArr++] = kRDca;
        arr[iArr++] = kRSDca;
        arr[iArr++] = kRDcaXY;
        arr[iArr++] = kRDcaZ;
        arr[iArr++] = ktpc;
        arr[iArr++] = ktof;
        
        arr[iArr++] = pMom.M();
        arr[iArr++] = pMom.Perp();
        arr[iArr++] = pMom.PseudoRapidity();
        arr[iArr++] = pMom.Rapidity();
        arr[iArr++] = pMom.Phi();
        arr[iArr++] = pDca;
        
        arr[iArr++] = pRMom.M();
        arr[iArr++] = pRMom.Perp();
        arr[iArr++] = pRMom.PseudoRapidity();
        arr[iArr++] = pRMom.Rapidity();
        arr[iArr++] = pRMom.Phi();
        arr[iArr++] = pRDca;
        arr[iArr++] = pRSDca;
        arr[iArr++] = pRDcaXY;
        arr[iArr++] = pRDcaZ;
        arr[iArr++] = pitpc;
        arr[iArr++] = pitof;
        
        arr[iArr++] = mHftMatchRatio;
        arr[iArr++] = mHftMatchRatio;
        
        arr[iArr++] = primaryPt;
        arr[iArr++] = primaryY;
        arr[iArr++] = primaryPhi;
        
        nt_sig->Fill(arr);
        //if (ipart % 1000 == 1) nt_sig->AutoSave("SaveSelf");
    }
    if(mUseHist) {
        if(fabs(b->Rapidity())>=1.0) return;
        //cut
        float weight = mHftMatchRatio * mHftMatchRatio * ktof * pitof * ktpc * pitpc * ptWg;
        int ptIndex = getD0PtIndex(rMom.Perp());
        if(ptIndex < 0) return;
        bool isD0 = ( pRDca > anaCuts::pDca[ptIndex]
                     && kRDca > anaCuts::kDca[ptIndex]
                     && dca12 < anaCuts::dcaDaughters[ptIndex]
                     && decayLength > anaCuts::decayLength[ptIndex]
                     && cosTheta > anaCuts::cosTheta[ptIndex]
                     //&& dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex]
                     && kRPt > 0.6 && pRPt > 0.6
                     && fabs(kREta) < 1 && fabs(pREta) < 1
                     );
        bool isD0_prompt = isD0 && cosTheta>0.95 && dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex];
        if(isD0) {
            hD0pT_Cut_doubleCount->Fill(D0RPt,weight);
            h2massSigDoubleCount->Fill(D0RPt,rMom.M(),weight);
            if(rMom.M()>1.845 && rMom.M()<1.885) { //about 2Sigma
                h2D0DcaSigDoubleCount->Fill(D0RPt,dcaD0ToPv/1.e4,weight);
            }
        }
        if(isD0_prompt) {
            h2massSigDoubleCount_D0->Fill(D0RPt,rMom.M(),weight);
        }
    }
}

void getKinematics(TLorentzVector& b, double const mass)
{
    float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
    float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
    float const phi = TMath::TwoPi() * gRandom->Rndm();
    
    primaryPt = pt;
    primaryY = y;
    primaryPhi = phi;
    
    float const mT = sqrt(mass * mass + pt * pt);
    float const pz = mT * sinh(y);
    float const E = mT * cosh(y);
    
    b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
    
    return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 newPos(pos);
    newPos.SetZ(0);
    
    TVector3 newP(p);
    newP.SetZ(0);
    
    TVector3 newVertex(vertex);
    newVertex.SetZ(0);
    
    TVector3 posDiff = newPos - newVertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
    return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    if (sin(p.Theta()) == 0) return 0;
    else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
    TVector3 posDiff = pos2 - pos1;
    TVector3 pu1 = p1.Unit();
    TVector3 pu2 = p2.Unit();
    double pu1Pu2 = pu1.Dot(pu2);
    double g = posDiff.Dot(pu1);
    double k = posDiff.Dot(pu2);
    double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
    double s1 = g + s2 * pu1Pu2;
    TVector3 posDca1 = pos1 + pu1 * s1;
    TVector3 posDca2 = pos2 + pu2 * s2;
    v0 = 0.5 * (posDca1 + posDca2);
    return (posDca1 - posDca2).Mag();
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
    float const pt = b.Perp();
    float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));
    
    TLorentzVector sMom;
    sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
    return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
    float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
    float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
    float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));
    
    return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

int getPtIndexDca(double pT)
{
    for (int i = 0; i < nPtBinsDca; i++)
    {
        if ((pT >= ptEdgeDca[i]) && (pT < ptEdgeDca[i + 1]))
            return i;
    }
    return nPtBinsDca - 1 ;
}

int getEtaIndexDca(double Eta)
{
    for (int i = 0; i < nEtasDca; i++)
    {
        if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
            return i;
    }
    return nEtasDca - 1 ;
}

int getVzIndexDca(double Vz)
{
    for (int i = 0; i < nVzsDca; i++)
    {
        if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
            return i;
    }
    return nVzsDca - 1 ;
}


int getPtIndexHftRatio(double pT)
{
    for (int i = 0; i < nPtBinsHftRatio; i++)
    {
        if ((pT >= ptEdgeHftRatio[i]) && (pT < ptEdgeHftRatio[i + 1]))
            return i;
    }
    return nPtBinsHftRatio - 1 ;
}

int getEtaIndexHftRatio(double Eta)
{
    for (int i = 0; i < nEtasHftRatio; i++)
    {
        if ((Eta >= EtaEdgeHftRatio[i]) && (Eta < EtaEdgeHftRatio[i + 1]))
            return i;
    }
    return nEtasHftRatio - 1 ;
}

int getVzIndexHftRatio(double Vz)
{
    for (int i = 0; i < nVzsHftRatio; i++)
    {
        if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
            return i;
    }
    return nVzsHftRatio - 1 ;
}


TVector3 smearPosData(int const iParticleIndex, double const vz, int cent, TLorentzVector const& rMom, TVector3 const& pos)
{
    int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
    int const iVzIndex = getVzIndexDca(vz);
    // int const iPhiIndex = getPhiIndexDca(rMom.Phi());
    int const iPtIndex = getPtIndexDca(rMom.Perp());
    
    double sigmaPosZ = 0;
    double sigmaPosXY = 0;
    
    if (cent == 8) cent = 7;
    //All the centrality position smear was based on 0-10% centrality input, so here the cent==0
    
    // h2Dca[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom2(sigmaPosXY,sigmaPosZ);
    h2Dca[iParticleIndex][iPtIndex]->GetRandom2(sigmaPosXY, sigmaPosZ);
    sigmaPosZ *= 1.e4;
    sigmaPosXY *= 1.e4;
    /*if (h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
     {
     do sigmaPosZ = h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
     while (fabs(sigmaPosZ) > 1.e3);
     }
     
     if (h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
     {
     do sigmaPosXY = h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
     while (fabs(sigmaPosXY) > 1.e3);
     }
     */
    
    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
    newPos -= momPerp.Unit() * sigmaPosXY;
    
    return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 getVertex(int const centrality)
{
    double rdmVz=0;
    
    /*if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
     else
     {
     do rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
     while (fabs(rdmVz) > gVzCut);
     }*/
    
    return TVector3(0., 0., rdmVz);
}

//___________
void bookObjects()
{
    cout << "Loading input hft match ratio ..." << endl;
    
    cout << "Loading input Dca2D ..." << endl;
    TFile fDca2D("2DProjection_DcaXyZ_sPHENIX.root");
    for (int iParticle = 0; iParticle < nParticles; ++iParticle)
    {
        for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
        {
            h2Dca[iParticle][iPt] = (TH2F*)(fDca2D.Get(Form("mh2DcaPtPart_%i_%i", iParticle, iPt)));
            h2Dca[iParticle][iPt]->SetDirectory(0);
        }
    }
    fDca2D.Close();
    
    cout << "Loading momentum resolution ..." << endl;
    fPionMomResolution = new TF1("fPionMomResolution","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]*x*x)",0,20);
    fPionMomResolution->SetParameters(0.0120118, -1.01542e-09, 0.00169927);
    fKaonMomResolution = new TF1("fKaonMomResolution","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]*x*x)",0,20);
    fKaonMomResolution->SetParameters(0.0113136, 0.00672005, 0.00179494);
    fProtonMomResolution = new TF1("fProtonMomResolution","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]*x*x)",0,20);
    fProtonMomResolution->SetParameters(0.0101259, 0.0124398, 0.00188631);
    
    cout << "Loading k/pi from hijing sample ..." << endl;
    //pi/k/p number
    TFile fDaughters("PiKPNumber.root");
    h3PiKPNumber = (TH3I*)fDaughters.Get("h3PiKPNumber60_80");
    h3PiKPNumber->SetDirectory(0);
    fDaughters.Close();
    
    //pi/k/p pT weight
    TFile fDaughterPtWeight("input_DaughterPtWg.root"); //0-12GeV
    hpiPtWg = (TH1F*)fDaughterPtWeight.Get("hpiPtWg");
    hpiPtWg->SetDirectory(0);
    hkPtWg = (TH1F*)fDaughterPtWeight.Get("hkPtWg");
    hkPtWg->SetDirectory(0);
    hpPtWg = (TH1F*)fDaughterPtWeight.Get("hpPtWg");
    hpPtWg->SetDirectory(0);
    fDaughterPtWeight.Close();
    
    cout << "Loading Tpc tracking efficiency ..." << endl;
    fTpcK = new TF1("fTpcK","[0]*exp(-pow(x/[1],[2]))",0.2,20.);
    fTpcK->SetParameters(0.951716, 0.169177, -0.511248);
    fTpcPi = new TF1("fTpcPi","[0]*exp(-pow(x/[1],[2]))",0.2,20.);
    fTpcPi->SetParameters(1.06019, 0.0664062, -0.297934);
    fTpcP = new TF1("fTpcP","[0]*exp(-pow(x/[1],[2]))",0.2,20.);
    fTpcP->SetParameters(0.998745, 0.0812092, -0.370134);
    
    cout << "Loading Tof matching efficiency ..." << endl;
    fTofK = new TF1("fTofK","[0]/(pow(x-[1],2)+[2])-[4]/(exp(x-[3])+[5])+[6]",0.2,3);
    fTofK->SetParameters(-0.00377627*1.6689, 0.198926, 0.00805162, -9.01765, -0.653181*1.6689, -0.523133, 0.58466*1.6689);
    fTofPi = new TF1("fTofPi","[0]/(pow(x-[1],2)+[2])-[4]/(exp(x-[3])+[5])+[6]",0.2,3);
    fTofPi->SetParameters(-0.0743809*1.6689, -2.89051, -8.97966, -0.199471, 1.6689*(2.26781e-07), -1.69786, 0.602073*1.6689);
    fTofP = new TF1("fTofP","[0]/(pow(x-[1],2)+[2])-[4]/(exp(x-[3])+[5])+[6]",0.2,3);
    fTofP->SetParameters(-0.0110193*1.6689, 0.257745, 0.0203688, -0.196571, -(4.60251e-05)*1.6689, -1.70964, 0.596767*1.6689);
    
    cout << "Loading D0/B pT weight ..." << endl;
    TFile fD0PtWeight("D0Pt_Wg.root"); //0-12GeV
    hD0PtWgFonll = (TH1F*)fD0PtWeight.Get("hD0Pt_wg_60_80");
    hD0PtWgFonll->SetDirectory(0);
    hBPtWgFonll = (TH1F*)fD0PtWeight.Get("hBPt_wg_60_80");
    hBPtWgFonll->SetDirectory(0);
    fD0PtWeight.Close();
    
    cout << "Done with loading all files ..." << endl;
    
    result = new TFile(outFileName.Data(), "recreate");
    result->SetCompressionLevel(1);
    result->cd();
}
//___________
void write()
{
    result->cd();
    if(mUseHist) WriteHist(result);
    if(mUseTree) {
        nt_sig->Write();
        nt_bg->Write();
    }
    result->Close();
}
int getD0PtIndex(float const pt)
{
    int bin = -1;
    for (int i = 0; i < anaCuts::nPtBins; i++)
    {
        if ((pt >= anaCuts::PtEdge[i]) && (pt < anaCuts::PtEdge[i + 1]))
            bin = i;
    }
    return bin;
}
