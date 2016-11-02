/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
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

using namespace std;

void particlePionFill(int const kf, TLorentzVector b, TVector3 vertex, TVector3& rPionPos, TLorentzVector& rPionMom);
void particleKaonFill(int const kf, TLorentzVector b, TVector3 vertex, TVector3& rKaonPos, TLorentzVector& rKaonMom);
void particleProtonFill(int const kf, TLorentzVector b, TVector3 vertex, TVector3& rProtonPos, TLorentzVector& rProtonMom);
void particleKaonFillwoPID(TVector3 kPos, TLorentzVector kMom, TVector3& rKaonPos, TLorentzVector& rKaonMom);
void particlePionFillwoPID(TVector3 piPos, TLorentzVector piMom, TVector3& rPionPos, TLorentzVector& rPionMom);
void particleProtonFillwoPID(TVector3 pPos, TLorentzVector pMom, TVector3& rProtonPos, TLorentzVector& rProtonMom);
void rcD0Fill(TVector3 vertex, TVector3 rKaonPos, TLorentzVector rKaonMom, TVector3 rPionPos, TLorentzVector rPionMom);
void TmpnDaughtersFill(int const nkaons, int const npions, int const nprotons);
void getKinematics(TLorentzVector& b, double const mass);
void getKinematicsPionPlus(TLorentzVector& b, double const mass);
void getKinematicsPionMinus(TLorentzVector& b, double const mass);
void getKinematicsKaonPlus(TLorentzVector& b, double const mass);
void getKinematicsKaonMinus(TLorentzVector& b, double const mass);
void getKinematicsProtonPlus(TLorentzVector& b, double const mass);
void getKinematicsProtonMinus(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
// int getPhiIndexDca(double);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
// int getPhiIndexHftRatio(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TNtuple* TMPnt;
TNtuple* ntPion;
TNtuple* ntKaon;
TNtuple* ntProton;
TFile* result;
TH3F* h3daughters;

TF1* fPionMomResolution = NULL;
TF1* fKaonMomResolution = NULL; //currentlly use the same momentum resolution as PION
TF1* fProtonMomResolution = NULL;

TH1D* fpiplusPtHijing = NULL;
TH1D* fkplusPtHijing = NULL;
TH1D* fpplusPtHijing = NULL;

TH1F* fpiplusEtaHijing = NULL;
TH1F* fkplusEtaHijing = NULL;
TH1F* fpplusEtaHijing = NULL;

TH1F* fpiplusPhiHijing = NULL;
TH1F* fkplusPhiHijing = NULL;
TH1F* fpplusPhiHijing = NULL;

//Minus
TH1D* fpiminusPtHijing = NULL;
TH1D* fkminusPtHijing = NULL;
TH1D* fpminusPtHijing = NULL;

TH1F* fpiminusEtaHijing = NULL;
TH1F* fkminusEtaHijing = NULL;
TH1F* fpminusEtaHijing = NULL;

TH1F* fpiminusPhiHijing = NULL;
TH1F* fkminusPhiHijing = NULL;
TH1F* fpminusPhiHijing = NULL;

//
double piplus = 1.1479498;//e7
double piminus = 1.1774518;//e7
double kplus = 9.7553;//e5
double kminus = 8.68088;//e5
double pplus = 7.51646;//e5
double pminus = 5.6805;//e5

TH1F* hpions = NULL;
TH1F* hkaons = NULL;
TH1F* hprotons = NULL;
TH3F* h3kaonspionsprotons = NULL;
const Int_t nParticles = 3;
const Int_t nCentHftRatio = 1;

// HFT ratio binning
const Int_t nEtasHftRatio = 1;
const Int_t nVzsHftRatio = 1;
const Int_t nPtBinsHftRatio = 17;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] = { -1.0, 1.0 };
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = { -6.0e4,  6.0e4};
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] = { 0.3 ,  0.4 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 2.6 , 3.0 , 4. , 10. };

// DCA binning
int const nVzsDca = 1;
float const VzEdgeDca[nVzsDca + 1] = { -6.e4, 6.e4};

int const nEtasDca = 1;
float const EtaEdgeDca[nEtasDca + 1] = { -1.0, 1.0};

const Int_t nPtBinsDca = 17;
const Double_t ptEdgeDca[nPtBinsDca + 1] =
{ 0.3 ,  0.4 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 2.6 , 3.0 , 4. , 10. };

TH1D* h1Vz[nCentHftRatio];

// TH1D* hHftRatio1[nParticles][nEtasHftRatio][nVzsHftRatio][nCentHftRatio];
TF1*  fMatch; //HFT match ratio/ consistent as 0.91
int const nCentDca = 1;
// TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nPtBinsDca];
TH2D* h2Dca[nParticles][nPtBinsDca];

TH1D* hTpcPi[nCentHftRatio];
TH1D* hTpcK[nCentHftRatio];
TH1D* hTpcP[nCentHftRatio];

string outFileName = "SimulationD0BG.toyMc.root";
std::pair<float, float> const momentumRange(0, 10);

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCentHftRatio] = { 4.};

//============== main  program ==================
void toyD0BackGround(int nEvts = 1000)
{
   // TMemStat mem;
   // mem.Enable();
   TStopwatch*   stopWatch = new TStopwatch();
   stopWatch->Start();
   gRandom->SetSeed();
   bookObjects();

   for (int iEvt = 0; iEvt < nEvts; iEvt++)
   {
      if (!(iEvt % 1000))
         cout << "____________ iEvt= " << iEvt / static_cast<float>(nEvts) << " ________________" << endl;

      double tmpx, tmpy, tmpz; //Maybe proton kaon pion are related
      h3kaonspionsprotons->GetRandom3(tmpx, tmpy, tmpz);
      // cout<<tmpx << tmpy << tmpz <<endl;
      // cout<<(int) tmpx <<" "<<(int)(tmpx+0.5) <<endl;
      const int nkaonsAll     =  floor(tmpx + 0.5);
      const int npionsAll     =  floor(tmpy + 0.5);
      const int nprotonsAll   =  floor(tmpz + 0.5);
      // TClonesArray ptl("TParticle", 10);
      // cout<<"nprotons=" <<nprotons <<" ;nkaons=" <<nkaons <<" ;npions="<<npions<<endl;
      const int nkaonsPos = floor(nkaonsAll * kplus / (kplus + kminus));
      const int nkaonsNeg = nkaonsAll - nkaonsPos;
      const int npionsPos = floor(npionsAll * piplus / (piplus + piminus));
      const int npionsNeg = npionsAll - npionsPos;
      const int nprotonsPos = floor(nprotonsAll * pplus / (pplus + pminus));
      const int nprotonsNeg = nprotonsAll - nprotonsPos;

      const int n_Alls = nkaonsAll + npionsAll + nprotonsAll;

      int const centrality = floor(nCentHftRatio * gRandom->Rndm());
      TVector3 const vertex = getVertex(centrality);

      //Positive
      // TLorentzVector *b_pi[npions] = new TLorentzVector;
      TLorentzVector* b_pi_pos = new TLorentzVector[npionsPos] ;
      TLorentzVector* b_k_pos = new TLorentzVector[nkaonsPos] ;
      TLorentzVector* b_p_pos = new TLorentzVector[nprotonsPos] ;

      TLorentzVector* r_pi_pos = new TLorentzVector[npionsPos] ;
      TLorentzVector* r_k_pos  = new TLorentzVector[nkaonsPos] ;
      TLorentzVector* r_p_pos = new TLorentzVector[nprotonsPos] ;

      TVector3 r_KaonPos_pos[nkaonsPos];
      TVector3 r_PionPos_pos[npionsPos];
      TVector3 r_ProtonPos_pos[nprotonsPos];

      //Negative
      TLorentzVector* b_pi_neg = new TLorentzVector[npionsNeg] ;
      TLorentzVector* b_k_neg = new TLorentzVector[nkaonsNeg] ;
      TLorentzVector* b_p_neg = new TLorentzVector[nprotonsNeg] ;

      TLorentzVector* r_pi_neg = new TLorentzVector[npionsNeg] ;
      TLorentzVector* r_k_neg  = new TLorentzVector[nkaonsNeg] ;
      TLorentzVector* r_p_neg = new TLorentzVector[nprotonsNeg] ;

      TVector3 r_KaonPos_neg[nkaonsNeg];
      TVector3 r_PionPos_neg[npionsNeg];
      TVector3 r_ProtonPos_neg[nprotonsNeg];

      vector<int> Index;//
      Index.clear();
      int index = 0;
      vector<int> Charge;//
      Charge.clear();

      TVector3 r_KaonPos_All[n_Alls];
      TVector3 r_PionPos_All[n_Alls];

      TLorentzVector* b_pi_All = new TLorentzVector[n_Alls] ;
      TLorentzVector* b_k_All = new TLorentzVector[n_Alls] ;

      TLorentzVector* r_pi_All = new TLorentzVector[n_Alls] ;
      TLorentzVector* r_k_All  = new TLorentzVector[n_Alls] ;

      // 211 pi+   -321 K-   2212 proton
      // cout<< " n_Alls = " << n_Alls << " , nkaonsPos = " << nkaonsPos << " , nkaonsNeg = " << nkaonsNeg << " , npionsPos = " << npionsPos << " , npionsNeg = " << npionsNeg << " , nprotonsPos = " << nprotonsPos << " , nprotonsNeg = " << nprotonsNeg << endl;

      //Fill kaons
      for (int ik = 0; ik < nkaonsPos; ik++)
      {
         getKinematicsKaonPlus(b_k_pos[ik],  M_KAON_PLUS);
         particleKaonFill(321, b_k_pos[ik],  vertex, r_KaonPos_pos[ik], r_k_pos[ik]);

         Index.push_back(index);
         Charge.push_back(1);

         particleKaonFillwoPID(r_KaonPos_pos[ik], r_k_pos[ik], r_KaonPos_All[index], r_k_All[index]);
         particlePionFillwoPID(r_KaonPos_pos[ik], r_k_pos[ik], r_PionPos_All[index], r_pi_All[index]);

         index++;
      }
      // cout<< "after nkaonsPos, index = " << index <<endl;

      for (int ik = 0; ik < nkaonsNeg; ik++)
      {
         getKinematicsKaonMinus(b_k_neg[ik],  M_KAON_MINUS);
         particleKaonFill(-321, b_k_neg[ik],  vertex, r_KaonPos_neg[ik], r_k_neg[ik]);

         Index.push_back(index);
         Charge.push_back(-1);

         particleKaonFillwoPID(r_KaonPos_neg[ik], r_k_neg[ik], r_KaonPos_All[index], r_k_All[index]);
         particlePionFillwoPID(r_KaonPos_neg[ik], r_k_neg[ik], r_PionPos_All[index], r_pi_All[index]);

         index++;
      }
      // cout<< "after nkaonsNeg, index = " << index <<endl;


      //Fill Pion
      for (int ipi = 0; ipi < npionsPos; ipi++)
      {
         getKinematicsPionPlus(b_pi_pos[ipi],  M_PION_PLUS);
         particlePionFill(211, b_pi_pos[ipi],  vertex, r_PionPos_pos[ipi], r_pi_pos[ipi]);

         Index.push_back(index);
         Charge.push_back(1);

         particleKaonFillwoPID(r_PionPos_pos[ipi], r_pi_pos[ipi], r_KaonPos_All[index], r_k_All[index]);
         particlePionFillwoPID(r_PionPos_pos[ipi], r_pi_pos[ipi], r_PionPos_All[index], r_pi_All[index]);

         index++;
      }
      // cout<< "after npionsPos, index = " << index <<endl;

      for (int ipi = 0; ipi < npionsNeg; ipi++)
      {
         getKinematicsPionMinus(b_pi_neg[ipi],  M_PION_MINUS);
         particlePionFill(-211, b_pi_neg[ipi],  vertex, r_PionPos_neg[ipi], r_pi_neg[ipi]);

         Index.push_back(index);
         Charge.push_back(-1);

         particleKaonFillwoPID(r_PionPos_neg[ipi], r_pi_neg[ipi], r_KaonPos_All[index], r_k_All[index]);
         particlePionFillwoPID(r_PionPos_neg[ipi], r_pi_neg[ipi], r_PionPos_All[index], r_pi_All[index]);

         index++;
      }
      // cout<< "after npionsNeg, index = " << index <<endl;

      //Fill protons
      for (int ip = 0; ip < nprotonsPos; ip++)
      {
         getKinematicsProtonPlus(b_p_pos[ip],  M_PROTON);
         particleProtonFill(2212, b_p_pos[ip],  vertex, r_ProtonPos_pos[ip], r_p_pos[ip]);

         Index.push_back(index);
         Charge.push_back(1);

         particleKaonFillwoPID(r_ProtonPos_pos[ip], r_p_pos[ip], r_KaonPos_All[index], r_k_All[index]);
         particlePionFillwoPID(r_ProtonPos_pos[ip], r_p_pos[ip], r_PionPos_All[index], r_pi_All[index]);

         index++;
      }
      // cout<< "after nprotonsPos, index = " << index <<endl;

      for (int ip = 0; ip < nprotonsNeg; ip++)
      {
         getKinematicsProtonMinus(b_p_neg[ip],  M_ANTIPROTON);
         particleProtonFill(-2212, b_p_neg[ip],  vertex, r_ProtonPos_neg[ip], r_p_neg[ip]);

         Index.push_back(index);
         Charge.push_back(-1);

         particleKaonFillwoPID(r_ProtonPos_neg[ip], r_p_neg[ip], r_KaonPos_All[index], r_k_All[index]);
         particlePionFillwoPID(r_ProtonPos_neg[ip], r_p_neg[ip], r_PionPos_All[index], r_pi_All[index]);

         index++;
      }
      // cout<< "after nprotonsNeg, index = " << index <<endl;

      TmpnDaughtersFill(nkaonsAll, npionsAll, nprotonsAll);

      // cout<<"HERE D is OK!"<< endl;
      // Next is for reconstruct
      int iSelf = 0;
      int iTmp = 0;
      for (int ik = 0; ik < n_Alls; ik++)
      {
         if (r_k_All[ik].Perp() < 0.3) continue;

         if (fabs(dcaSigned(r_k_All[ik].Vect(), r_KaonPos_All[ik], vertex)) < 55) continue;

         for (int ipi = 0; ipi < n_Alls; ipi++)
         {
            if (r_pi_All[ipi].Perp() < 0.3) continue;

            if (fabs(dcaSigned(r_pi_All[ipi].Vect(), r_PionPos_All[ipi], vertex)) < 55) continue;

            // if(fabs(r_k_All[ik].Perp()-r_pi_All[ipi].Perp())<1e-6 && fabs(r_k_All[ik].Eta()-r_pi_All[ipi].Eta())<1e-6 && fabs(r_k_All[ik].Phi()-r_pi_All[ipi].Phi())<1e-6 ) //tricky way to avoid missed-particle make pair with himself
            if (ik == ipi) //same condition as previous line
            {
               iSelf++;
               continue;
            }

            if (Charge[ik] != Charge[ipi]) continue; //remove unlike sign// only save same sign for BG

            rcD0Fill(vertex, r_KaonPos_All[ik], r_k_All[ik], r_PionPos_All[ipi], r_pi_All[ipi]);
            iTmp++;

         }
      }

      // cout <<"npionsAll = " << npionsAll << " , nkaonsAll = " << nkaonsAll << " , nprotonsAll = " << nprotonsAll << " , iSelf = " << iSelf << " , iTmp = " << iTmp << endl;
      // cout<<"HERE E is OK!"<< endl;


      if (iEvt % 1000 == 1)
      {
         // ntPion->AutoSave("SaveSelf");
         // ntKaon->AutoSave("SaveSelf");
         // ntProton->AutoSave("SaveSelf");
         nt->AutoSave("SaveSelf");
         TMPnt->AutoSave("SaveSelf");
      }
      // cout<<"HERE F is OK!"<< endl;

      delete[] b_pi_pos;
      delete[] r_pi_pos;
      delete[] b_k_pos;
      delete[] r_k_pos;
      delete[] b_p_pos;
      delete[] r_p_pos;

      delete[] b_pi_neg;
      delete[] r_pi_neg;
      delete[] b_k_neg;
      delete[] r_k_neg;
      delete[] b_p_neg;
      delete[] r_p_neg;

      delete[] b_pi_All;
      delete[] r_pi_All;
      delete[] b_k_All;
      delete[] r_k_All;
      // cout<<"HERE G is OK!"<< endl;
   }

   write();

// mem.Show();
   stopWatch->Stop();
   stopWatch->Print();
}

void particlePionFill(int const kf, TLorentzVector piMom, TVector3 vertex, TVector3& rPionPos, TLorentzVector& rPionMom)
{

   // smear momentum
   TLorentzVector const piRMom = smearMom(piMom, fPionMomResolution);

   int const centrality = floor(nCentHftRatio * gRandom->Rndm());//centrality is always 0 here
   // smear position
   TVector3 const piRPos = smearPosData(0, vertex.z(), centrality, piRMom, vertex);

   // reconstruct
   float const piDca = dca(piMom.Vect(), vertex, vertex);
   float const piPt  = piMom.Perp();
   float const piRPt = piRMom.Perp();
   float const piREta = piRMom.Eta();
   float const piRPhi = piRMom.Phi();
   float const piRDca = dcaSigned(piRMom.Vect(), piRPos, vertex);
   float const piRDcaXY = dcaXY(piRMom.Vect(), piRPos, vertex);
   float const piRDcaZ = dcaZ(piRMom.Vect(), piRPos, vertex);

   rPionPos = piRPos;
   rPionMom = piRMom;

   // save
   // float arr[110];
   // int iArr = 0;
   // arr[iArr++] = piDca;
   // arr[iArr++] = piPt;
   // arr[iArr++] = piRPt;
   // arr[iArr++] = piREta;
   // arr[iArr++] = piRPhi;
   // arr[iArr++] = piRDca;
   // arr[iArr++] = piRDcaXY;
   // arr[iArr++] = piRDcaZ;
   // arr[iArr++] = kf;
   // arr[iArr++] = tpcReconstructed(0, 1 , 0, piRMom);
   //
   // ntPion->Fill(arr);
}

void particleProtonFill(int const kf, TLorentzVector pMom, TVector3 vertex, TVector3& rProtonPos, TLorentzVector& rProtonMom)
{

   // smear momentum
   TLorentzVector const pRMom = smearMom(pMom, fProtonMomResolution);

   int const centrality = floor(nCentHftRatio * gRandom->Rndm());//centrality is always 0 here
   // smear position
   TVector3 const pRPos = smearPosData(0, vertex.z(), centrality, pRMom, vertex);

   // reconstruct
   float const pDca = dca(pMom.Vect(), vertex, vertex);
   float const pPt  = pMom.Perp();
   float const pRPt = pRMom.Perp();
   float const pREta = pRMom.Eta();
   float const pRPhi = pRMom.Phi();
   float const pRDca = dcaSigned(pRMom.Vect(), pRPos, vertex);
   float const pRDcaXY = dcaXY(pRMom.Vect(), pRPos, vertex);
   float const pRDcaZ = dcaZ(pRMom.Vect(), pRPos, vertex);

   rProtonPos = pRPos;
   rProtonMom = pRMom;

   // save
   // float arr[110];
   // int iArr = 0;
   // arr[iArr++] = pDca;
   // arr[iArr++] = pPt;
   // arr[iArr++] = pRPt;
   // arr[iArr++] = pREta;
   // arr[iArr++] = pRPhi;
   // arr[iArr++] = pRDca;
   // arr[iArr++] = pRDcaXY;
   // arr[iArr++] = pRDcaZ;
   // arr[iArr++] = kf;
   // arr[iArr++] = tpcReconstructed(0, 1 , 0, pRMom);
   //
   // ntProton->Fill(arr);
}

void particleKaonFill(int const kf, TLorentzVector kMom, TVector3 vertex, TVector3& rKaonPos, TLorentzVector& rKaonMom)
{

   // smear momentum
   TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution);

   int const centrality = floor(nCentHftRatio * gRandom->Rndm());
   // smear position
   TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, vertex);

   // reconstruct
   // this is only for kind of test, no need to produce for large sample
   float const kDca = dca(kMom.Vect(), vertex, vertex);
   float const kPt  = kMom.Perp();
   float const kRPt = kRMom.Perp();
   float const kREta = kRMom.Eta();
   float const kRPhi = kRMom.Phi();
   float const kRDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
   float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
   float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);

   rKaonPos = kRPos;
   rKaonMom = kRMom;

   // save
   // float arr[110];
   // int iArr = 0;
   // arr[iArr++] = kDca;
   // arr[iArr++] = kPt;
   // arr[iArr++] = kRPt;
   // arr[iArr++] = kREta;
   // arr[iArr++] = kRPhi;
   // arr[iArr++] = kRDca;
   // arr[iArr++] = kRDcaXY;
   // arr[iArr++] = kRDcaZ;
   // arr[iArr++] = kf;
   // arr[iArr++] = tpcReconstructed(1, -1 , 0, kRMom);
   //
   // ntKaon->Fill(arr);
}

//_______________________________ Kaon sample
void particleKaonFillwoPID(TVector3 kPos, TLorentzVector kMom, TVector3& rKaonPos, TLorentzVector& rKaonMom)
{
   rKaonPos = kPos;
   rKaonMom.SetPtEtaPhiM(kMom.Perp(), kMom.Eta() , kMom.Phi(), M_KAON_MINUS);
}
//_______________________________ Pion sample 2
void particlePionFillwoPID(TVector3 piPos, TLorentzVector piMom, TVector3& rPionPos, TLorentzVector& rPionMom)
{
   rPionPos = piPos;
   rPionMom.SetPtEtaPhiM(piMom.Perp(), piMom.Eta() , piMom.Phi(), M_PION_MINUS);
}
//_______________________________ Pion sample 2
void particleProtonFillwoPID(TVector3 pPos, TLorentzVector pMom, TVector3& rProtonPos, TLorentzVector& rProtonMom)
{
   rProtonPos = pPos;
   rProtonMom.SetPtEtaPhiM(pMom.Perp(), pMom.Eta() , pMom.Phi(), M_PROTON);
}
//_______________________________ RC D0
void rcD0Fill(TVector3 vertex, TVector3 rKaonPos, TLorentzVector rKaonMom, TVector3 rPionPos, TLorentzVector rPionMom)
{
   // reconstruct D0
   TLorentzVector const rMom =  rKaonMom + rPionMom;
   // if (rMom.Perp() < 0.99) return;//default cut is 1GeV/c
   if (rMom.M() > 2.1 || rMom.M() < 1.7) return;

   float const kRSDca = dcaSigned(rKaonMom.Vect(), rKaonPos, vertex);
   // float const kRDcaXY = dcaXY(rKaonMom.Vect(), rKaonPos, vertex);
   // float const kRDcaZ = dcaZ(rKaonMom.Vect(), rKaonPos, vertex);
   float const piRSDca = dcaSigned(rPionMom.Vect(), rPionPos, vertex);
   // float const piRDcaXY = dcaXY(rPionMom.Vect(), rPionPos, vertex);
   // float const piRDcaZ = dcaZ(rPionMom.Vect(), rPionPos, vertex);
   if (fabs(kRSDca) < 55 || fabs(piRSDca) < 55) return;

   TVector3 v12, v0;
   float const dca12 = dca1To2(rKaonMom.Vect(), rKaonPos, rPionMom.Vect(), rPionPos, v12);
   if (dca12 > 85.) return;

   float dcaDaughters = dca12 ;
   if (dcaDaughters > 85.) return;

   v0 = v12;
   float const decayLength = (v0 - vertex).Mag();
   if (decayLength < 140.) return;
   float const dcaToPv = dca(rMom.Vect(), v0, vertex);
   // if (dcaToPv > 65.) return;//this cut was removed for b->D decay
   float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
   if (cosTheta < 0.95) return; //default low limit cut is 0.97

   TLorentzVector kRMomRest = rKaonMom;
   TVector3 beta;
   beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
   kRMomRest.Boost(-beta);
   float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());

   // save
   float arr[110];
   int iArr = 0;
   arr[iArr++] = rKaonMom.Perp();
   arr[iArr++] = rKaonMom.PseudoRapidity();
   arr[iArr++] = tpcReconstructed(1, -1 , 0, rKaonMom);
   arr[iArr++] = rPionMom.Perp();
   arr[iArr++] = rPionMom.PseudoRapidity();
   arr[iArr++] = tpcReconstructed(0, 1 , 0, rPionMom);

   arr[iArr++] = rMom.M();
   arr[iArr++] = rMom.Perp();
   arr[iArr++] = rMom.Rapidity();
   arr[iArr++] = cosTheta;
   arr[iArr++] = decayLength;
   arr[iArr++] = dcaToPv;
   arr[iArr++] = kRSDca;
   arr[iArr++] = piRSDca;
   arr[iArr++] = dcaDaughters;
   arr[iArr++] = cosThetaStar;

   arr[iArr++] = matchHft(1, vertex.z(), 0, rKaonMom);
   arr[iArr++] = matchHft(0, vertex.z(), 0, rPionMom);

   nt->Fill(arr);
}

void TmpnDaughtersFill(int const nkaons, int const npions, int const nprotons)
{
   // save
   float arr[110];
   int iArr = 0;
   arr[iArr++] = nkaons;
   arr[iArr++] = npions;
   arr[iArr++] = nprotons;

   TMPnt->Fill(arr);
   h3daughters->Fill(nkaons, npions, nprotons);
}

void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   float const phi = TMath::TwoPi() * gRandom->Rndm();

   float const mT = sqrt(mass * mass + pt * pt);
   float const pz = mT * sinh(y);
   float const E = mT * cosh(y);

   // b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
   b.SetPtEtaPhiM(pt, y , phi, mass);
}

void getKinematicsPionPlus(TLorentzVector& b, double const mass)
{
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   // float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const phi = TMath::TwoPi() * gRandom->Rndm();
   float const pt = fpiplusPtHijing->GetRandom();
   float const y = fpiplusEtaHijing->GetRandom();
   float const phi = fpiplusPhiHijing->GetRandom();

   b.SetPtEtaPhiM(pt, y , phi, mass);
}

void getKinematicsPionMinus(TLorentzVector& b, double const mass)
{
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   // float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const phi = TMath::TwoPi() * gRandom->Rndm();
   float const pt = fpiminusPtHijing->GetRandom();
   float const y = fpiminusEtaHijing->GetRandom();
   float const phi = fpiminusPhiHijing->GetRandom();

   b.SetPtEtaPhiM(pt, y , phi, mass);
}

void getKinematicsKaonPlus(TLorentzVector& b, double const mass)
{
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   // float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const phi = TMath::TwoPi() * gRandom->Rndm();
   float const pt = fkplusPtHijing->GetRandom();
   float const y = fkplusEtaHijing->GetRandom();
   float const phi = fkplusPhiHijing->GetRandom();

   b.SetPtEtaPhiM(pt, y , phi, mass);
}

void getKinematicsKaonMinus(TLorentzVector& b, double const mass)
{
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   // float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const phi = TMath::TwoPi() * gRandom->Rndm();
   float const pt = fkminusPtHijing->GetRandom();
   float const y = fkminusEtaHijing->GetRandom();
   float const phi = fkminusPhiHijing->GetRandom();

   b.SetPtEtaPhiM(pt, y , phi, mass);
}

void getKinematicsProtonPlus(TLorentzVector& b, double const mass)
{
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   // float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const phi = TMath::TwoPi() * gRandom->Rndm();
   float const pt = fpplusPtHijing->GetRandom();
   float const y = fpplusEtaHijing->GetRandom();
   float const phi = fpplusPhiHijing->GetRandom();

   b.SetPtEtaPhiM(pt, y , phi, mass);
}

void getKinematicsProtonMinus(TLorentzVector& b, double const mass)
{
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   // float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const phi = TMath::TwoPi() * gRandom->Rndm();
   float const pt = fpminusPtHijing->GetRandom();
   float const y = fpminusEtaHijing->GetRandom();
   float const phi = fpminusPhiHijing->GetRandom();

   b.SetPtEtaPhiM(pt, y , phi, mass);
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
   double rdmVz;

   if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
   else
   {
      do rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
      while (fabs(rdmVz) > gVzCut);
   }

   return TVector3(0., 0., rdmVz);
}

bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom)
{
   TH1D* h = NULL;

   if (iParticleIndex == 0)
   {
      h = hTpcPi[cent];
   }
   else if (iParticleIndex == 1)
   {
      h = hTpcK[cent];
   }
   else if (iParticleIndex == 2)
   {
      h = hTpcP[cent];
   }

   // ...

   int const bin = h->FindBin(mom.Perp());

   return gRandom->Rndm() < h->GetBinContent(bin);
}

bool matchHft(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& mom)
{
   int const iEtaIndex = getEtaIndexHftRatio(mom.PseudoRapidity());
   int const iVzIndex = getVzIndexHftRatio(vz);
   // int const iPhiIndex = getPhiIndexHftRatio(mom.Phi());

   // int const bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][cent]->FindBin(mom.Perp());
   // return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][cent]->GetBinContent(bin);
   return gRandom->Rndm() < fMatch->Eval(mom.Perp());
}
//___________
void bookObjects()
{
   cout << "Loading input momentum resolution ..." << endl;

   TFile f("Function_Momentum_resolution_From_sPHENIX.root");
   fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPionPlus");
   fKaonMomResolution = (TF1*)f.Get("fPion")->Clone("fKaonMinus");
   fProtonMomResolution = (TF1*)f.Get("fPion")->Clone("fProtonPlus");
   f.Close();

   cout << "Loading input spectra ..." << endl;

   TFile fpikpHijing("hPionKaonProtonSpectra_D0Sample_HijingPlusFunction2.root");
   fpiplusPtHijing = (TH1D*)fpikpHijing.Get("mhpionsplusSpectraNew")->Clone("mhpionsplusSpectraNew");
   fpiplusPtHijing->SetDirectory(0);
   fkplusPtHijing = (TH1D*)fpikpHijing.Get("mhkaonsplusSpectraNew")->Clone("mhkaonsplusSpectraNew");
   fkplusPtHijing->SetDirectory(0);
   fpplusPtHijing = (TH1D*)fpikpHijing.Get("mhprotonsplusSpectraNew")->Clone("mhprotonsplusSpectraNew");
   fpplusPtHijing->SetDirectory(0);

   fpiplusEtaHijing = (TH1F*)fpikpHijing.Get("mhpionsplusEta")->Clone("mhpionsplusEta");
   fpiplusEtaHijing->SetDirectory(0);
   fkplusEtaHijing = (TH1F*)fpikpHijing.Get("mhkaonsplusEta")->Clone("mhkaonsplusEta");
   fkplusEtaHijing->SetDirectory(0);
   fpplusEtaHijing = (TH1F*)fpikpHijing.Get("mhprotonsplusEta")->Clone("mhprotonsplusEta");
   fpplusEtaHijing->SetDirectory(0);

   fpiplusPhiHijing = (TH1F*)fpikpHijing.Get("mhpionsplusPhi")->Clone("mhpionsplusPhi");
   fpiplusPhiHijing->SetDirectory(0);
   fkplusPhiHijing = (TH1F*)fpikpHijing.Get("mhkaonsplusPhi")->Clone("mhkaonsplusPhi");
   fkplusPhiHijing->SetDirectory(0);
   fpplusPhiHijing = (TH1F*)fpikpHijing.Get("mhprotonsplusPhi")->Clone("mhprotonsplusPhi");
   fpplusPhiHijing->SetDirectory(0);

   //minus
   fpiminusPtHijing = (TH1D*)fpikpHijing.Get("mhpionsminusSpectraNew")->Clone("mhpionsminusSpectraNew");
   fpiminusPtHijing->SetDirectory(0);
   fkminusPtHijing = (TH1D*)fpikpHijing.Get("mhkaonsminusSpectraNew")->Clone("mhkaonsminusSpectraNew");
   fkminusPtHijing->SetDirectory(0);
   fpminusPtHijing = (TH1D*)fpikpHijing.Get("mhprotonsminusSpectraNew")->Clone("mhprotonsminusSpectraNew");
   fpminusPtHijing->SetDirectory(0);

   fpiminusEtaHijing = (TH1F*)fpikpHijing.Get("mhpionsminusEta")->Clone("mhpionsminusEta");
   fpiminusEtaHijing->SetDirectory(0);
   fkminusEtaHijing = (TH1F*)fpikpHijing.Get("mhkaonsminusEta")->Clone("mhkaonsminusEta");
   fkminusEtaHijing->SetDirectory(0);
   fpminusEtaHijing = (TH1F*)fpikpHijing.Get("mhprotonsminusEta")->Clone("mhprotonsminusEta");
   fpminusEtaHijing->SetDirectory(0);

   fpiminusPhiHijing = (TH1F*)fpikpHijing.Get("mhpionsminusPhi")->Clone("mhpionsminusPhi");
   fpiminusPhiHijing->SetDirectory(0);
   fkminusPhiHijing = (TH1F*)fpikpHijing.Get("mhkaonsminusPhi")->Clone("mhkaonsminusPhi");
   fkminusPhiHijing->SetDirectory(0);
   fpminusPhiHijing = (TH1F*)fpikpHijing.Get("mhprotonsminusPhi")->Clone("mhprotonsminusPhi");
   fpminusPhiHijing->SetDirectory(0);

   fpikpHijing.Close();

   cout << "Loading input particle numbers ..." << endl;

   TFile fDaughters("hDaughters_D0Sample.root");
   hpions = (TH1F*)fDaughters.Get("mhpions")->Clone("hpions");
   hpions->SetDirectory(0);
   // cout <<"__TTT__"<< hpions->GetEntries() << endl;
   hkaons = (TH1F*)fDaughters.Get("mhkaons")->Clone("hkaons");
   hkaons->SetDirectory(0);
   hprotons = (TH1F*)fDaughters.Get("mhprotons")->Clone("hprotons");
   hprotons->SetDirectory(0);
   h3kaonspionsprotons = (TH3F*)fDaughters.Get("mh3kaonspionsprotons")->Clone("h3kaonspionsprotons");
   h3kaonspionsprotons->SetDirectory(0);
   fDaughters.Close();

   cout << "Loading Vz distributions ..." << endl;

   TFile fVertex("hVz_D0Sample.root");

   for (int ii = 0; ii < nCentHftRatio; ++ii)
   {
      h1Vz[ii]      = (TH1D*)fVertex.Get(Form("mh1Vz"))->Clone(Form("mh1Vz_%d", ii));
      h1Vz[ii]->SetDirectory(0);
   }

   fVertex.Close();

   cout << "Loading input HFT ratios and DCA ..." << endl;
   TFile fHftRatio1("Function_HFTMatch.root");
   fMatch = (TF1*)fHftRatio1.Get("fMatch")->Clone("fMatch");
   fHftRatio1.Close();

   TFile fDca1("2DProjection_DcaXyZ_16Oct30.root");

   for (int iParticle = 0; iParticle < nParticles; ++iParticle)
   {
      // DCA
      for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
      {
         h2Dca[iParticle][iPt] = (TH2D*)(fDca1.Get(Form("mh2DcaPtPart_%i_%i", iParticle, iPt)));
         h2Dca[iParticle][iPt]->SetDirectory(0);
      }
      cout << "Finished loading Dca: " <<  endl;
   }

   fDca1.Close();

   cout << " Loading TPC tracking efficiencies " << endl;

   TFile fTpcPi("Eff_Pion_embedding_FromHijing_1StopAll_16Oct30.root");
   TFile fTpcK("Eff_Kaon_embedding_FromHijing_1StopAll_16Oct30.root");
   TFile fTpcP("Eff_Proton_embedding_FromHijing_1StopAll_16Oct30.root");

   for (int iCent = 0; iCent < nCentHftRatio; ++iCent)
   {
      hTpcPi[iCent] = (TH1D*)fTpcPi.Get(Form("h1Ratiocent_%i", iCent))->Clone(Form("h1PiRatiocent_%i", iCent));
      hTpcPi[iCent]->SetDirectory(0);
      hTpcK[iCent] = (TH1D*)fTpcK.Get(Form("h1Ratiocent_%i", iCent))->Clone(Form("h1KRatiocent_%i", iCent));
      hTpcK[iCent]->SetDirectory(0);
      hTpcP[iCent] = (TH1D*)fTpcP.Get(Form("h1Ratiocent_%i", iCent))->Clone(Form("h1PRatiocent_%i", iCent));
      hTpcP[iCent]->SetDirectory(0);
   }

   fTpcPi.Close();
   fTpcK.Close();
   fTpcP.Close();

   cout << "Done with loading all files ..." << endl;

   result = new TFile(outFileName.c_str(), "recreate");
   result->SetCompressionLevel(1);
   result->cd();

   int BufSize = (int)pow(2., 16.);

   nt = new TNtuple("nt", "", "kPt:kEta:kTpc:piPt:piEta:piTpc:rD0M:rD0Pt:rD0Y:cosTheta:decayLength:dcaToPv:kDca:piDca:dcaDaughters:cosThetaStar:kHft:piHft", BufSize); //RC LC
   // nt->SetAutoSave(-500000); // autosave every 1 Mbytes
   TMPnt = new TNtuple("TMPnt", "", "nkaons:npions:nprotons"); //tmp ntuple for pion kaon protons numbers
   h3daughters = new TH3F("h3daughters", "", 100, 0, 100, 1000, 0, 1000, 100, 0, 100);

   ntPion = new TNtuple("ntPion", "", "piDca:piPt:piRPt:piREta:piRPhi:piRDca:piRDcaXY:piRDcaZ:pikf:piTpc"); // Rc Pion1
   ntKaon = new TNtuple("ntKaon", "", "kDca:kPt:kRPt:kREta:kRPhi:kRDca:kRDcaXY:kRDcaZ:kkf:kTpc"); // Rc Kaon
   ntProton = new TNtuple("ntProton", "", "pDca:pPt:pRPt:pREta:pRPhi:pRDca:pRDcaXY:pRDcaZ:pkf:piTpc"); // Rc Pion1
}
//___________
void write()
{
   result->cd();
   // ntPion->Write();
   // ntKaon->Write();
   // ntProton->Write();
   nt->Write();
   TMPnt->Write();
   h3daughters->Write();
   result->Close();
}
