/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *
 * RUN the CODE:
 * root -b -l -q toyMcEffZeroDecayLength.C++
 *
 * *********************************************************************
*/

#include <iostream>
#include <fstream>

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

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
// float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& posDca1, TVector3& posDca2, TVector3& v0);
TVector3 getVertex(int centrality);
bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
bool reconstructD0(int const centrality, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndex(double);
int getEtaIndex(double);
int getVzIndex(double);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
// int getPhiIndexHftRatio(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fKaonMomResolution = NULL;
TF1* fPionMomResolution = NULL;
TF1* fWeightFunction = NULL;
TGraph* grEff[3];
const Int_t nParticles = 2;
const Int_t nCent = 1;

// HFT ratio binning
const Int_t nEtasHftRatio = 1;
const Int_t nVzsHftRatio = 1;
const Int_t nPtBinsHftRatio = 17;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] = { -1.0,  1.0 };
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = { -6.0e4, 6.0e4};
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] = 
  { 0.3 ,  0.4 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 2.6 , 3.0 , 4. , 10. };
// const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] = {-3.14159 , -2.51327 , -1.88496 , -1.25664 , -0.628319 , 0 , 0.628319 , 1.25664 , 1.88496 , 2.51327 , 3.14159};

// DCA binning
int const nVzs = 1;
float const VzEdge[nVzs + 1] = { -6.e4,  6.e4};

int const nEtas = 1;
float const EtaEdge[nEtas + 1] = { -1.0, 1.0};

const Int_t nPtBins = 17;
const Double_t ptEdge[nPtBins + 1] =
  { 0.3 ,  0.4 , 0.5 , 0.6 , 0.7,  0.8 , 0.9 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2. , 2.4 , 2.6 , 3.0 , 4. , 10. };

TH1D* h1Vz[nCent];

// TH1D* hHftRatio1[nParticles][nEtasHftRatio][nVzsHftRatio][nCent];
TF1*  fMatch; //HFT match ratio/ consistent as 0.91
// TH1D* h1DcaZ1[nParticles][nEtas][nVzs][nCent][nPtBins];
// TH1D* h1DcaXY1[nParticles][nEtas][nVzs][nCent][nPtBins];
int const nCentDca = 1;
// TH2D* h2Dca[nParticles][nEtas][nVzs][nCentDca][nPtBins];
TH2D* h2Dca[nParticles][nPtBins];

TH1D* hTpcPiPlus[nCent];
TH1D* hTpcKMinus[nCent];

string outFileName = "D0.toyMc.Hijing.root";
std::pair<int, int> const decayChannels(747, 807);
std::pair<float, float> const momentumRange(0, 10);

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.0;
float const acceptanceEta= 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCent] = { 4.};

//============== main  program ==================
void toyMcEffZeroDecayLength(int npart = 200000)
{
   TStopwatch*   stopWatch = new TStopwatch();
   stopWatch->Start();
   gRandom->SetSeed();
   bookObjects();

   pydecay = TPythia6Decayer::Instance();
   pydecay->Init();

   setDecayChannels(763); // D0 --> Kpi

   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   for (int ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 10000))
         cout << "____________ ipart = " << ipart/static_cast<float>(npart) << " ________________" << endl;

      getKinematics(*b_d, M_D_0);
      // cout<<"M_D_0="<<M_D_0<<endl;

      decayAndFill(421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      // decayAndFill(-421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      if (ipart%1000 == 1) nt->AutoSave("SaveSelf");
   }

   write();
   stopWatch->Stop();   
   stopWatch->Print();
}

void setDecayChannels(int const mdme)
{
   for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
   TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
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

   fill(kf, b, weight, kMom, pMom, v00);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00)
{
   int const centrality = floor(nCent * gRandom->Rndm());
   // int const centrality = floor(0);

   TVector3 const vertex = getVertex(centrality);
   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

   v00 += vertex;

   // smear momentum
   TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution);
   TLorentzVector const pRMom = smearMom(pMom, fPionMomResolution);

   // smear position
   TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, v00);
   TVector3 const pRPos = smearPosData(0, vertex.z(), centrality, pRMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(pMom, pRMom, v00);

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
   TVector3 posDca1;
   TVector3 posDca2;
   float const dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, posDca1, posDca2, v0);
   float const decayLength = (v0 - vertex).Mag();
   float const dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
   float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
   float const angle12 = kRMom.Vect().Angle(pRMom.Vect());

   float const mposDca1 = (posDca1 - vertex).Mag();
   float const mposDca2 = (posDca2 - vertex).Mag();

   TLorentzVector kRMomRest = kRMom;
   TVector3 beta;
   beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
   kRMomRest.Boost(-beta);
   float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());

   int const charge = kf > 0? 1: -1;

   // save
   float arr[110];
   int iArr = 0;
   arr[iArr++] = centrality;
   arr[iArr++] = vertex.X();
   arr[iArr++] = vertex.Y();
   arr[iArr++] = vertex.Z();
   arr[iArr++] = getVzIndex(vertex.Z());

   arr[iArr++] = kf;
   arr[iArr++] = weight;
   arr[iArr++] = b->M();
   arr[iArr++] = b->Perp();
   arr[iArr++] = b->PseudoRapidity();
   arr[iArr++] = b->Rapidity();
   arr[iArr++] = b->Phi();
   arr[iArr++] = v00.X();
   arr[iArr++] = v00.Y();
   arr[iArr++] = v00.Z();

   arr[iArr++] = rMom.M();
   arr[iArr++] = rMom.Perp();
   arr[iArr++] = rMom.PseudoRapidity();
   arr[iArr++] = rMom.Rapidity();
   arr[iArr++] = rMom.Phi();
   arr[iArr++] = v0.X();
   arr[iArr++] = v0.Y();
   arr[iArr++] = v0.Z();
   arr[iArr++] = reconstructD0(centrality, rMom);

   arr[iArr++] = posDca1.X();
   arr[iArr++] = posDca1.Y();
   arr[iArr++] = posDca1.Z();
   arr[iArr++] = posDca2.X();
   arr[iArr++] = posDca2.Y();
   arr[iArr++] = posDca2.Z();
   arr[iArr++] = mposDca1;
   arr[iArr++] = mposDca2;

   arr[iArr++] = dca12;
   arr[iArr++] = decayLength;
   arr[iArr++] = dcaD0ToPv;
   arr[iArr++] = cosTheta;
   arr[iArr++] = angle12;
   arr[iArr++] = cosThetaStar;

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
   arr[iArr++] = kRPos.X();
   arr[iArr++] = kRPos.Y();
   arr[iArr++] = kRPos.Z();
   arr[iArr++] = kRDca;
   arr[iArr++] = kRSDca;
   arr[iArr++] = kRDcaXY;
   arr[iArr++] = kRDcaZ;
   arr[iArr++] = getEtaIndex(kRMom.PseudoRapidity());
   arr[iArr++] = getPtIndex(kRMom.Perp());
   arr[iArr++] = tpcReconstructed(1,-1*charge,centrality,kRMom);

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
   arr[iArr++] = pRPos.X();
   arr[iArr++] = pRPos.Y();
   arr[iArr++] = pRPos.Z();
   arr[iArr++] = pRDca;
   arr[iArr++] = pRSDca;
   arr[iArr++] = pRDcaXY;
   arr[iArr++] = pRDcaZ;
   arr[iArr++] = getEtaIndex(pRMom.PseudoRapidity());
   arr[iArr++] = getPtIndex(pRMom.Perp());
   arr[iArr++] = tpcReconstructed(0,charge,centrality,pRMom);

   arr[iArr++] = matchHft(1, vertex.z(), centrality, kRMom);
   arr[iArr++] = matchHft(0, vertex.z(), centrality, pRMom);

   nt->Fill(arr);
}

void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   // float const eta = gRandom->Uniform(-acceptanceEta, acceptanceEta);
   float const phi = TMath::TwoPi() * gRandom->Rndm();

   float const mT = sqrt(mass * mass + pt * pt);
   float const pz = mT * sinh(y);
   float const E = mT * cosh(y);

   b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
   // b.SetPtEtaPhiM(pt, eta, phi, mass); //change to eta cut
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
   TVector3 posDiff = pos - vertex;
   return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  float sign = posDiff.x()*p.y()-posDiff.y()*p.x() > 0 ? +1 : -1; 
 
  return sign*p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
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
  float sign = posDiff.x()*p.y()-posDiff.y()*p.x() > 0 ? +1 : -1; 
  return sign*newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  if(sin(p.Theta())==0) return 0;
  else return (-posDiff.x()*cos(p.Phi())-posDiff.y()*sin(p.Phi()))*cos(p.Theta())/sin(p.Theta())+posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& posDca1, TVector3& posDca2, TVector3& v0)
{
   TVector3 posDiff = pos2 - pos1;
   TVector3 pu1 = p1.Unit();
   TVector3 pu2 = p2.Unit();
   double pu1Pu2 = pu1.Dot(pu2);
   double g = posDiff.Dot(pu1);
   double k = posDiff.Dot(pu2);
   double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
   double s1 = g + s2 * pu1Pu2;
   // TVector3 posDca1 = pos1 + pu1 * s1;
   // TVector3 posDca2 = pos2 + pu2 * s2;
   posDca1 = pos1 + pu1 * s1;
   posDca2 = pos2 + pu2 * s2;
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

int getPtIndex(double pT)
{
   for (int i = 0; i < nPtBins; i++)
   {
      if ((pT >= ptEdge[i]) && (pT < ptEdge[i + 1]))
         return i;
   }
   return nPtBins - 1 ;
}

int getEtaIndex(double Eta)
{
   for (int i = 0; i < nEtas; i++)
   {
      if ((Eta >= EtaEdge[i]) && (Eta < EtaEdge[i + 1]))
         return i;
   }
   return nEtas - 1 ;
}

int getVzIndex(double Vz)
{
   for (int i = 0; i < nVzs; i++)
   {
      if ((Vz >= VzEdge[i]) && (Vz < VzEdge[i + 1]))
         return i;
   }
   return nVzs - 1 ;
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

// int getPhiIndexHftRatio(double Phi)
// {
//    for (int i = 0; i < nPhisHftRatio; i++)
//    {
//       if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
//          return i;
//    }
//    return nPhisHftRatio - 1 ;
// }

TVector3 smearPosData(int const iParticleIndex, double const vz, int cent, TLorentzVector const& rMom, TVector3 const& pos)
{
   int const iEtaIndex = getEtaIndex(rMom.PseudoRapidity());
   int const iVzIndex = getVzIndex(vz);
   int const iPtIndex = getPtIndex(rMom.Perp());

   double sigmaPosZ = 0;
   double sigmaPosXY = 0;

   if(cent == 8) cent = 7;

   gRandom->SetSeed();
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
     while ( fabs(rdmVz) > gVzCut);
   }

   return TVector3(0., 0., rdmVz);
}

bool reconstructD0(int const centrality, TLorentzVector const& mom)
{
   TGraph* gr = NULL;

   if (centrality < 4) gr = grEff[0];
   else if (centrality < 7) gr = grEff[1];
   else gr = grEff[2];

   return gRandom->Rndm() < gr->Eval(mom.Perp());
}

bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom)
{
  TH1D* h = NULL;

  if(iParticleIndex == 0)
  {
    h = hTpcPiPlus[cent];
  }
  else
  {
    h = hTpcKMinus[cent];
  }

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
   fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
   fKaonMomResolution = (TF1*)f.Get("fPion")->Clone("fKaon");
   f.Close();

   cout << "Loading input spectra ..." << endl;
   TFile fPP("pp200_spectra.root");
   fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("f1Levy");
   fPP.Close();

   TFile fVertex("hVz_D0Sample.root");

   for (int ii = 0; ii < nCent; ++ii)
   {
      h1Vz[ii]      = (TH1D*)(fVertex.Get(Form("mh1Vz")));
      h1Vz[ii]->SetDirectory(0);
   }

   fVertex.Close();

   cout << "Loading input HFT ratios and DCA ..." << endl;
   TFile fHftRatio1("Function_HFTMatch.root");
   fMatch = (TF1*)fHftRatio1.Get("fMatch")->Clone("fMatch");
   fHftRatio1.Close();

   TFile fDca1("2DProjection_DcaXyZ_sPHENIX.root");
   for (int iParticle = 0; iParticle < nParticles; ++iParticle)
   {
      // DCA
      for (int iPt = 0; iPt < nPtBins; ++iPt)
      {
         h2Dca[iParticle][iPt] = (TH2D*)(fDca1.Get(Form("mh2DcaPtPart_%i_%i", iParticle, iPt)));
         h2Dca[iParticle][iPt]->SetDirectory(0);
      }
      cout << "Finished loading Dca: " <<  endl;
   }

   fDca1.Close();

   cout << " Loading TPC tracking efficiencies " << endl;

   TFile fTpcPiPlus("Eff_PionPlus_embedding_FromHijing_1StopAll_16Oct21.root");
   TFile fTpcKMinus("Eff_KaonMinus_embedding_FromHijing_1StopAll_16Oct21.root");

   for(int iCent = 0; iCent< nCent; ++iCent)
   {
     hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(Form("h1Ratiocent_%i",iCent));
     hTpcPiPlus[iCent]->SetDirectory(0);
     hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(Form("h1Ratiocent_%i",iCent));
     hTpcKMinus[iCent]->SetDirectory(0);
   }

   fTpcPiPlus.Close();
   fTpcKMinus.Close();

   cout << "Done with loading all files ..." << endl;

   grEff[0] = new TGraph("eff_4080.csv", "%lg %lg", ",");
   grEff[1] = new TGraph("eff_1040.csv", "%lg %lg", ",");
   grEff[2] = new TGraph("eff_010.csv", "%lg %lg", ",");

   result = new TFile(outFileName.c_str(), "recreate");
   result->SetCompressionLevel(1);
   result->cd();

   int BufSize = (int)pow(2., 16.);
   // int Split = 1;
   nt = new TNtuple("nt", "", "cent:vx:vy:vz:vzIdx:"
                    "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                    "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:reco:" // Rc D0
                    "pos1x:pos1y:pos1z:pos2x:pos2y:pos2z:pos1:pos2:" // Rc pair
                    "dca12:decayLength:dcaD0ToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
                    "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                    "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kRSDca:kRDcaXY:kRDcaZ:kEtaIdx:kPtIdx:kTpc:" // Rc Kaon
                    "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                    "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:pRSDca:pRDcaXY:pRDcaZ:pEtaIdx:pPtIdx:pTpc:" // Rc Pion1
                    "kHft:pHft",BufSize);
   // nt->SetAutoSave(-500000); // autosave every 1 Mbytes
}
//___________
void write()
{
   result->cd();
   nt->Write();
   result->Close();
}
