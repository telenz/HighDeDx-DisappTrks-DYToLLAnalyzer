#ifndef FUNCTIONS_H
#define FUNCTIONS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "chiClass.h"
#include "declarationsOfClasses.h"
#include "TH3.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------

struct evt::GenParticle_s chipGenParticle;
struct evt::GenParticle_s chimGenParticle;

std::vector<evt::GenParticle_s> chipmGenParticle;


std::vector<double> etaCSC, phiCSC, etaEcal, phiEcal;  

bool zeroChip = true;
bool zeroChim = true;
bool zeroJet  = true;

int nChi0 = 0;


int mismatchedGenChiToTrack = 0;
int matchedGenChiToTrack    = 0;

int nChiInSimTrack          = 0;
int nChiInSimVertex         = 0;

//--------------------------------------------------------------------------------------------------

double trackCaloIsolation(struct Track_s *track){


  double caloTowerIso05RhoCorr = 
    std::max(0.,track->caloHadDeltaRp5 + track->caloEMDeltaRp5 - sdouble_value*TMath::Pi()*0.5*0.5); 

  return caloTowerIso05RhoCorr;
}

//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedTau(struct Track_s* track){


  for(unsigned int i=0; i<Tau.size(); i++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0; 

    if(Tau[i].pt<20)                                    continue;
    if(abs(Tau[i].eta)>2.3)                             continue;
    if(Tau[i].byLooseCombinedIsolationDeltaBetaCorr<=0) continue;
    if(Tau[i].decayModeFinding<=0)                      continue;
    if(Tau[i].againstElectronLoose<=0)                  continue;
    if(Tau[i].againstMuonTight<=0)                      continue;

    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Tau[i].phi));
    dEta = std::abs(track->eta - Tau[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedElectron(struct Track_s* track){


  for(unsigned int i=0; i<Electron.size(); i++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0; 

    if(Electron[i].mvaNonTrigV0<0) continue;
    if(Electron[i].pt<10)          continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Electron[i].phi));
    dEta = std::abs(track->eta - Electron[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedMuon(struct Track_s* track){

  for(unsigned int i=0; i<Muon.size(); i++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0;

    if(Muon[i].pt<10.)   continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Muon[i].phi));
    dEta = std::abs(track->eta - Muon[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedJet(struct evt::Track_s track, std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(track.phi-jetColl[i].phi));
    double dEta = std::abs(track.eta - jetColl[i].eta);
    double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2)); 

    if(dR<0.5) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
/*****************************************
  Find chargino in GenParticle collection 
*****************************************/
void findChiInGenParticleCollection(){

  zeroChip = true;
  zeroChim = true;

  chipmGenParticle.clear();
  ChiTrack.clear();
  ChiTrack.resize(2);

  int idx=0;

  for(int i=0; i<evt::nGenParticle; i++){

    if(abs(evt::GenParticle[i].pdgId)==1000024){

      if(evt::GenParticle[i].pdgId>0 && zeroChip){

	chipGenParticle = evt::GenParticle[i];
	chipmGenParticle.push_back(evt::GenParticle[i]);
	fillChiTrackWithGenParticleVariables(&ChiTrack[idx], &evt::GenParticle[i]);
	idx+=1;
	zeroChip = false;
      }
      else if(evt::GenParticle[i].pdgId<0 && zeroChim){

	chimGenParticle = evt::GenParticle[i];
	chipmGenParticle.push_back(evt::GenParticle[i]);
	fillChiTrackWithGenParticleVariables(&ChiTrack[idx], &evt::GenParticle[i]);
	idx+=1;
	zeroChim = false;
      }	      
    }
    if(!zeroChip && !zeroChim) break;
  }

  //if(zeroChip || zeroChim) cout<<"To few charginos in GenParticle collection!"<<endl;
}
//--------------------------------------------------------------------------------------------------
// Get only the Chi from the full Track Collection
std::vector<Track_s> findChiInRecoTrackCollection(std::vector<Track_s>& trkCollection, Hist* hist){

  std::vector<Track_s> chiTrackCollection;
  chiTrackCollection.clear();

  int idxMin  = -100;
 
  for(unsigned int j=0; j<ChiTrack.size();j++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0;
    double dRmin=10000.;

    for(unsigned int i=0; i<trkCollection.size(); i++){
      dPhi = std::abs(TVector2::Phi_mpi_pi(trkCollection[i].phi - ChiTrack[j].genphi));
      dEta = std::abs(trkCollection[i].eta - ChiTrack[j].geneta);
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );

      if(dR<dRmin){
	dRmin=dR;
	idxMin=i;
      }
    }
    
    if(dRmin<0.01){
      chiTrackCollection.push_back(trkCollection[idxMin]);
      fillChiTrackWithRecoTrackVariables(&ChiTrack[j], &trkCollection[idxMin]);
      ChiTrack[j].matched = true;
      matchedGenChiToTrack += 1;
    }
    hist->htrackDeltaRSimRecoTracks -> Fill(dRmin,weight);
    hist->hSimTrackType             -> Fill(std::abs(ChiTrack[j].SimTracktype),weight);

  }

  return chiTrackCollection;
}
//--------------------------------------------------------------------------------------------------
struct evt::GenParticle_s  findLeadingJetInGenParticleCollectionWithPtGt30(){

  struct evt::GenParticle_s leadingJetGenParticle;
  zeroJet = true;
  for(int i=0; i<evt::nGenParticle; i++){

    if(abs(evt::GenParticle[i].pdgId)<=6){
      if(evt::GenParticle[i].pt<30) continue;
      leadingJetGenParticle = evt::GenParticle[i];
      zeroJet = false;
      break;
    }
	       
  }
  return leadingJetGenParticle;
}

//--------------------------------------------------------------------------------------------------
bool areTwoJetsBackToBack(std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    for(unsigned int j=i+1; j<jetColl.size(); j++){

      double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-jetColl[j].phi));
      if(abs(dPhi)>=2.5) return true;     
      
    }
  }

  return false;

}
//--------------------------------------------------------------------------------------------------
bool isMetInJetDirection(std::vector<evt::Jet_s> jetColl, double metPhi){

  for(unsigned int i=0; i<jetColl.size(); i++){

    if(i>1) break;
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-metPhi));
    if(dPhi<=0.5) return true;     
    
  }

  return false;

}


//--------------------------------------------------------------------------------------------------
bool leadingJetRequirementsFullfilled(struct evt::Jet_s* leadingJet, TH1D* countsEventCuts){

  // jetId: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  if(leadingJet==0)                                return false;
  if(leadingJet->pt<=110.)                         return false;
  countsEventCuts->Fill("leadingJetPtGt110GeV", evt::weight);
  if(std::abs(leadingJet->eta)>=2.4)               return false;
  countsEventCuts->Fill("absLeadJetEtaLt2p4", evt::weight);
  if(leadingJet->chargedHadronEnergyFraction<=0.2) return false;
  countsEventCuts->Fill("CHEFgt0p2", evt::weight);
  if(leadingJet->chargedEmEnergyFraction>=0.5)     return false;
  countsEventCuts->Fill("CHEmEFle0p5", evt::weight);
  if(leadingJet->neutralHadronEnergyFraction>=0.7) return false;
  countsEventCuts->Fill("NHEFle0p7", evt::weight);
  if(leadingJet->neutralEmEnergyFraction>=0.7)     return false;
  countsEventCuts->Fill("NEmEFle0p7", evt::weight);

  return true;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Jet_s>  getSubleadingJetCollection(){

  std::vector<evt::Jet_s> jetCollection;
  jetCollection.clear();
  for(unsigned int i=0; i<evt::Jet.size(); i++){

    //bool isCharginoCandidate = false;

    if(evt::Jet[i].pt<=20.)                          continue;
    if(std::abs(evt::Jet[i].eta)>=4.5)               continue;
    //if(evt::Jet[i].neutralHadronEnergyFraction>=0.7) continue;
    //if(evt::Jet[i].chargedEmEnergyFraction>=0.5)     continue;

    jetCollection.push_back(evt::Jet[i]);
  }

  return jetCollection;
}
//--------------------------------------------------------------------------------------------------
bool isGoodVertex(){
    
  if(evt::Vertex[0].z > 24.)                                                     return false;
  if(std::sqrt(std::pow(evt::Vertex[0].x,2) + std::pow(evt::Vertex[0].y,2)) > 2) return false;
  if(evt::Vertex[0].ndof < 4)                                                    return false;
  return true;
}
//--------------------------------------------------------------------------------------------------
bool getTrkIsMatchedDeadEcal(struct evt::Track_s *track){

  for(unsigned int i=0; i<etaEcal.size(); i++){

    double dPhi  = 0;
    double dEta  = 0;
    double dR    = 0;  
        
    dPhi = std::abs(TVector2::Phi_mpi_pi(phiEcal[i] - track->phi));
    dEta = std::abs(etaEcal[i] - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.05) return true;
  }  
  return false;
}
//--------------------------------------------------------------------------------------------------
bool getTrkIsMatchedBadCSC(struct evt::Track_s *track){

  for(unsigned int i=0; i<etaCSC.size(); i++){

    double dPhi  = 0;
    double dEta  = 0;
    double dR    = 0;  
          
    dPhi = std::abs(TVector2::Phi_mpi_pi(phiCSC[i] - track->phi));
    dEta = std::abs(etaCSC[i] - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.25) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isWithinIntermoduleGapsOfECAL(struct evt::Track_s *track){

  if(track->eta<-1.14018   && track->eta>-1.1439)     return true;
  if(track->eta<-0.791884  && track->eta>-0.796051)   return true;
  if(track->eta<-0.44356   && track->eta>-0.447911)   return true;
  if(track->eta<0.00238527 && track->eta>-0.00330793) return true;
  if(track->eta<0.446183   && track->eta>0.441949)    return true;
  if(track->eta<0.793955   && track->eta>0.789963)    return true;
  if(track->eta<1.14164    && track->eta>1.13812)     return true;
  
  return false;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCuts(std::vector<evt::Track_s> inputColl, bool keepCriteria)
{

  std::vector<evt::Track_s> outputColl;

  for(unsigned int i=0; i<inputColl.size(); i++){
    if(!keepCriteria) continue;
    outputColl.push_back(inputColl[i]);
  }

  return outputColl;

}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCandidateCuts(std::vector<evt::Track_s> trackCollection, TH1D* countsTrackCriteria)
{

  std::vector<evt::Track_s> outputColl;
  bool firstTrack1 = true;
  bool firstTrack2 = true;
  bool firstTrack3 = true;
  bool firstTrack4 = true;

  for(unsigned int i=0; i<trackCollection.size(); i++){
    if(firstTrack1){
      countsTrackCriteria->Fill("beforeTrackCriteria", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    if(trackCollection[i].pt<=10.)                                            continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("PtGreater10GeV", weight);  
      firstTrack2 = false;
    }
    //.................................................................................//
    if(std::abs(trackCollection[i].eta)>2.4)                                  continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("EtaLess2p4", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(!trackCollection[i].trackHighPurity)                                   continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("highPurity", weight);
      firstTrack4 = false;
    }

    //.................................................................................//
    outputColl.push_back(trackCollection[i]);
  }
  return outputColl;
  
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCleaningCuts(std::vector<evt::Track_s> trackCollection, TH1D* countsTrackCriteria)
{

  std::vector<evt::Track_s> outputColl;

  bool firstTrack1  = true;
  bool firstTrack2  = true;
  bool firstTrack3  = true;
  bool firstTrack4  = true;
  bool firstTrack5  = true;
  bool firstTrack6  = true;
  bool firstTrack7  = true;
  bool firstTrack8  = true;
  bool firstTrack9  = true;
  bool firstTrack10 = true;
  bool firstTrack11 = true;
  bool firstTrack12 = true;


  for(unsigned int i=0; i<trackCollection.size(); i++){

    //.................................................................................//
    //if(std::abs(trackCollection[i].eta)>1.42 && std::abs(trackCollection[i].eta)<1.65)        continue;
    if(firstTrack1){
      countsTrackCriteria->Fill("EtaLess1p42Gt1p65", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    //if(std::abs(trackCollection[i].eta)>0.15 && std::abs(trackCollection[i].eta)<0.35)        continue;
    if(firstTrack1){
      countsTrackCriteria->Fill("EtaLess0p15Gt0p35", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    //if(std::abs(trackCollection[i].eta)>1.55 && std::abs(trackCollection[i].eta)<1.85)        continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("EtaLess1p55Gt1p85", weight);
      firstTrack2 = false;
    }
    //.................................................................................//
    if(getTrkIsMatchedDeadEcal(&trackCollection[i]))                                          continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("isMatchedDeadEcal", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(isWithinIntermoduleGapsOfECAL(&trackCollection[i]))                                    continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("notWithinECALGap", weight);
      firstTrack4 = false;
    }
    //.................................................................................//
    if(getTrkIsMatchedBadCSC(&trackCollection[i]))                                            continue;
    if(firstTrack5){
      countsTrackCriteria->Fill("isMatchedBadCSC", weight);
      firstTrack5 = false;
    }
    //.................................................................................//
    double _dvx = trackCollection[i].vx - Vertex[0].x;
    double _dvy = trackCollection[i].vy - Vertex[0].y;
    double d0 = abs( - _dvx*trackCollection[i].py + _dvy*trackCollection[i].px)/trackCollection[i].pt;
    if(abs(d0)>0.02)                                                                          continue;
    if(firstTrack6){
      countsTrackCriteria->Fill("d0Less0p2mm", weight);
      firstTrack6 = false;
    }
    //.................................................................................//
    double _dvz = trackCollection[i].vz - Vertex[0].z;
    double dZ = _dvz - ( _dvx*trackCollection[i].px + _dvy*trackCollection[i].py)/trackCollection[i].pt * (trackCollection[i].pz/trackCollection[i].pt);
    if(abs(dZ)>0.5)                                                                           continue;
    if(firstTrack7){
      countsTrackCriteria->Fill("dZLess5mm", weight);
      firstTrack7 = false;
    }
    //.................................................................................//
    //if(trackCollection[i].numberOfValidHits<7)                                                continue;
    if(firstTrack8){
      countsTrackCriteria->Fill("NOfValidHitsGreater7", weight);
      firstTrack8 = false;
    }
    //.................................................................................//
    if(trackCollection[i].hitPattern_trackerLayersWithoutMeasurement>0)                       continue;
    if(firstTrack9){
      countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", weight);
      firstTrack9 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackerExpectedHitsInner_numberOfLostHits>0)                        continue;
    if(firstTrack10){
      countsTrackCriteria->Fill("NOfLostHitsInnerEq0", weight);
      firstTrack10 = false;
    }
    //.................................................................................//
    //if(trackCollection[i].hitPattern_trackerLayersWithoutMeasurement==0 && trackCollection[i].trackerExpectedHitsInner_numberOfLostHits==0)            continue;
    if(firstTrack11){
      countsTrackCriteria->Fill("missingMiddleAndInnerHits", weight);
      firstTrack11 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackRelIso03>=0.1)                                                 continue;
    if(firstTrack12){
      countsTrackCriteria->Fill("TrackIsolationDeltaR0p3Less0p1", weight);
      firstTrack12 = false;
    }
    //.................................................................................//

    outputColl.push_back(trackCollection[i]);
    //.................................................................................//
  
  }
  
  
  return outputColl;
  
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackParticleMatchingCuts(std::vector<evt::Track_s> trackCollection, std::vector<Jet_s>& jetColl, TH1D* countsTrackCriteria)
{

  std::vector<evt::Track_s> outputColl;

  bool firstTrack1  = true;
  bool firstTrack2  = true;
  bool firstTrack3  = true;
  bool firstTrack4  = true;

  for(unsigned int i=0; i<trackCollection.size(); i++){

    //.................................................................................//
    if(isTrackReconstructedJet(trackCollection[i], jetColl))                                  continue;
    if(firstTrack1){
      countsTrackCriteria->Fill("InJetCollectionR0p5", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedTau(&trackCollection[i]))                                          continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("InTauCollectionR0p15", weight);
      firstTrack2 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedElectron(&trackCollection[i]))                                     continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("InElectronCollectionR0p15", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedMuon(&trackCollection[i]))                                         continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("InMuonCollectionR0p15", weight);
      firstTrack4 = false;
    }
    //.................................................................................//
    outputColl.push_back(trackCollection[i]);
    //.................................................................................//

  }
  
  return outputColl;
  
}
//--------------------------------------------------------------------------------------------------
std::vector<Track_s> getFakeTracksInTrackCollection(const std::vector<Track_s>& inputCollection){

  std::vector<Track_s> outputCollection;
  outputCollection.clear();
  double dRchip = 0;
  double dRchim = 0;

  for(unsigned int i=0; i<inputCollection.size(); i++){

    dRchip=10000;
    dRchim=10000;

    if(!zeroChip){
      double dPhichip = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chipGenParticle.phi));
      double dEtachip = std::abs(inputCollection[i].eta - chipGenParticle.eta);
      dRchip   = std::sqrt(pow(dPhichip,2) + pow(dEtachip,2));
    }
    if(!zeroChim){
      double dPhichim = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chimGenParticle.phi));
      double dEtachim = std::abs(inputCollection[i].eta - chimGenParticle.eta);
      dRchim   = std::sqrt(pow(dPhichim,2) + pow(dEtachim,2));
    }
    if(dRchim>0.01 && dRchip>0.01) outputCollection.push_back(inputCollection[i]);
  }
  return outputCollection;
  
}
//--------------------------------------------------------------------------------------------------
// Match Tracks to a generator particle in the GenCollection
void matchTrackToGenParticle(std::vector<Track_s>& inputCollection){

  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;
  


  for(unsigned int i=0; i<inputCollection.size(); i++){

    inputCollection[i].pdgId=0;
    inputCollection[i].status=-1;
    inputCollection[i].beta=10;
    double dRmin = 0.01;
    int    idx   = -1;

    
    for(unsigned int j=0; j<GenParticle.size(); j++){

      if(GenParticle[j].status==2) continue;

      dEta = std::abs(inputCollection[i].eta - GenParticle[j].eta);
      dPhi = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - GenParticle[j].phi));
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );
      if(dR<dRmin){
	dRmin = dR;
	idx = j;
      }
    }
    if(dRmin<0.01){
	inputCollection[i].pdgId = GenParticle[idx].pdgId;
	inputCollection[i].status= GenParticle[idx].status;
	inputCollection[i].beta  = GenParticle[idx].p/GenParticle[idx].energy;
	inputCollection[i].genPt = GenParticle[idx].pt;
	inputCollection[i].genE  = GenParticle[idx].energy;
	inputCollection[i].genEt = GenParticle[idx].et;
    }
  }
}
//--------------------------------------------------------------------------------------------------
double getLifetime(string filename){

  unsigned found         = filename.find("ctau");
  string TargetLifetime  = filename.substr(found + 5);
  found                  = TargetLifetime.find("cm");
  TargetLifetime         = TargetLifetime.substr(0,found);
  //cout<<"Lifetime  = "<<TargetLifetime<<endl;

  return atoi(TargetLifetime.c_str());

}
//--------------------------------------------------------------------------------------------------
bool isMuonTight(struct evt::MuonPFlow_s* muon){
  
  if( !muon->isGlobalMuon )                                          return false;
  if( !muon->isPFMuon )                                              return false;
  if( muon->globalTrack_chi2/muon->globalTrack_ndof>10 )             return false;
  if( muon->globalTrack_hitPattern_numberOfValidMuonHits<1 )         return false;
  if( muon->numberOfMatchedStations<2 )                              return false;
  if( muon->innerTrack_hitPattern_trackerLayersWithMeasurement<6 )   return false;
  if( muon->innerTrack_hitPattern_numberOfValidPixelHits<1 )         return false;
  if( muon->dB>=0.2 )                                                return false;
  if( std::abs( muon->vertex_z - evt::Vertex[0].z ) >= 0.5 )         return false;

  return true;
}
//--------------------------------------------------------------------------------------------------
float relPFIso(struct evt::MuonPFlow_s* muon){

  float relIso = 
    ( muon->chargedHadronIso + 
      std::max( 0.0, muon->neutralHadronIso + muon->photonIso - 0.5*muon->puChargedHadronIso ) ) / muon->pt; 

  return relIso;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::MuonPFlow_s>  getTightMuonsInEvent(){

  std::vector<evt::MuonPFlow_s> muonCollection;
  muonCollection.clear();
  
  for(unsigned int i=0; i<evt::MuonPFlow.size(); i++){

    if( evt::MuonPFlow[i].pt<=25 )                     continue;
    if( std::abs(evt::MuonPFlow[i].eta)>2.5 )          continue;
    if( !isMuonTight(&evt::MuonPFlow[i]) )             continue;
    if( relPFIso(&MuonPFlow[i])>0.12 )                   continue;

    muonCollection.push_back(evt::MuonPFlow[i]);
  }

  return muonCollection;
}
//--------------------------------------------------------------------------------------------------
float isGoodMVAElectron(struct evt::ElectronPFlow_s* electron){

  double sceta = fabs(electron->superCluster_eta);

  if(electron->pt >= 20.0) {

    if((0.0 <= sceta) && (sceta <= 0.8)) {
      if(electron->mvaTrigV0 > 0.94) {
	return true;
      }
    } else if((0.8 < sceta) && (sceta <= 1.479)) {
      if(electron->mvaTrigV0 > 0.85) {
	return true;
      }
    } else if((1.479 < sceta) && (sceta <= 2.5) ) {
      if(electron->mvaTrigV0 > 0.92) {
	return true;
      }
    }
	  
  }else if(electron->pt > 10.0) {
    if((0.0 <= sceta) && (sceta <= 0.8)) {
      if(electron->mvaTrigV0 > 0.00) {
	return true;
      }
    } else if((0.8 < sceta) && (sceta <= 1.479)) {
      if(electron->mvaTrigV0 > 0.10) {
	return true;
      }
    } else if((1.479 < sceta) && (sceta <= 2.5) ) {
      if(electron->mvaTrigV0 > 0.62) {
	return true;
      }
    } 
  }
	

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isElectronTight(struct evt::ElectronPFlow_s* electron){
  
  if( !electron->passConversionVeto )                                          return false;
  if( electron->gsfTrack_trackerExpectedHitsInner_numberOfLostHits >0 )        return false;
  if( !isGoodMVAElectron(electron) )                                           return false;

  return true;
}
//--------------------------------------------------------------------------------------------------
float relPFIsoRho(struct evt::ElectronPFlow_s* electron){

  float relIsoRho =
    ( electron->chargedHadronIso + 
      std::max( 0.0, electron->neutralHadronIso + electron->photonIso - sdoublePF_value*electron->Aeff04 ) ) 
    / electron->pt; 

  return relIsoRho;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::ElectronPFlow_s>  getTightElectronsInEvent(){

  std::vector<evt::ElectronPFlow_s> electronCollection;
  electronCollection.clear();
  
  for(unsigned int i=0; i<evt::ElectronPFlow.size(); i++){

    if( evt::ElectronPFlow[i].pt<=25 )                     continue;
    if( std::abs(evt::ElectronPFlow[i].eta)>2.5 )          continue;
    if( !isElectronTight(&evt::ElectronPFlow[i]) )         continue;
    if( relPFIsoRho(&ElectronPFlow[i])>0.15 )              continue;

    electronCollection.push_back(evt::ElectronPFlow[i]);
  }

  return electronCollection;
}

//--------------------------------------------------------------------------------------------------
TLorentzVector lorentzVectorE(float pt, float eta, float phi, float energy){
  
  TLorentzVector v4;
  v4.SetPtEtaPhiE(pt,eta,phi,energy);
  
  return v4;
};

TLorentzVector lorentzVectorM(float pt, float eta, float phi, float mass){
  
  TLorentzVector v4;
  v4.SetPtEtaPhiM(pt,eta,phi,mass);
  
  return v4;
};
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Match Tracks to a generator particle in the GenCollection
void matchMuonToGenParticle(std::vector<MuonPFlow_s>& inputCollection){

  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;
  


  for(unsigned int i=0; i<inputCollection.size(); i++){

    inputCollection[i].pdgId=0;
    inputCollection[i].status=-1;
    double dRmin = 0.01;
    int    idx   = -1;

    
    for(unsigned int j=0; j<GenParticle.size(); j++){

      dEta = std::abs(inputCollection[i].eta - GenParticle[j].eta);
      dPhi = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - GenParticle[j].phi));
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );
      if(dR<dRmin){
	dRmin = dR;
	idx = j;
      }
    }
    if(dRmin<0.01){
	inputCollection[i].pdgId = GenParticle[idx].pdgId;
	inputCollection[i].status= GenParticle[idx].status;
    }
  }
}
#endif

