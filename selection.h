#ifndef SELECTION_H
#define SELECTION_H
//-----------------------------------------------------------------------------
#include "histogramClass.h"
#include "functions.h"
#include "triggerFunctions.h"
#include "declarationsOfClasses.h"
#include "Math/LorentzVector.h"
#include <iostream>
#include <vector>
using namespace std;
//-----------------------------------------------------------------------------

Event::Event(TString histName, outputFile ofile_):
hist(histName, ofile_)
{ 

  countsTrackCriteria = new TH1D("countsTrackCriteria","countsTrackCriteria",1,0,1);
  countsTrackCriteria->SetBit(TH1::kCanRebin);
  countsTrackCriteria->SetStats(0);
  countsEventCuts = new TH1D("countsEventCuts","countsEventCuts",1,0,1);
  countsEventCuts->SetBit(TH1::kCanRebin);
  countsEventCuts->SetStats(0);
      
  onlyChi                = false;
  noChi                  = false;
  triggerRequirements    = false;
  trackPreselection      = false;
  qcdSupression          = false;
  trackCandidateCutFinal = false;

  TrackPtRequirement         = false;
  NumOfLostOuterRequirement  = false;
  CaloIsolationRequirement   = false;
  DeDxRequirement            = false;
  
  tightMuonCut                 = false;
  tightElectronCut             = false;
  TagAndProbeElectronCut       = false;
  TagAndProbeMuonCut           = false;
  TagAndProbeTauCut            = false;
  
  invertTrackPtRequirement         = false;
  invertCaloIsolationRequirement   = false;
  invertNumOfLostOuterRequirement  = false;
  invertDeDxRequirement            = false;
}

int Event::Selection()
{
    
  hist.variables.clearVectors();

  TrackColl.clear();
  JetColl.clear();
  MuonColl.clear();
  ElectronColl.clear();
      
  TrackColl=evt::Track;
  JetColl=evt::Jet;
  MuonColl=evt::MuonPFlow;
  ElectronColl=evt::ElectronPFlow;

  double dPhiMax=0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Blinding!!!!!
  if(edmEventHelper_isRealData){
   
    for(unsigned int i=0; i<TrackColl.size(); i++){
      //if(TrackColl[i].pt>=70 && TrackColl[i].ASmi>=0.4 /*&& trackCaloIsolation(&TrackColl[i])<10*/) return 0;
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  countsEventCuts->Fill("noCuts", weight);
  
  
  subleadingJetColl = getSubleadingJetCollection();
  // Special Collection
  if(onlyChi)
    {
      TrackColl = findChiInRecoTrackCollection(TrackColl,&hist);
    }
  else if(noChi)
    {
      TrackColl = getFakeTracksInTrackCollection(TrackColl);
    }
  
  //.................................................................................//
  //%%%%%%%%% Trigger Requirements %%%%%%%%%%%%%
  if(triggerRequirements){
    // 1.) Trigger Cut

    if(edmEventHelper_isRealData){
      int edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v           = getTriggerResult("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v");
      int edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v          = getTriggerResult("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v");
      int edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                           = getTriggerResult("edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v");
      int Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v  = getTriggerPrescales("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v");
      int Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v = getTriggerPrescales("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v");
      int Prescale_edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                  = getTriggerPrescales("edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v");
      hist.hMonoCentralPFJet80_PFMETnoMu95            -> Fill(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v,weight);
      hist.hMonoCentralPFJet80_PFMETnoMu105           -> Fill(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v,weight);
      hist.hMET120_HBHENoiseCleaned                   -> Fill(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v,weight);
      hist.hMonoCentralPFJet80_PFMETnoMu95_prescale   -> Fill(Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v,weight);
      hist.hMonoCentralPFJet80_PFMETnoMu105_prescale  -> Fill(Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v,weight);
      hist.hMET120_HBHENoiseCleaned_prescale          -> Fill(Prescale_edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v,weight);

      
      
      if(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v  == 1 ||
	 edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v == 1 ||
	 edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                  == 1)
	{}
      else{ 
	cout<<"something wrong with triggers!"<<endl;
	return 0; }
    }
    countsEventCuts->Fill("triggerCut_OnlyData", weight);
	     
    // 2.) MET cut
    if(MET_pt<=0.) return 0;
    countsEventCuts->Fill("metCut", weight);
    
    // 3.) Leading Jet Cut
    //if(!leadingJetRequirementsFullfilled(&JetColl[0], countsEventCuts)) return 0;
  }
  //.................................................................................//
  //%%%%%%%%% QCD supression %%%%%%%%%%%%%
  for(unsigned int i=0; i<subleadingJetColl.size(); i++){

    for(unsigned int j=i+1; j< subleadingJetColl.size(); j++){
	    
      double dPhi = std::abs(TVector2::Phi_mpi_pi( subleadingJetColl[i].phi- subleadingJetColl[j].phi));
      if(dPhi>dPhiMax) dPhiMax = dPhi;
      hist.hDeltaPhi->Fill(dPhi,weight);     
    }
  }
  hist.hDeltaPhiMaxbeforeCut->Fill(dPhiMax,weight);
  if(qcdSupression){
    //if(areTwoJetsBackToBack(subleadingJetColl)) return 0;
  }
  hist.hDeltaPhiMax->Fill(dPhiMax,weight);
  countsEventCuts->Fill("DeltaPhiCut", weight);

  // DeltaPhi(Jet,MET) cut
  if(qcdSupression)
    {
      //if(isMetInJetDirection(subleadingJetColl,MET_phi)) return 0;
    }
  countsEventCuts->Fill("1MetwithDeltaPhiMin2Jetsgt0p5", weight);

  //.................................................................................//
  //%%%%%%%%% Track Preselection %%%%%%%%%%%%%
  if(trackPreselection){
    // 4.)
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("NonEmptyTrkColl", weight);
    TrackColl = trackCandidateCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackCandCut", weight);
    
    // 5.)
    TrackColl = trackCleaningCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackCleaningCut", weight);

    // 6.)
    TrackColl = trackParticleMatchingCuts(TrackColl,subleadingJetColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackParticleMatchingCut", weight);
  }
  //.................................................................................//
  //%%%%%%%%% Tight Lepton Selection %%%%%%%%%%%%%
  if(tightMuonCut){

    MuonColl = getTightMuonsInEvent();
    if(MuonColl.size() != 2)                    return 0;
    if(MuonColl[0].charge == MuonColl[1].charge)                    return 0;
    TLorentzVector muon1 = lorentzVectorE(MuonColl[0].pt, MuonColl[0].eta, MuonColl[0].phi, MuonColl[0].energy);
    TLorentzVector muon2 = lorentzVectorE(MuonColl[1].pt, MuonColl[1].eta, MuonColl[1].phi, MuonColl[1].energy);
    double invMassMuons = (muon1+muon2).M();
    if(invMassMuons>100 || invMassMuons<80)  return 0;
    countsEventCuts->Fill("TwoTightMuonsInEvent", weight);
    matchMuonToGenParticle(MuonColl);
    //cout<<"invariant mass from muons= "<<invMassMuons<<endl;

  }
  
  if(tightElectronCut){

    ElectronColl = getTightElectronsInEvent();
    if(ElectronColl.size() != 2)                    return 0;
    if(ElectronColl[0].charge == ElectronColl[1].charge)                    return 0;
    TLorentzVector electron1 = lorentzVectorE(ElectronColl[0].pt, ElectronColl[0].eta, ElectronColl[0].phi, ElectronColl[0].energy);
    TLorentzVector electron2 = lorentzVectorE(ElectronColl[1].pt, ElectronColl[1].eta, ElectronColl[1].phi, ElectronColl[1].energy);
    double invMassElectrons = (electron1+electron2).M();
    if(invMassElectrons>100 || invMassElectrons<80) return 0;
    countsEventCuts->Fill("TwoTightElectronsInEvent", weight);
    //cout<<"invariant mass from electrons= "<<invMassElectrons<<endl;

  }

  if(TagAndProbeElectronCut){

    ElectronColl = getTightElectronsInEvent();
    if(ElectronColl.size() == 0)     return 0;

    std::vector<Track_s> TrackCollAux = TrackColl;
    TrackColl.clear();
    for(unsigned int j=0; j<ElectronColl.size(); j++){
      for(unsigned int i=0; i<TrackCollAux.size(); i++){

	if(ElectronColl[j].charge == TrackCollAux[i].charge)             continue;
	TLorentzVector electron = lorentzVectorM(ElectronColl[j].pt, ElectronColl[j].eta, ElectronColl[j].phi, 0.510998910 * pow(10,-3) );
	TLorentzVector track    = lorentzVectorM(TrackCollAux[i].pt, TrackCollAux[i].eta, TrackCollAux[i].phi, 0.510998910 * pow(10,-3) );
	double invMass          = (electron+track).M();
	if(invMass>100 || invMass<80)       continue;
	TrackColl.push_back(TrackCollAux[i]);
      }
    }

    if(TrackColl.size() == 0)    return 0;
    countsEventCuts->Fill("TagAndProbeConditionsFullFilled", weight);
  }

  if(TagAndProbeMuonCut){

    MuonColl = getTightMuonsInEvent();
    if(MuonColl.size() == 0)     return 0;

    std::vector<Track_s> TrackCollAux = TrackColl;
    TrackColl.clear();
    for(unsigned int j=0; j<MuonColl.size(); j++){
      for(unsigned int i=0; i<TrackCollAux.size(); i++){

	if(MuonColl[j].charge == TrackCollAux[i].charge)             continue;
	TLorentzVector muon     = lorentzVectorM(MuonColl[j].pt    , MuonColl[j].eta    , MuonColl[j].phi    ,  105.6583715 * pow(10,-3) );
	TLorentzVector track    = lorentzVectorM(TrackCollAux[i].pt, TrackCollAux[i].eta, TrackCollAux[i].phi,  105.6583715 * pow(10,-3) );
	double invMass          = (muon+track).M();
	if(invMass>100 || invMass<80)       continue;
	TrackColl.push_back(TrackCollAux[i]);
      }
    }

    if(TrackColl.size() == 0)    return 0;
    countsEventCuts->Fill("TagAndProbeConditionsFullFilled", weight);
  }

  if(TagAndProbeTauCut){

    MuonColl = getTightMuonsInEvent();
    if(MuonColl.size() == 0)     return 0;

    std::vector<Track_s> TrackCollAux = TrackColl;
    TrackColl.clear();
    for(unsigned int j=0; j<MuonColl.size(); j++){

      TLorentzVector muon     = lorentzVectorE(MuonColl[j].pt, MuonColl[j].eta, MuonColl[j].phi, MuonColl[j].energy);
      TLorentzVector met      = lorentzVectorE(MET_pt, MET_eta, MET_phi, MET_energy);
      double mt               = (muon+met).Mt();
      if(mt>40)               continue;

      for(unsigned int i=0; i<TrackCollAux.size(); i++){
	
	if(MuonColl[j].charge == TrackCollAux[i].charge)             continue;
	TLorentzVector track    = lorentzVectorM(TrackCollAux[i].pt, TrackCollAux[i].eta, TrackCollAux[i].phi, 1776.82 * pow(10,-3) );
	double invMass          = (muon+track).M();
	if(invMass>75 || invMass<40)       continue;
	TrackColl.push_back(TrackCollAux[i]);
      }
    }

    if(TrackColl.size() == 0)    return 0;
    countsEventCuts->Fill("TagAndProbeConditionsFullFilled", weight);
  }
  //.................................................................................//
  //%%%%%%%%% Final track cuts  BEGIN %%%%%%%%%%%%%
  // Final Track Cuts
  if(trackCandidateCutFinal)
    {
      TrackColl = finalTrackCuts(TrackColl,countsTrackCriteria);
      if(TrackColl.size()==0) return 0;
    }
  countsEventCuts->Fill("finalTrackCuts", weight);
  //%%%%%%%%% Final track cuts  END %%%%%%%%%%%%%%%
  //.................................................................................//
  //matchTrackToSimTrack(TrackColl);     
  matchTrackToGenParticle(TrackColl);
  
  hist.FillTrackVariables(TrackColl,weight);
  hist.FillCharginoHistograms(ChiTrack,weight);
  hist.hMet->Fill(evt::MET_pt,weight);
  hist.hLuminosityBlock->Fill(evt::edmEventHelper_luminosityBlock, weight);
  if(JetColl.size()!=0){
    hist.h1stjetpt->Fill(JetColl[0].pt,weight);
    hist.variables.LeadingJetPt = JetColl[0].pt;
  }
  hist.variables.met       = MET_pt;
  hist.variables.nJets     = subleadingJetColl.size();
  hist.variables.VertexZ   = evt::Vertex[0].z;
  hist.hnPFJetsub->Fill(subleadingJetColl.size(),weight);

  for(unsigned int i=0; i<TrackColl.size(); i++){
    hist.htrackpt1stjetpt->Fill(TrackColl[i].pt,JetColl[0].pt,weight);  
  }

  // Fill muon variables
  for(int i=0; i<MuonColl.size(); i++){

    hist.variables.muonPt.push_back(MuonColl[i].pt);
    hist.variables.muonEt.push_back(MuonColl[i].et);
    hist.variables.muonE.push_back(MuonColl[i].energy);
    hist.variables.muonPhi.push_back(MuonColl[i].phi);
    hist.variables.muonEta.push_back(MuonColl[i].eta);
    hist.variables.muonChi2.push_back(MuonColl[i].globalTrack_chi2);
    hist.variables.muonndof.push_back(MuonColl[i].globalTrack_ndof);
    hist.variables.muondB.push_back(MuonColl[i].dB);
    hist.variables.muond0.push_back(MuonColl[i].globalTrack_d0);
    hist.variables.muonVertex_z.push_back(MuonColl[i].vertex_z);
    hist.variables.muonNOfValidMuonHits.push_back(MuonColl[i].globalTrack_hitPattern_numberOfValidMuonHits);
    hist.variables.muonTrackerLayersWithMeas.push_back(MuonColl[i].innerTrack_hitPattern_trackerLayersWithMeasurement);
    hist.variables.muonNOfValidPixelHits.push_back(MuonColl[i].innerTrack_hitPattern_numberOfValidPixelHits);
    hist.variables.muonNOfMatchedMuonStations.push_back(MuonColl[i].numberOfMatchedStations);
    hist.variables.muonChargedHadronIso.push_back(MuonColl[i].chargedHadronIso);
    hist.variables.muonNeutralHadronIso.push_back(MuonColl[i].neutralHadronIso);
    hist.variables.muonPhotonIso.push_back(MuonColl[i].photonIso);
    hist.variables.muonpuChargedHadronIso.push_back(MuonColl[i].puChargedHadronIso);
    hist.variables.muonPdgId.push_back(MuonColl[i].pdgId);
    hist.variables.muonStatus.push_back(MuonColl[i].status);
  }


  //-----------------------------------------------

  if(chipmGenParticle.size()>0) hist.FillGenParticleHistograms(chipmGenParticle,weight);


  hist.tree->Fill();

  return 0;

};


std::vector<Track_s> Event::finalTrackCuts(std::vector<Track_s> trackCollection, TH1D* countsTrackCriteria){

  std::vector<Track_s> outputCollection;
  if(trackCollection.size()==0) return trackCollection;

  TrackColl.clear();

  for(unsigned int i=0; i<1; i++){

    //.................................................................................//
    if(TrackPtRequirement){
      if(invertTrackPtRequirement){
	if(trackCollection[i].pt>35.)                                             continue;
      }
      else{
	if(trackCollection[i].pt<=35.)                                            continue;
      }
    }
    countsTrackCriteria->Fill("PtGreater70GeV", weight);
    //.................................................................................//
    if(CaloIsolationRequirement){
      if(invertCaloIsolationRequirement){
	if(trackCaloIsolation(&trackCollection[i])<=5)                           continue;
      }
      else{
	if(trackCaloIsolation(&trackCollection[i])>5)                            continue;
      }
    }
    countsTrackCriteria->Fill("CaloIsolation0p5", weight);
    //.................................................................................//
    if(NumOfLostOuterRequirement){
      if(invertNumOfLostOuterRequirement){
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits>=1)           continue;
      }
      else{
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits<1)            continue;
      }
    }
    countsTrackCriteria->Fill("NOfLostHitsOuterGe1", weight);
    //.................................................................................//
    /*
    if(DeDxRequirement){
      if(invertDeDxRequirement){
	if(trackCollection[i].ASmi>=0.2){
	  continue;
	}
      }
      else{
	if(trackCollection[i].ASmi<0.2)                                                              continue;
      }
    }
    countsTrackCriteria->Fill("DeDxASmiGe0p2", weight);
    */
    //.................................................................................//

    outputCollection.push_back(trackCollection[i]);
  }


  return outputCollection;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//--------------------------------------------------------------------------------------------------
#endif

