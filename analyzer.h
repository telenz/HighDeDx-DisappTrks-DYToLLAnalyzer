#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Fri Jun 12 11:39:02 2015 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
// -- System

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>

#include "analyzerutil.h"
#include "treestream.h"
#include "pdg.h"

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"

namespace evt {
//-----------------------------------------------------------------------------
// --- Declare variables
//-----------------------------------------------------------------------------
std::vector<float>	ElectronPFlow_Aeff04(200,0);
std::vector<int>	ElectronPFlow_charge(200,0);
std::vector<float>	ElectronPFlow_chargedHadronIso(200,0);
std::vector<float>	ElectronPFlow_energy(200,0);
std::vector<float>	ElectronPFlow_et(200,0);
std::vector<float>	ElectronPFlow_eta(200,0);
std::vector<int>	ElectronPFlow_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(200,0);
std::vector<float>	ElectronPFlow_mvaNonTrigV0(200,0);
std::vector<float>	ElectronPFlow_mvaTrigV0(200,0);
std::vector<float>	ElectronPFlow_neutralHadronIso(200,0);
std::vector<int>	ElectronPFlow_passConversionVeto(200,0);
std::vector<float>	ElectronPFlow_phi(200,0);
std::vector<float>	ElectronPFlow_photonIso(200,0);
std::vector<float>	ElectronPFlow_pt(200,0);
std::vector<float>	ElectronPFlow_puChargedHadronIso(200,0);
std::vector<float>	ElectronPFlow_pz(200,0);
std::vector<float>	ElectronPFlow_superCluster_eta(200,0);
std::vector<float>	Electron_Aeff04(200,0);
std::vector<int>	Electron_charge(200,0);
std::vector<float>	Electron_chargedHadronIso(200,0);
std::vector<float>	Electron_energy(200,0);
std::vector<float>	Electron_et(200,0);
std::vector<float>	Electron_eta(200,0);
std::vector<int>	Electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(200,0);
std::vector<float>	Electron_mvaNonTrigV0(200,0);
std::vector<float>	Electron_mvaTrigV0(200,0);
std::vector<float>	Electron_neutralHadronIso(200,0);
std::vector<int>	Electron_passConversionVeto(200,0);
std::vector<float>	Electron_phi(200,0);
std::vector<float>	Electron_photonIso(200,0);
std::vector<float>	Electron_pt(200,0);
std::vector<float>	Electron_puChargedHadronIso(200,0);
std::vector<float>	Electron_pz(200,0);
std::vector<float>	Electron_superCluster_eta(200,0);
std::vector<int>	GenParticle_charge(2000,0);
std::vector<float>	GenParticle_energy(2000,0);
std::vector<float>	GenParticle_et(2000,0);
std::vector<float>	GenParticle_eta(2000,0);
std::vector<float>	GenParticle_mass(2000,0);
std::vector<float>	GenParticle_p(2000,0);
std::vector<int>	GenParticle_pdgId(2000,0);
std::vector<float>	GenParticle_phi(2000,0);
std::vector<float>	GenParticle_pt(2000,0);
std::vector<float>	GenParticle_pz(2000,0);
std::vector<int>	GenParticle_status(2000,0);
std::vector<float>	Jet_chargedEmEnergyFraction(200,0);
std::vector<float>	Jet_chargedHadronEnergyFraction(200,0);
std::vector<float>	Jet_energy(200,0);
std::vector<float>	Jet_et(200,0);
std::vector<float>	Jet_eta(200,0);
std::vector<float>	Jet_neutralEmEnergyFraction(200,0);
std::vector<float>	Jet_neutralHadronEnergyFraction(200,0);
std::vector<float>	Jet_phi(200,0);
std::vector<float>	Jet_pt(200,0);
std::vector<float>	Jet_pz(200,0);
float	MET_energy;
float	MET_et;
float	MET_eta;
float	MET_phi;
float	MET_pt;
float	MET_pz;
std::vector<int>	MuonPFlow_charge(200,0);
std::vector<float>	MuonPFlow_chargedHadronIso(200,0);
std::vector<float>	MuonPFlow_dB(200,0);
std::vector<float>	MuonPFlow_energy(200,0);
std::vector<float>	MuonPFlow_et(200,0);
std::vector<float>	MuonPFlow_eta(200,0);
std::vector<float>	MuonPFlow_globalTrack_chi2(200,0);
std::vector<float>	MuonPFlow_globalTrack_d0(200,0);
std::vector<int>	MuonPFlow_globalTrack_hitPattern_numberOfValidMuonHits(200,0);
std::vector<float>	MuonPFlow_globalTrack_ndof(200,0);
std::vector<float>	MuonPFlow_innerTrack_dz(200,0);
std::vector<int>	MuonPFlow_innerTrack_hitPattern_numberOfValidPixelHits(200,0);
std::vector<int>	MuonPFlow_innerTrack_hitPattern_trackerLayersWithMeasurement(200,0);
std::vector<int>	MuonPFlow_isGlobalMuon(200,0);
std::vector<int>	MuonPFlow_isPFMuon(200,0);
std::vector<int>	MuonPFlow_isStandAloneMuon(200,0);
std::vector<int>	MuonPFlow_isTrackerMuon(200,0);
std::vector<float>	MuonPFlow_neutralHadronIso(200,0);
std::vector<int>	MuonPFlow_numberOfMatchedStations(200,0);
std::vector<float>	MuonPFlow_phi(200,0);
std::vector<float>	MuonPFlow_photonIso(200,0);
std::vector<float>	MuonPFlow_pt(200,0);
std::vector<float>	MuonPFlow_puChargedHadronIso(200,0);
std::vector<float>	MuonPFlow_pz(200,0);
std::vector<float>	MuonPFlow_vertex_z(200,0);
std::vector<int>	Muon_charge(200,0);
std::vector<float>	Muon_chargedHadronIso(200,0);
std::vector<float>	Muon_dB(200,0);
std::vector<float>	Muon_energy(200,0);
std::vector<float>	Muon_et(200,0);
std::vector<float>	Muon_eta(200,0);
std::vector<float>	Muon_globalTrack_chi2(200,0);
std::vector<float>	Muon_globalTrack_d0(200,0);
std::vector<int>	Muon_globalTrack_hitPattern_numberOfValidMuonHits(200,0);
std::vector<float>	Muon_globalTrack_ndof(200,0);
std::vector<float>	Muon_innerTrack_dz(200,0);
std::vector<int>	Muon_innerTrack_hitPattern_numberOfValidPixelHits(200,0);
std::vector<int>	Muon_innerTrack_hitPattern_trackerLayersWithMeasurement(200,0);
std::vector<int>	Muon_isGlobalMuon(200,0);
std::vector<int>	Muon_isPFMuon(200,0);
std::vector<int>	Muon_isStandAloneMuon(200,0);
std::vector<int>	Muon_isTrackerMuon(200,0);
std::vector<float>	Muon_neutralHadronIso(200,0);
std::vector<int>	Muon_numberOfMatchedStations(200,0);
std::vector<float>	Muon_phi(200,0);
std::vector<float>	Muon_photonIso(200,0);
std::vector<float>	Muon_pt(200,0);
std::vector<float>	Muon_puChargedHadronIso(200,0);
std::vector<float>	Muon_pz(200,0);
std::vector<float>	Muon_vertex_z(200,0);
std::vector<int>	PileupSummaryInfo_getBunchCrossing(10,0);
std::vector<int>	PileupSummaryInfo_getPU_NumInteractions(10,0);
std::vector<float>	PileupSummaryInfo_getTrueNumInteractions(10,0);
std::vector<float>	Tau_againstElectronLoose(200,0);
std::vector<float>	Tau_againstMuonTight(200,0);
std::vector<float>	Tau_byLooseCombinedIsolationDeltaBetaCorr(200,0);
std::vector<int>	Tau_charge(200,0);
std::vector<float>	Tau_decayModeFinding(200,0);
std::vector<float>	Tau_energy(200,0);
std::vector<float>	Tau_et(200,0);
std::vector<float>	Tau_eta(200,0);
std::vector<float>	Tau_phi(200,0);
std::vector<float>	Tau_pt(200,0);
std::vector<float>	Tau_pz(200,0);
std::vector<float>	Track_caloEMDeltaRp3(2000,0);
std::vector<float>	Track_caloEMDeltaRp4(2000,0);
std::vector<float>	Track_caloEMDeltaRp5(2000,0);
std::vector<float>	Track_caloHadDeltaRp3(2000,0);
std::vector<float>	Track_caloHadDeltaRp4(2000,0);
std::vector<float>	Track_caloHadDeltaRp5(2000,0);
std::vector<float>	Track_charge(2000,0);
std::vector<float>	Track_chi2(2000,0);
std::vector<float>	Track_eta(2000,0);
std::vector<unsigned short>	Track_hitPattern_trackerLayersWithoutMeasurement(2000,0);
std::vector<float>	Track_ndof(2000,0);
std::vector<unsigned short>	Track_numberOfValidHits(2000,0);
std::vector<float>	Track_phi(2000,0);
std::vector<float>	Track_pt(2000,0);
std::vector<float>	Track_ptError(2000,0);
std::vector<float>	Track_px(2000,0);
std::vector<float>	Track_py(2000,0);
std::vector<float>	Track_pz(2000,0);
std::vector<int>	Track_trackHighPurity(2000,0);
std::vector<float>	Track_trackRelIso03(2000,0);
std::vector<unsigned short>	Track_trackerExpectedHitsInner_numberOfLostHits(2000,0);
std::vector<unsigned short>	Track_trackerExpectedHitsOuter_numberOfHits(2000,0);
std::vector<float>	Track_vx(2000,0);
std::vector<float>	Track_vy(2000,0);
std::vector<float>	Track_vz(2000,0);
std::vector<float>	Vertex_ndof(100,0);
std::vector<float>	Vertex_position_rho(100,0);
std::vector<float>	Vertex_x(100,0);
std::vector<float>	Vertex_y(100,0);
std::vector<float>	Vertex_z(100,0);
int	edmEventHelper_bunchCrossing;
int	edmEventHelper_event;
int	edmEventHelper_isRealData;
int	edmEventHelper_luminosityBlock;
int	edmEventHelper_orbitNumber;
int	edmEventHelper_run;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v21;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v22;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v23;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v24;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v25;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8;
int	edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v1;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v10;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v11;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v12;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v13;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v14;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v15;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v16;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v17;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v18;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v19;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v2;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v20;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v3;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v4;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v5;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v6;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v7;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v8;
int	edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v9;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v1;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v10;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v11;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v12;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v13;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v14;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v15;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v16;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v17;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v18;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v19;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v2;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v20;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v3;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v4;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v5;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v6;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v7;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v8;
int	edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v9;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v1;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v10;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v11;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v12;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v13;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v14;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v15;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v16;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v17;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v18;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v19;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v2;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v20;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v3;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v4;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v5;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v6;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v7;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v8;
int	edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v9;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v1;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v10;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v11;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v12;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v13;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v14;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v15;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v16;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v17;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v18;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v19;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v2;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v20;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v3;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v4;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v5;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v6;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v7;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v8;
int	edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v9;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v21;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v22;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v23;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v24;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v25;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8;
int	edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v1;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v10;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v11;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v12;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v13;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v14;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v15;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v16;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v17;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v18;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v19;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v2;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v20;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v3;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v4;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v5;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v6;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v7;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v8;
int	edmTriggerResultsHelper_value_HLT_Ele27_WP80_v9;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v1;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v10;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v11;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v12;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v13;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v14;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v15;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v16;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v17;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v18;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v19;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v2;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v20;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v3;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v4;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v5;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v6;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v7;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v8;
int	edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v9;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v1;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v10;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v11;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v12;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v13;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v14;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v15;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v16;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v17;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v18;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v19;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v2;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v20;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v3;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v4;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v5;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v6;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v7;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v8;
int	edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v9;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v1;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v10;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v11;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v12;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v13;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v14;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v15;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v16;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v17;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v18;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v19;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v2;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v20;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v3;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v4;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v5;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v6;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v7;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v8;
int	edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v9;
int	nElectron;
int	nElectronPFlow;
int	nGenParticle;
int	nJet;
int	nMuon;
int	nMuonPFlow;
int	nPileupSummaryInfo;
int	nTau;
int	nTrack;
int	nVertex;
double	sdoublePF_value;
float	sdouble_value;
float	weight;

//-----------------------------------------------------------------------------
// --- indexmap keeps track of which objects have been flagged for selection
// --- IMPORTANT: initialize must be called every event to clear selection
std::map<std::string, std::vector<int> > indexmap;
void initialize()
{
  for(std::map<std::string, std::vector<int> >::iterator
    item=indexmap.begin(); 
    item != indexmap.end();
	++item)
	item->second.clear();
}

void select(std::string objname)
{
  indexmap[objname] = std::vector<int>();
}

void select(std::string objname, int index)
{
  try
    {
      indexmap[objname].push_back(index);
    }
  catch (...)
    {
      std::cout << "*** perhaps you failed to call select for " 
                << objname << std::endl;
      assert(0);
    }
}

//-----------------------------------------------------------------------------
// --- Structs can be filled by calling fillObjects()
// --- after the call to stream.read(...)
//-----------------------------------------------------------------------------
struct Electron_s
{
  float	energy;
  float	et;
  float	pz;
  float	pt;
  float	phi;
  float	eta;
  int	charge;
  float	mvaNonTrigV0;
  float	mvaTrigV0;
  int	passConversionVeto;
  int	gsfTrack_trackerExpectedHitsInner_numberOfLostHits;
  float	superCluster_eta;
  float	chargedHadronIso;
  float	neutralHadronIso;
  float	photonIso;
  float	puChargedHadronIso;
  float	Aeff04;
};
std::vector<Electron_s> Electron(200);

std::ostream& operator<<(std::ostream& os, const Electron_s& o)
{
  char r[1024];
  os << "Electron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "mvaNonTrigV0", (double)o.mvaNonTrigV0); os << r;
  sprintf(r, "  %-32s: %f\n", "mvaTrigV0", (double)o.mvaTrigV0); os << r;
  sprintf(r, "  %-32s: %f\n", "passConversionVeto", (double)o.passConversionVeto); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_trackerExpectedHitsInner_numberOfLostHits", (double)o.gsfTrack_trackerExpectedHitsInner_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_eta", (double)o.superCluster_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "puChargedHadronIso", (double)o.puChargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "Aeff04", (double)o.Aeff04); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct ElectronPFlow_s
{
  float	energy;
  float	et;
  float	pz;
  float	pt;
  float	phi;
  float	eta;
  int	charge;
  float	mvaNonTrigV0;
  float	mvaTrigV0;
  int	passConversionVeto;
  int	gsfTrack_trackerExpectedHitsInner_numberOfLostHits;
  float	superCluster_eta;
  float	chargedHadronIso;
  float	neutralHadronIso;
  float	photonIso;
  float	puChargedHadronIso;
  float	Aeff04;
};
std::vector<ElectronPFlow_s> ElectronPFlow(200);

std::ostream& operator<<(std::ostream& os, const ElectronPFlow_s& o)
{
  char r[1024];
  os << "ElectronPFlow" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "mvaNonTrigV0", (double)o.mvaNonTrigV0); os << r;
  sprintf(r, "  %-32s: %f\n", "mvaTrigV0", (double)o.mvaTrigV0); os << r;
  sprintf(r, "  %-32s: %f\n", "passConversionVeto", (double)o.passConversionVeto); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_trackerExpectedHitsInner_numberOfLostHits", (double)o.gsfTrack_trackerExpectedHitsInner_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_eta", (double)o.superCluster_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "puChargedHadronIso", (double)o.puChargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "Aeff04", (double)o.Aeff04); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct GenParticle_s
{
  int	charge;
  float	p;
  float	energy;
  float	et;
  float	pz;
  float	pt;
  float	phi;
  float	eta;
  float	mass;
  int	pdgId;
  int	status;
};
std::vector<GenParticle_s> GenParticle(2000);

std::ostream& operator<<(std::ostream& os, const GenParticle_s& o)
{
  char r[1024];
  os << "GenParticle" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "status", (double)o.status); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Jet_s
{
  float	energy;
  float	et;
  float	pt;
  float	pz;
  float	phi;
  float	eta;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
};
std::vector<Jet_s> Jet(200);

std::ostream& operator<<(std::ostream& os, const Jet_s& o)
{
  char r[1024];
  os << "Jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Muon_s
{
  float	energy;
  float	et;
  float	pz;
  float	pt;
  float	phi;
  float	eta;
  int	charge;
  int	isPFMuon;
  int	isGlobalMuon;
  int	isTrackerMuon;
  int	isStandAloneMuon;
  float	globalTrack_chi2;
  float	globalTrack_ndof;
  float	globalTrack_d0;
  float	dB;
  float	vertex_z;
  float	innerTrack_dz;
  int	globalTrack_hitPattern_numberOfValidMuonHits;
  int	innerTrack_hitPattern_trackerLayersWithMeasurement;
  int	innerTrack_hitPattern_numberOfValidPixelHits;
  int	numberOfMatchedStations;
  float	chargedHadronIso;
  float	neutralHadronIso;
  float	photonIso;
  float	puChargedHadronIso;
  int	pdgId;
  int	status;

};
std::vector<Muon_s> Muon(200);

std::ostream& operator<<(std::ostream& os, const Muon_s& o)
{
  char r[1024];
  os << "Muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "isPFMuon", (double)o.isPFMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isStandAloneMuon", (double)o.isStandAloneMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_chi2", (double)o.globalTrack_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_ndof", (double)o.globalTrack_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_d0", (double)o.globalTrack_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "vertex_z", (double)o.vertex_z); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_dz", (double)o.innerTrack_dz); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_hitPattern_numberOfValidMuonHits", (double)o.globalTrack_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_trackerLayersWithMeasurement", (double)o.innerTrack_hitPattern_trackerLayersWithMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_numberOfValidPixelHits", (double)o.innerTrack_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfMatchedStations", (double)o.numberOfMatchedStations); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "puChargedHadronIso", (double)o.puChargedHadronIso); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct MuonPFlow_s
{
  float	energy;
  float	et;
  float	pz;
  float	pt;
  float	phi;
  float	eta;
  int	charge;
  int	isPFMuon;
  int	isGlobalMuon;
  int	isTrackerMuon;
  int	isStandAloneMuon;
  float	globalTrack_chi2;
  float	globalTrack_ndof;
  float	globalTrack_d0;
  float	dB;
  float	vertex_z;
  float	innerTrack_dz;
  int	globalTrack_hitPattern_numberOfValidMuonHits;
  int	innerTrack_hitPattern_trackerLayersWithMeasurement;
  int	innerTrack_hitPattern_numberOfValidPixelHits;
  int	numberOfMatchedStations;
  float	chargedHadronIso;
  float	neutralHadronIso;
  float	photonIso;
  float	puChargedHadronIso;
  int	pdgId;
  int	status;

};
std::vector<MuonPFlow_s> MuonPFlow(200);

std::ostream& operator<<(std::ostream& os, const MuonPFlow_s& o)
{
  char r[1024];
  os << "MuonPFlow" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "isPFMuon", (double)o.isPFMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isStandAloneMuon", (double)o.isStandAloneMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_chi2", (double)o.globalTrack_chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_ndof", (double)o.globalTrack_ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_d0", (double)o.globalTrack_d0); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "vertex_z", (double)o.vertex_z); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_dz", (double)o.innerTrack_dz); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_hitPattern_numberOfValidMuonHits", (double)o.globalTrack_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_trackerLayersWithMeasurement", (double)o.innerTrack_hitPattern_trackerLayersWithMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_numberOfValidPixelHits", (double)o.innerTrack_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfMatchedStations", (double)o.numberOfMatchedStations); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronIso", (double)o.chargedHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronIso", (double)o.neutralHadronIso); os << r;
  sprintf(r, "  %-32s: %f\n", "photonIso", (double)o.photonIso); os << r;
  sprintf(r, "  %-32s: %f\n", "puChargedHadronIso", (double)o.puChargedHadronIso); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PileupSummaryInfo_s
{
  int	getBunchCrossing;
  int	getPU_NumInteractions;
  float	getTrueNumInteractions;
};
std::vector<PileupSummaryInfo_s> PileupSummaryInfo(10);

std::ostream& operator<<(std::ostream& os, const PileupSummaryInfo_s& o)
{
  char r[1024];
  os << "PileupSummaryInfo" << std::endl;
  sprintf(r, "  %-32s: %f\n", "getBunchCrossing", (double)o.getBunchCrossing); os << r;
  sprintf(r, "  %-32s: %f\n", "getPU_NumInteractions", (double)o.getPU_NumInteractions); os << r;
  sprintf(r, "  %-32s: %f\n", "getTrueNumInteractions", (double)o.getTrueNumInteractions); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Tau_s
{
  float	energy;
  float	et;
  float	pz;
  float	pt;
  float	phi;
  float	eta;
  int	charge;
  float	byLooseCombinedIsolationDeltaBetaCorr;
  float	decayModeFinding;
  float	againstElectronLoose;
  float	againstMuonTight;
};
std::vector<Tau_s> Tau(200);

std::ostream& operator<<(std::ostream& os, const Tau_s& o)
{
  char r[1024];
  os << "Tau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "byLooseCombinedIsolationDeltaBetaCorr", (double)o.byLooseCombinedIsolationDeltaBetaCorr); os << r;
  sprintf(r, "  %-32s: %f\n", "decayModeFinding", (double)o.decayModeFinding); os << r;
  sprintf(r, "  %-32s: %f\n", "againstElectronLoose", (double)o.againstElectronLoose); os << r;
  sprintf(r, "  %-32s: %f\n", "againstMuonTight", (double)o.againstMuonTight); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Track_s
{
  float	pt;
  float	ptError;
  float	px;
  float	py;
  float	pz;
  float	phi;
  float	eta;
  float	vx;
  float	vy;
  float	vz;
  float	chi2;
  float	ndof;
  float	charge;
  unsigned short	numberOfValidHits;
  unsigned short	hitPattern_trackerLayersWithoutMeasurement;
  unsigned short	trackerExpectedHitsInner_numberOfLostHits;
  unsigned short	trackerExpectedHitsOuter_numberOfHits;
  int	trackHighPurity;
  float	trackRelIso03;
  float	caloEMDeltaRp3;
  float	caloHadDeltaRp3;
  float	caloEMDeltaRp4;
  float	caloHadDeltaRp4;
  float	caloEMDeltaRp5;
  float	caloHadDeltaRp5;
  int pdgId;
  int status;
  double genPt;
  double genE;
  double genEt;
  double beta;

};
std::vector<Track_s> Track(2000);

std::ostream& operator<<(std::ostream& os, const Track_s& o)
{
  char r[1024];
  os << "Track" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "ptError", (double)o.ptError); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "vx", (double)o.vx); os << r;
  sprintf(r, "  %-32s: %f\n", "vy", (double)o.vy); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  sprintf(r, "  %-32s: %f\n", "chi2", (double)o.chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfValidHits", (double)o.numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "hitPattern_trackerLayersWithoutMeasurement", (double)o.hitPattern_trackerLayersWithoutMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "trackerExpectedHitsInner_numberOfLostHits", (double)o.trackerExpectedHitsInner_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "trackerExpectedHitsOuter_numberOfHits", (double)o.trackerExpectedHitsOuter_numberOfHits); os << r;
  sprintf(r, "  %-32s: %f\n", "trackHighPurity", (double)o.trackHighPurity); os << r;
  sprintf(r, "  %-32s: %f\n", "trackRelIso03", (double)o.trackRelIso03); os << r;
  sprintf(r, "  %-32s: %f\n", "caloEMDeltaRp3", (double)o.caloEMDeltaRp3); os << r;
  sprintf(r, "  %-32s: %f\n", "caloHadDeltaRp3", (double)o.caloHadDeltaRp3); os << r;
  sprintf(r, "  %-32s: %f\n", "caloEMDeltaRp4", (double)o.caloEMDeltaRp4); os << r;
  sprintf(r, "  %-32s: %f\n", "caloHadDeltaRp4", (double)o.caloHadDeltaRp4); os << r;
  sprintf(r, "  %-32s: %f\n", "caloEMDeltaRp5", (double)o.caloEMDeltaRp5); os << r;
  sprintf(r, "  %-32s: %f\n", "caloHadDeltaRp5", (double)o.caloHadDeltaRp5); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Vertex_s
{
  float	ndof;
  float	x;
  float	y;
  float	z;
  float	position_rho;
};
std::vector<Vertex_s> Vertex(100);

std::ostream& operator<<(std::ostream& os, const Vertex_s& o)
{
  char r[1024];
  os << "Vertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  sprintf(r, "  %-32s: %f\n", "position_rho", (double)o.position_rho); os << r;
  return os;
}
//-----------------------------------------------------------------------------

inline void fillElectron()
{
  Electron.resize(Electron_energy.size());
  for(unsigned int i=0; i < Electron.size(); ++i)
    {
      Electron[i].energy	= Electron_energy[i];
      Electron[i].et	= Electron_et[i];
      Electron[i].pz	= Electron_pz[i];
      Electron[i].pt	= Electron_pt[i];
      Electron[i].phi	= Electron_phi[i];
      Electron[i].eta	= Electron_eta[i];
      Electron[i].charge	= Electron_charge[i];
      Electron[i].mvaNonTrigV0	= Electron_mvaNonTrigV0[i];
      Electron[i].mvaTrigV0	= Electron_mvaTrigV0[i];
      Electron[i].passConversionVeto	= Electron_passConversionVeto[i];
      Electron[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= Electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      Electron[i].superCluster_eta	= Electron_superCluster_eta[i];
      Electron[i].chargedHadronIso	= Electron_chargedHadronIso[i];
      Electron[i].neutralHadronIso	= Electron_neutralHadronIso[i];
      Electron[i].photonIso	= Electron_photonIso[i];
      Electron[i].puChargedHadronIso	= Electron_puChargedHadronIso[i];
      Electron[i].Aeff04	= Electron_Aeff04[i];
    }
}

inline void fillElectronPFlow()
{
  ElectronPFlow.resize(ElectronPFlow_energy.size());
  for(unsigned int i=0; i < ElectronPFlow.size(); ++i)
    {
      ElectronPFlow[i].energy	= ElectronPFlow_energy[i];
      ElectronPFlow[i].et	= ElectronPFlow_et[i];
      ElectronPFlow[i].pz	= ElectronPFlow_pz[i];
      ElectronPFlow[i].pt	= ElectronPFlow_pt[i];
      ElectronPFlow[i].phi	= ElectronPFlow_phi[i];
      ElectronPFlow[i].eta	= ElectronPFlow_eta[i];
      ElectronPFlow[i].charge	= ElectronPFlow_charge[i];
      ElectronPFlow[i].mvaNonTrigV0	= ElectronPFlow_mvaNonTrigV0[i];
      ElectronPFlow[i].mvaTrigV0	= ElectronPFlow_mvaTrigV0[i];
      ElectronPFlow[i].passConversionVeto	= ElectronPFlow_passConversionVeto[i];
      ElectronPFlow[i].gsfTrack_trackerExpectedHitsInner_numberOfLostHits	= ElectronPFlow_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i];
      ElectronPFlow[i].superCluster_eta	= ElectronPFlow_superCluster_eta[i];
      ElectronPFlow[i].chargedHadronIso	= ElectronPFlow_chargedHadronIso[i];
      ElectronPFlow[i].neutralHadronIso	= ElectronPFlow_neutralHadronIso[i];
      ElectronPFlow[i].photonIso	= ElectronPFlow_photonIso[i];
      ElectronPFlow[i].puChargedHadronIso	= ElectronPFlow_puChargedHadronIso[i];
      ElectronPFlow[i].Aeff04	= ElectronPFlow_Aeff04[i];
    }
}

inline void fillGenParticle()
{
  GenParticle.resize(GenParticle_charge.size());
  for(unsigned int i=0; i < GenParticle.size(); ++i)
    {
      GenParticle[i].charge	= GenParticle_charge[i];
      GenParticle[i].p	= GenParticle_p[i];
      GenParticle[i].energy	= GenParticle_energy[i];
      GenParticle[i].et	= GenParticle_et[i];
      GenParticle[i].pz	= GenParticle_pz[i];
      GenParticle[i].pt	= GenParticle_pt[i];
      GenParticle[i].phi	= GenParticle_phi[i];
      GenParticle[i].eta	= GenParticle_eta[i];
      GenParticle[i].mass	= GenParticle_mass[i];
      GenParticle[i].pdgId	= GenParticle_pdgId[i];
      GenParticle[i].status	= GenParticle_status[i];
    }
}

inline void fillJet()
{
  Jet.resize(Jet_energy.size());
  for(unsigned int i=0; i < Jet.size(); ++i)
    {
      Jet[i].energy	= Jet_energy[i];
      Jet[i].et	= Jet_et[i];
      Jet[i].pt	= Jet_pt[i];
      Jet[i].pz	= Jet_pz[i];
      Jet[i].phi	= Jet_phi[i];
      Jet[i].eta	= Jet_eta[i];
      Jet[i].chargedHadronEnergyFraction	= Jet_chargedHadronEnergyFraction[i];
      Jet[i].neutralHadronEnergyFraction	= Jet_neutralHadronEnergyFraction[i];
      Jet[i].chargedEmEnergyFraction	= Jet_chargedEmEnergyFraction[i];
      Jet[i].neutralEmEnergyFraction	= Jet_neutralEmEnergyFraction[i];
    }
}

inline void fillMuon()
{
  Muon.resize(Muon_energy.size());
  for(unsigned int i=0; i < Muon.size(); ++i)
    {
      Muon[i].energy	= Muon_energy[i];
      Muon[i].et	= Muon_et[i];
      Muon[i].pz	= Muon_pz[i];
      Muon[i].pt	= Muon_pt[i];
      Muon[i].phi	= Muon_phi[i];
      Muon[i].eta	= Muon_eta[i];
      Muon[i].charge	= Muon_charge[i];
      Muon[i].isPFMuon	= Muon_isPFMuon[i];
      Muon[i].isGlobalMuon	= Muon_isGlobalMuon[i];
      Muon[i].isTrackerMuon	= Muon_isTrackerMuon[i];
      Muon[i].isStandAloneMuon	= Muon_isStandAloneMuon[i];
      Muon[i].globalTrack_chi2	= Muon_globalTrack_chi2[i];
      Muon[i].globalTrack_ndof	= Muon_globalTrack_ndof[i];
      Muon[i].globalTrack_d0	= Muon_globalTrack_d0[i];
      Muon[i].dB	= Muon_dB[i];
      Muon[i].vertex_z	= Muon_vertex_z[i];
      Muon[i].innerTrack_dz	= Muon_innerTrack_dz[i];
      Muon[i].globalTrack_hitPattern_numberOfValidMuonHits	= Muon_globalTrack_hitPattern_numberOfValidMuonHits[i];
      Muon[i].innerTrack_hitPattern_trackerLayersWithMeasurement	= Muon_innerTrack_hitPattern_trackerLayersWithMeasurement[i];
      Muon[i].innerTrack_hitPattern_numberOfValidPixelHits	= Muon_innerTrack_hitPattern_numberOfValidPixelHits[i];
      Muon[i].numberOfMatchedStations	= Muon_numberOfMatchedStations[i];
      Muon[i].chargedHadronIso	= Muon_chargedHadronIso[i];
      Muon[i].neutralHadronIso	= Muon_neutralHadronIso[i];
      Muon[i].photonIso	= Muon_photonIso[i];
      Muon[i].puChargedHadronIso	= Muon_puChargedHadronIso[i];
    }
}

inline void fillMuonPFlow()
{
  MuonPFlow.resize(MuonPFlow_energy.size());
  for(unsigned int i=0; i < MuonPFlow.size(); ++i)
    {
      MuonPFlow[i].energy	= MuonPFlow_energy[i];
      MuonPFlow[i].et	= MuonPFlow_et[i];
      MuonPFlow[i].pz	= MuonPFlow_pz[i];
      MuonPFlow[i].pt	= MuonPFlow_pt[i];
      MuonPFlow[i].phi	= MuonPFlow_phi[i];
      MuonPFlow[i].eta	= MuonPFlow_eta[i];
      MuonPFlow[i].charge	= MuonPFlow_charge[i];
      MuonPFlow[i].isPFMuon	= MuonPFlow_isPFMuon[i];
      MuonPFlow[i].isGlobalMuon	= MuonPFlow_isGlobalMuon[i];
      MuonPFlow[i].isTrackerMuon	= MuonPFlow_isTrackerMuon[i];
      MuonPFlow[i].isStandAloneMuon	= MuonPFlow_isStandAloneMuon[i];
      MuonPFlow[i].globalTrack_chi2	= MuonPFlow_globalTrack_chi2[i];
      MuonPFlow[i].globalTrack_ndof	= MuonPFlow_globalTrack_ndof[i];
      MuonPFlow[i].globalTrack_d0	= MuonPFlow_globalTrack_d0[i];
      MuonPFlow[i].dB	= MuonPFlow_dB[i];
      MuonPFlow[i].vertex_z	= MuonPFlow_vertex_z[i];
      MuonPFlow[i].innerTrack_dz	= MuonPFlow_innerTrack_dz[i];
      MuonPFlow[i].globalTrack_hitPattern_numberOfValidMuonHits	= MuonPFlow_globalTrack_hitPattern_numberOfValidMuonHits[i];
      MuonPFlow[i].innerTrack_hitPattern_trackerLayersWithMeasurement	= MuonPFlow_innerTrack_hitPattern_trackerLayersWithMeasurement[i];
      MuonPFlow[i].innerTrack_hitPattern_numberOfValidPixelHits	= MuonPFlow_innerTrack_hitPattern_numberOfValidPixelHits[i];
      MuonPFlow[i].numberOfMatchedStations	= MuonPFlow_numberOfMatchedStations[i];
      MuonPFlow[i].chargedHadronIso	= MuonPFlow_chargedHadronIso[i];
      MuonPFlow[i].neutralHadronIso	= MuonPFlow_neutralHadronIso[i];
      MuonPFlow[i].photonIso	= MuonPFlow_photonIso[i];
      MuonPFlow[i].puChargedHadronIso	= MuonPFlow_puChargedHadronIso[i];
    }
}

inline void fillPileupSummaryInfo()
{
  PileupSummaryInfo.resize(PileupSummaryInfo_getBunchCrossing.size());
  for(unsigned int i=0; i < PileupSummaryInfo.size(); ++i)
    {
      PileupSummaryInfo[i].getBunchCrossing	= PileupSummaryInfo_getBunchCrossing[i];
      PileupSummaryInfo[i].getPU_NumInteractions	= PileupSummaryInfo_getPU_NumInteractions[i];
      PileupSummaryInfo[i].getTrueNumInteractions	= PileupSummaryInfo_getTrueNumInteractions[i];
    }
}

inline void fillTau()
{
  Tau.resize(Tau_energy.size());
  for(unsigned int i=0; i < Tau.size(); ++i)
    {
      Tau[i].energy	= Tau_energy[i];
      Tau[i].et	= Tau_et[i];
      Tau[i].pz	= Tau_pz[i];
      Tau[i].pt	= Tau_pt[i];
      Tau[i].phi	= Tau_phi[i];
      Tau[i].eta	= Tau_eta[i];
      Tau[i].charge	= Tau_charge[i];
      Tau[i].byLooseCombinedIsolationDeltaBetaCorr	= Tau_byLooseCombinedIsolationDeltaBetaCorr[i];
      Tau[i].decayModeFinding	= Tau_decayModeFinding[i];
      Tau[i].againstElectronLoose	= Tau_againstElectronLoose[i];
      Tau[i].againstMuonTight	= Tau_againstMuonTight[i];
    }
}

inline void fillTrack()
{
  Track.resize(Track_pt.size());
  for(unsigned int i=0; i < Track.size(); ++i)
    {
      Track[i].pt	= Track_pt[i];
      Track[i].ptError	= Track_ptError[i];
      Track[i].px	= Track_px[i];
      Track[i].py	= Track_py[i];
      Track[i].pz	= Track_pz[i];
      Track[i].phi	= Track_phi[i];
      Track[i].eta	= Track_eta[i];
      Track[i].vx	= Track_vx[i];
      Track[i].vy	= Track_vy[i];
      Track[i].vz	= Track_vz[i];
      Track[i].chi2	= Track_chi2[i];
      Track[i].ndof	= Track_ndof[i];
      Track[i].charge	= Track_charge[i];
      Track[i].numberOfValidHits	= Track_numberOfValidHits[i];
      Track[i].hitPattern_trackerLayersWithoutMeasurement	= Track_hitPattern_trackerLayersWithoutMeasurement[i];
      Track[i].trackerExpectedHitsInner_numberOfLostHits	= Track_trackerExpectedHitsInner_numberOfLostHits[i];
      Track[i].trackerExpectedHitsOuter_numberOfHits	= Track_trackerExpectedHitsOuter_numberOfHits[i];
      Track[i].trackHighPurity	= Track_trackHighPurity[i];
      Track[i].trackRelIso03	= Track_trackRelIso03[i];
      Track[i].caloEMDeltaRp3	= Track_caloEMDeltaRp3[i];
      Track[i].caloHadDeltaRp3	= Track_caloHadDeltaRp3[i];
      Track[i].caloEMDeltaRp4	= Track_caloEMDeltaRp4[i];
      Track[i].caloHadDeltaRp4	= Track_caloHadDeltaRp4[i];
      Track[i].caloEMDeltaRp5	= Track_caloEMDeltaRp5[i];
      Track[i].caloHadDeltaRp5	= Track_caloHadDeltaRp5[i];
    }
}

inline void fillVertex()
{
  Vertex.resize(Vertex_ndof.size());
  for(unsigned int i=0; i < Vertex.size(); ++i)
    {
      Vertex[i].ndof	= Vertex_ndof[i];
      Vertex[i].x	= Vertex_x[i];
      Vertex[i].y	= Vertex_y[i];
      Vertex[i].z	= Vertex_z[i];
      Vertex[i].position_rho	= Vertex_position_rho[i];
    }
}


void fillObjects()
{
  fillElectron();
  fillElectronPFlow();
  fillGenParticle();
  fillJet();
  fillMuon();
  fillMuonPFlow();
  fillPileupSummaryInfo();
  fillTau();
  fillTrack();
  fillVertex();
}

//-----------------------------------------------------------------------------
// --- Call saveSelectedObjects() just before call to addEvent if
// --- you wish to save only the selected objects
//-----------------------------------------------------------------------------
// Select objects for which the select function was called
void saveSelectedObjects()
{
  int n = 0;

  n = 0;
  try
    {
       n = indexmap["Electron"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Electron"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Electron_energy[i]	= Electron_energy[j];
          Electron_et[i]	= Electron_et[j];
          Electron_pz[i]	= Electron_pz[j];
          Electron_pt[i]	= Electron_pt[j];
          Electron_phi[i]	= Electron_phi[j];
          Electron_eta[i]	= Electron_eta[j];
          Electron_charge[i]	= Electron_charge[j];
          Electron_mvaNonTrigV0[i]	= Electron_mvaNonTrigV0[j];
          Electron_mvaTrigV0[i]	= Electron_mvaTrigV0[j];
          Electron_passConversionVeto[i]	= Electron_passConversionVeto[j];
          Electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i]	= Electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[j];
          Electron_superCluster_eta[i]	= Electron_superCluster_eta[j];
          Electron_chargedHadronIso[i]	= Electron_chargedHadronIso[j];
          Electron_neutralHadronIso[i]	= Electron_neutralHadronIso[j];
          Electron_photonIso[i]	= Electron_photonIso[j];
          Electron_puChargedHadronIso[i]	= Electron_puChargedHadronIso[j];
          Electron_Aeff04[i]	= Electron_Aeff04[j];
        }
      nElectron = n;
    }

  n = 0;
  try
    {
       n = indexmap["ElectronPFlow"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["ElectronPFlow"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          ElectronPFlow_energy[i]	= ElectronPFlow_energy[j];
          ElectronPFlow_et[i]	= ElectronPFlow_et[j];
          ElectronPFlow_pz[i]	= ElectronPFlow_pz[j];
          ElectronPFlow_pt[i]	= ElectronPFlow_pt[j];
          ElectronPFlow_phi[i]	= ElectronPFlow_phi[j];
          ElectronPFlow_eta[i]	= ElectronPFlow_eta[j];
          ElectronPFlow_charge[i]	= ElectronPFlow_charge[j];
          ElectronPFlow_mvaNonTrigV0[i]	= ElectronPFlow_mvaNonTrigV0[j];
          ElectronPFlow_mvaTrigV0[i]	= ElectronPFlow_mvaTrigV0[j];
          ElectronPFlow_passConversionVeto[i]	= ElectronPFlow_passConversionVeto[j];
          ElectronPFlow_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[i]	= ElectronPFlow_gsfTrack_trackerExpectedHitsInner_numberOfLostHits[j];
          ElectronPFlow_superCluster_eta[i]	= ElectronPFlow_superCluster_eta[j];
          ElectronPFlow_chargedHadronIso[i]	= ElectronPFlow_chargedHadronIso[j];
          ElectronPFlow_neutralHadronIso[i]	= ElectronPFlow_neutralHadronIso[j];
          ElectronPFlow_photonIso[i]	= ElectronPFlow_photonIso[j];
          ElectronPFlow_puChargedHadronIso[i]	= ElectronPFlow_puChargedHadronIso[j];
          ElectronPFlow_Aeff04[i]	= ElectronPFlow_Aeff04[j];
        }
      nElectronPFlow = n;
    }

  n = 0;
  try
    {
       n = indexmap["GenParticle"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["GenParticle"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          GenParticle_charge[i]	= GenParticle_charge[j];
          GenParticle_p[i]	= GenParticle_p[j];
          GenParticle_energy[i]	= GenParticle_energy[j];
          GenParticle_et[i]	= GenParticle_et[j];
          GenParticle_pz[i]	= GenParticle_pz[j];
          GenParticle_pt[i]	= GenParticle_pt[j];
          GenParticle_phi[i]	= GenParticle_phi[j];
          GenParticle_eta[i]	= GenParticle_eta[j];
          GenParticle_mass[i]	= GenParticle_mass[j];
          GenParticle_pdgId[i]	= GenParticle_pdgId[j];
          GenParticle_status[i]	= GenParticle_status[j];
        }
      nGenParticle = n;
    }

  n = 0;
  try
    {
       n = indexmap["Jet"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Jet"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Jet_energy[i]	= Jet_energy[j];
          Jet_et[i]	= Jet_et[j];
          Jet_pt[i]	= Jet_pt[j];
          Jet_pz[i]	= Jet_pz[j];
          Jet_phi[i]	= Jet_phi[j];
          Jet_eta[i]	= Jet_eta[j];
          Jet_chargedHadronEnergyFraction[i]	= Jet_chargedHadronEnergyFraction[j];
          Jet_neutralHadronEnergyFraction[i]	= Jet_neutralHadronEnergyFraction[j];
          Jet_chargedEmEnergyFraction[i]	= Jet_chargedEmEnergyFraction[j];
          Jet_neutralEmEnergyFraction[i]	= Jet_neutralEmEnergyFraction[j];
        }
      nJet = n;
    }

  n = 0;
  try
    {
       n = indexmap["Muon"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Muon"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Muon_energy[i]	= Muon_energy[j];
          Muon_et[i]	= Muon_et[j];
          Muon_pz[i]	= Muon_pz[j];
          Muon_pt[i]	= Muon_pt[j];
          Muon_phi[i]	= Muon_phi[j];
          Muon_eta[i]	= Muon_eta[j];
          Muon_charge[i]	= Muon_charge[j];
          Muon_isPFMuon[i]	= Muon_isPFMuon[j];
          Muon_isGlobalMuon[i]	= Muon_isGlobalMuon[j];
          Muon_isTrackerMuon[i]	= Muon_isTrackerMuon[j];
          Muon_isStandAloneMuon[i]	= Muon_isStandAloneMuon[j];
          Muon_globalTrack_chi2[i]	= Muon_globalTrack_chi2[j];
          Muon_globalTrack_ndof[i]	= Muon_globalTrack_ndof[j];
          Muon_globalTrack_d0[i]	= Muon_globalTrack_d0[j];
          Muon_dB[i]	= Muon_dB[j];
          Muon_vertex_z[i]	= Muon_vertex_z[j];
          Muon_innerTrack_dz[i]	= Muon_innerTrack_dz[j];
          Muon_globalTrack_hitPattern_numberOfValidMuonHits[i]	= Muon_globalTrack_hitPattern_numberOfValidMuonHits[j];
          Muon_innerTrack_hitPattern_trackerLayersWithMeasurement[i]	= Muon_innerTrack_hitPattern_trackerLayersWithMeasurement[j];
          Muon_innerTrack_hitPattern_numberOfValidPixelHits[i]	= Muon_innerTrack_hitPattern_numberOfValidPixelHits[j];
          Muon_numberOfMatchedStations[i]	= Muon_numberOfMatchedStations[j];
          Muon_chargedHadronIso[i]	= Muon_chargedHadronIso[j];
          Muon_neutralHadronIso[i]	= Muon_neutralHadronIso[j];
          Muon_photonIso[i]	= Muon_photonIso[j];
          Muon_puChargedHadronIso[i]	= Muon_puChargedHadronIso[j];
        }
      nMuon = n;
    }

  n = 0;
  try
    {
       n = indexmap["MuonPFlow"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["MuonPFlow"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          MuonPFlow_energy[i]	= MuonPFlow_energy[j];
          MuonPFlow_et[i]	= MuonPFlow_et[j];
          MuonPFlow_pz[i]	= MuonPFlow_pz[j];
          MuonPFlow_pt[i]	= MuonPFlow_pt[j];
          MuonPFlow_phi[i]	= MuonPFlow_phi[j];
          MuonPFlow_eta[i]	= MuonPFlow_eta[j];
          MuonPFlow_charge[i]	= MuonPFlow_charge[j];
          MuonPFlow_isPFMuon[i]	= MuonPFlow_isPFMuon[j];
          MuonPFlow_isGlobalMuon[i]	= MuonPFlow_isGlobalMuon[j];
          MuonPFlow_isTrackerMuon[i]	= MuonPFlow_isTrackerMuon[j];
          MuonPFlow_isStandAloneMuon[i]	= MuonPFlow_isStandAloneMuon[j];
          MuonPFlow_globalTrack_chi2[i]	= MuonPFlow_globalTrack_chi2[j];
          MuonPFlow_globalTrack_ndof[i]	= MuonPFlow_globalTrack_ndof[j];
          MuonPFlow_globalTrack_d0[i]	= MuonPFlow_globalTrack_d0[j];
          MuonPFlow_dB[i]	= MuonPFlow_dB[j];
          MuonPFlow_vertex_z[i]	= MuonPFlow_vertex_z[j];
          MuonPFlow_innerTrack_dz[i]	= MuonPFlow_innerTrack_dz[j];
          MuonPFlow_globalTrack_hitPattern_numberOfValidMuonHits[i]	= MuonPFlow_globalTrack_hitPattern_numberOfValidMuonHits[j];
          MuonPFlow_innerTrack_hitPattern_trackerLayersWithMeasurement[i]	= MuonPFlow_innerTrack_hitPattern_trackerLayersWithMeasurement[j];
          MuonPFlow_innerTrack_hitPattern_numberOfValidPixelHits[i]	= MuonPFlow_innerTrack_hitPattern_numberOfValidPixelHits[j];
          MuonPFlow_numberOfMatchedStations[i]	= MuonPFlow_numberOfMatchedStations[j];
          MuonPFlow_chargedHadronIso[i]	= MuonPFlow_chargedHadronIso[j];
          MuonPFlow_neutralHadronIso[i]	= MuonPFlow_neutralHadronIso[j];
          MuonPFlow_photonIso[i]	= MuonPFlow_photonIso[j];
          MuonPFlow_puChargedHadronIso[i]	= MuonPFlow_puChargedHadronIso[j];
        }
      nMuonPFlow = n;
    }

  n = 0;
  try
    {
       n = indexmap["PileupSummaryInfo"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PileupSummaryInfo"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PileupSummaryInfo_getBunchCrossing[i]	= PileupSummaryInfo_getBunchCrossing[j];
          PileupSummaryInfo_getPU_NumInteractions[i]	= PileupSummaryInfo_getPU_NumInteractions[j];
          PileupSummaryInfo_getTrueNumInteractions[i]	= PileupSummaryInfo_getTrueNumInteractions[j];
        }
      nPileupSummaryInfo = n;
    }

  n = 0;
  try
    {
       n = indexmap["Tau"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Tau"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Tau_energy[i]	= Tau_energy[j];
          Tau_et[i]	= Tau_et[j];
          Tau_pz[i]	= Tau_pz[j];
          Tau_pt[i]	= Tau_pt[j];
          Tau_phi[i]	= Tau_phi[j];
          Tau_eta[i]	= Tau_eta[j];
          Tau_charge[i]	= Tau_charge[j];
          Tau_byLooseCombinedIsolationDeltaBetaCorr[i]	= Tau_byLooseCombinedIsolationDeltaBetaCorr[j];
          Tau_decayModeFinding[i]	= Tau_decayModeFinding[j];
          Tau_againstElectronLoose[i]	= Tau_againstElectronLoose[j];
          Tau_againstMuonTight[i]	= Tau_againstMuonTight[j];
        }
      nTau = n;
    }

  n = 0;
  try
    {
       n = indexmap["Track"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Track"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Track_pt[i]	= Track_pt[j];
          Track_ptError[i]	= Track_ptError[j];
          Track_px[i]	= Track_px[j];
          Track_py[i]	= Track_py[j];
          Track_pz[i]	= Track_pz[j];
          Track_phi[i]	= Track_phi[j];
          Track_eta[i]	= Track_eta[j];
          Track_vx[i]	= Track_vx[j];
          Track_vy[i]	= Track_vy[j];
          Track_vz[i]	= Track_vz[j];
          Track_chi2[i]	= Track_chi2[j];
          Track_ndof[i]	= Track_ndof[j];
          Track_charge[i]	= Track_charge[j];
          Track_numberOfValidHits[i]	= Track_numberOfValidHits[j];
          Track_hitPattern_trackerLayersWithoutMeasurement[i]	= Track_hitPattern_trackerLayersWithoutMeasurement[j];
          Track_trackerExpectedHitsInner_numberOfLostHits[i]	= Track_trackerExpectedHitsInner_numberOfLostHits[j];
          Track_trackerExpectedHitsOuter_numberOfHits[i]	= Track_trackerExpectedHitsOuter_numberOfHits[j];
          Track_trackHighPurity[i]	= Track_trackHighPurity[j];
          Track_trackRelIso03[i]	= Track_trackRelIso03[j];
          Track_caloEMDeltaRp3[i]	= Track_caloEMDeltaRp3[j];
          Track_caloHadDeltaRp3[i]	= Track_caloHadDeltaRp3[j];
          Track_caloEMDeltaRp4[i]	= Track_caloEMDeltaRp4[j];
          Track_caloHadDeltaRp4[i]	= Track_caloHadDeltaRp4[j];
          Track_caloEMDeltaRp5[i]	= Track_caloEMDeltaRp5[j];
          Track_caloHadDeltaRp5[i]	= Track_caloHadDeltaRp5[j];
        }
      nTrack = n;
    }

  n = 0;
  try
    {
       n = indexmap["Vertex"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Vertex"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Vertex_ndof[i]	= Vertex_ndof[j];
          Vertex_x[i]	= Vertex_x[j];
          Vertex_y[i]	= Vertex_y[j];
          Vertex_z[i]	= Vertex_z[j];
          Vertex_position_rho[i]	= Vertex_position_rho[j];
        }
      nVertex = n;
    }
}

//-----------------------------------------------------------------------------
// --- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("patElectronHelper_patElectronsLoosePFlow.Aeff04", ElectronPFlow_Aeff04);
  stream.select("patElectronHelper_patElectronsLoosePFlow.charge", ElectronPFlow_charge);
  stream.select("patElectronHelper_patElectronsLoosePFlow.chargedHadronIso", ElectronPFlow_chargedHadronIso);
  stream.select("patElectronHelper_patElectronsLoosePFlow.energy", ElectronPFlow_energy);
  stream.select("patElectronHelper_patElectronsLoosePFlow.et", ElectronPFlow_et);
  stream.select("patElectronHelper_patElectronsLoosePFlow.eta", ElectronPFlow_eta);
  stream.select("patElectronHelper_patElectronsLoosePFlow.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", ElectronPFlow_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectronHelper_patElectronsLoosePFlow.mvaNonTrigV0", ElectronPFlow_mvaNonTrigV0);
  stream.select("patElectronHelper_patElectronsLoosePFlow.mvaTrigV0", ElectronPFlow_mvaTrigV0);
  stream.select("patElectronHelper_patElectronsLoosePFlow.neutralHadronIso", ElectronPFlow_neutralHadronIso);
  stream.select("patElectronHelper_patElectronsLoosePFlow.passConversionVeto", ElectronPFlow_passConversionVeto);
  stream.select("patElectronHelper_patElectronsLoosePFlow.phi", ElectronPFlow_phi);
  stream.select("patElectronHelper_patElectronsLoosePFlow.photonIso", ElectronPFlow_photonIso);
  stream.select("patElectronHelper_patElectronsLoosePFlow.pt", ElectronPFlow_pt);
  stream.select("patElectronHelper_patElectronsLoosePFlow.puChargedHadronIso", ElectronPFlow_puChargedHadronIso);
  stream.select("patElectronHelper_patElectronsLoosePFlow.pz", ElectronPFlow_pz);
  stream.select("patElectronHelper_patElectronsLoosePFlow.superCluster_eta", ElectronPFlow_superCluster_eta);
  stream.select("patElectronHelper_selectedPatElectrons.Aeff04", Electron_Aeff04);
  stream.select("patElectronHelper_selectedPatElectrons.charge", Electron_charge);
  stream.select("patElectronHelper_selectedPatElectrons.chargedHadronIso", Electron_chargedHadronIso);
  stream.select("patElectronHelper_selectedPatElectrons.energy", Electron_energy);
  stream.select("patElectronHelper_selectedPatElectrons.et", Electron_et);
  stream.select("patElectronHelper_selectedPatElectrons.eta", Electron_eta);
  stream.select("patElectronHelper_selectedPatElectrons.gsfTrack_trackerExpectedHitsInner_numberOfLostHits", Electron_gsfTrack_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("patElectronHelper_selectedPatElectrons.mvaNonTrigV0", Electron_mvaNonTrigV0);
  stream.select("patElectronHelper_selectedPatElectrons.mvaTrigV0", Electron_mvaTrigV0);
  stream.select("patElectronHelper_selectedPatElectrons.neutralHadronIso", Electron_neutralHadronIso);
  stream.select("patElectronHelper_selectedPatElectrons.passConversionVeto", Electron_passConversionVeto);
  stream.select("patElectronHelper_selectedPatElectrons.phi", Electron_phi);
  stream.select("patElectronHelper_selectedPatElectrons.photonIso", Electron_photonIso);
  stream.select("patElectronHelper_selectedPatElectrons.pt", Electron_pt);
  stream.select("patElectronHelper_selectedPatElectrons.puChargedHadronIso", Electron_puChargedHadronIso);
  stream.select("patElectronHelper_selectedPatElectrons.pz", Electron_pz);
  stream.select("patElectronHelper_selectedPatElectrons.superCluster_eta", Electron_superCluster_eta);
  stream.select("recoGenParticle_genParticlesReduced.charge", GenParticle_charge);
  stream.select("recoGenParticle_genParticlesReduced.energy", GenParticle_energy);
  stream.select("recoGenParticle_genParticlesReduced.et", GenParticle_et);
  stream.select("recoGenParticle_genParticlesReduced.eta", GenParticle_eta);
  stream.select("recoGenParticle_genParticlesReduced.mass", GenParticle_mass);
  stream.select("recoGenParticle_genParticlesReduced.p", GenParticle_p);
  stream.select("recoGenParticle_genParticlesReduced.pdgId", GenParticle_pdgId);
  stream.select("recoGenParticle_genParticlesReduced.phi", GenParticle_phi);
  stream.select("recoGenParticle_genParticlesReduced.pt", GenParticle_pt);
  stream.select("recoGenParticle_genParticlesReduced.pz", GenParticle_pz);
  stream.select("recoGenParticle_genParticlesReduced.status", GenParticle_status);
  stream.select("patJet_selectedPatJetsPFlow.chargedEmEnergyFraction", Jet_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.chargedHadronEnergyFraction", Jet_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.energy", Jet_energy);
  stream.select("patJet_selectedPatJetsPFlow.et", Jet_et);
  stream.select("patJet_selectedPatJetsPFlow.eta", Jet_eta);
  stream.select("patJet_selectedPatJetsPFlow.neutralEmEnergyFraction", Jet_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.neutralHadronEnergyFraction", Jet_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.phi", Jet_phi);
  stream.select("patJet_selectedPatJetsPFlow.pt", Jet_pt);
  stream.select("patJet_selectedPatJetsPFlow.pz", Jet_pz);
  stream.select("patMET_patMETsPFlow.energy", MET_energy);
  stream.select("patMET_patMETsPFlow.et", MET_et);
  stream.select("patMET_patMETsPFlow.eta", MET_eta);
  stream.select("patMET_patMETsPFlow.phi", MET_phi);
  stream.select("patMET_patMETsPFlow.pt", MET_pt);
  stream.select("patMET_patMETsPFlow.pz", MET_pz);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.charge", MuonPFlow_charge);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.chargedHadronIso", MuonPFlow_chargedHadronIso);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.dB", MuonPFlow_dB);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.energy", MuonPFlow_energy);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.et", MuonPFlow_et);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.eta", MuonPFlow_eta);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.globalTrack_chi2", MuonPFlow_globalTrack_chi2);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.globalTrack_d0", MuonPFlow_globalTrack_d0);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.globalTrack_hitPattern_numberOfValidMuonHits", MuonPFlow_globalTrack_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.globalTrack_ndof", MuonPFlow_globalTrack_ndof);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.innerTrack_dz", MuonPFlow_innerTrack_dz);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.innerTrack_hitPattern_numberOfValidPixelHits", MuonPFlow_innerTrack_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.innerTrack_hitPattern_trackerLayersWithMeasurement", MuonPFlow_innerTrack_hitPattern_trackerLayersWithMeasurement);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.isGlobalMuon", MuonPFlow_isGlobalMuon);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.isPFMuon", MuonPFlow_isPFMuon);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.isStandAloneMuon", MuonPFlow_isStandAloneMuon);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.isTrackerMuon", MuonPFlow_isTrackerMuon);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.neutralHadronIso", MuonPFlow_neutralHadronIso);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.numberOfMatchedStations", MuonPFlow_numberOfMatchedStations);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.phi", MuonPFlow_phi);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.photonIso", MuonPFlow_photonIso);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.pt", MuonPFlow_pt);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.puChargedHadronIso", MuonPFlow_puChargedHadronIso);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.pz", MuonPFlow_pz);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.vertex_z", MuonPFlow_vertex_z);
  stream.select("patMuon_selectedPatMuonsLoose.charge", Muon_charge);
  stream.select("patMuon_selectedPatMuonsLoose.chargedHadronIso", Muon_chargedHadronIso);
  stream.select("patMuon_selectedPatMuonsLoose.dB", Muon_dB);
  stream.select("patMuon_selectedPatMuonsLoose.energy", Muon_energy);
  stream.select("patMuon_selectedPatMuonsLoose.et", Muon_et);
  stream.select("patMuon_selectedPatMuonsLoose.eta", Muon_eta);
  stream.select("patMuon_selectedPatMuonsLoose.globalTrack_chi2", Muon_globalTrack_chi2);
  stream.select("patMuon_selectedPatMuonsLoose.globalTrack_d0", Muon_globalTrack_d0);
  stream.select("patMuon_selectedPatMuonsLoose.globalTrack_hitPattern_numberOfValidMuonHits", Muon_globalTrack_hitPattern_numberOfValidMuonHits);
  stream.select("patMuon_selectedPatMuonsLoose.globalTrack_ndof", Muon_globalTrack_ndof);
  stream.select("patMuon_selectedPatMuonsLoose.innerTrack_dz", Muon_innerTrack_dz);
  stream.select("patMuon_selectedPatMuonsLoose.innerTrack_hitPattern_numberOfValidPixelHits", Muon_innerTrack_hitPattern_numberOfValidPixelHits);
  stream.select("patMuon_selectedPatMuonsLoose.innerTrack_hitPattern_trackerLayersWithMeasurement", Muon_innerTrack_hitPattern_trackerLayersWithMeasurement);
  stream.select("patMuon_selectedPatMuonsLoose.isGlobalMuon", Muon_isGlobalMuon);
  stream.select("patMuon_selectedPatMuonsLoose.isPFMuon", Muon_isPFMuon);
  stream.select("patMuon_selectedPatMuonsLoose.isStandAloneMuon", Muon_isStandAloneMuon);
  stream.select("patMuon_selectedPatMuonsLoose.isTrackerMuon", Muon_isTrackerMuon);
  stream.select("patMuon_selectedPatMuonsLoose.neutralHadronIso", Muon_neutralHadronIso);
  stream.select("patMuon_selectedPatMuonsLoose.numberOfMatchedStations", Muon_numberOfMatchedStations);
  stream.select("patMuon_selectedPatMuonsLoose.phi", Muon_phi);
  stream.select("patMuon_selectedPatMuonsLoose.photonIso", Muon_photonIso);
  stream.select("patMuon_selectedPatMuonsLoose.pt", Muon_pt);
  stream.select("patMuon_selectedPatMuonsLoose.puChargedHadronIso", Muon_puChargedHadronIso);
  stream.select("patMuon_selectedPatMuonsLoose.pz", Muon_pz);
  stream.select("patMuon_selectedPatMuonsLoose.vertex_z", Muon_vertex_z);
  stream.select("PileupSummaryInfo_addPileupInfo.getBunchCrossing", PileupSummaryInfo_getBunchCrossing);
  stream.select("PileupSummaryInfo_addPileupInfo.getPU_NumInteractions", PileupSummaryInfo_getPU_NumInteractions);
  stream.select("PileupSummaryInfo_addPileupInfo.getTrueNumInteractions", PileupSummaryInfo_getTrueNumInteractions);
  stream.select("patTau_selectedPatTaus.againstElectronLoose", Tau_againstElectronLoose);
  stream.select("patTau_selectedPatTaus.againstMuonTight", Tau_againstMuonTight);
  stream.select("patTau_selectedPatTaus.byLooseCombinedIsolationDeltaBetaCorr", Tau_byLooseCombinedIsolationDeltaBetaCorr);
  stream.select("patTau_selectedPatTaus.charge", Tau_charge);
  stream.select("patTau_selectedPatTaus.decayModeFinding", Tau_decayModeFinding);
  stream.select("patTau_selectedPatTaus.energy", Tau_energy);
  stream.select("patTau_selectedPatTaus.et", Tau_et);
  stream.select("patTau_selectedPatTaus.eta", Tau_eta);
  stream.select("patTau_selectedPatTaus.phi", Tau_phi);
  stream.select("patTau_selectedPatTaus.pt", Tau_pt);
  stream.select("patTau_selectedPatTaus.pz", Tau_pz);
  stream.select("recoTrackHelper_generalTracksReduced.caloEMDeltaRp3", Track_caloEMDeltaRp3);
  stream.select("recoTrackHelper_generalTracksReduced.caloEMDeltaRp4", Track_caloEMDeltaRp4);
  stream.select("recoTrackHelper_generalTracksReduced.caloEMDeltaRp5", Track_caloEMDeltaRp5);
  stream.select("recoTrackHelper_generalTracksReduced.caloHadDeltaRp3", Track_caloHadDeltaRp3);
  stream.select("recoTrackHelper_generalTracksReduced.caloHadDeltaRp4", Track_caloHadDeltaRp4);
  stream.select("recoTrackHelper_generalTracksReduced.caloHadDeltaRp5", Track_caloHadDeltaRp5);
  stream.select("recoTrackHelper_generalTracksReduced.charge", Track_charge);
  stream.select("recoTrackHelper_generalTracksReduced.chi2", Track_chi2);
  stream.select("recoTrackHelper_generalTracksReduced.eta", Track_eta);
  stream.select("recoTrackHelper_generalTracksReduced.hitPattern_trackerLayersWithoutMeasurement", Track_hitPattern_trackerLayersWithoutMeasurement);
  stream.select("recoTrackHelper_generalTracksReduced.ndof", Track_ndof);
  stream.select("recoTrackHelper_generalTracksReduced.numberOfValidHits", Track_numberOfValidHits);
  stream.select("recoTrackHelper_generalTracksReduced.phi", Track_phi);
  stream.select("recoTrackHelper_generalTracksReduced.pt", Track_pt);
  stream.select("recoTrackHelper_generalTracksReduced.ptError", Track_ptError);
  stream.select("recoTrackHelper_generalTracksReduced.px", Track_px);
  stream.select("recoTrackHelper_generalTracksReduced.py", Track_py);
  stream.select("recoTrackHelper_generalTracksReduced.pz", Track_pz);
  stream.select("recoTrackHelper_generalTracksReduced.trackHighPurity", Track_trackHighPurity);
  stream.select("recoTrackHelper_generalTracksReduced.trackRelIso03", Track_trackRelIso03);
  stream.select("recoTrackHelper_generalTracksReduced.trackerExpectedHitsInner_numberOfLostHits", Track_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("recoTrackHelper_generalTracksReduced.trackerExpectedHitsOuter_numberOfHits", Track_trackerExpectedHitsOuter_numberOfHits);
  stream.select("recoTrackHelper_generalTracksReduced.vx", Track_vx);
  stream.select("recoTrackHelper_generalTracksReduced.vy", Track_vy);
  stream.select("recoTrackHelper_generalTracksReduced.vz", Track_vz);
  stream.select("recoVertex_offlinePrimaryVertices.ndof", Vertex_ndof);
  stream.select("recoVertex_offlinePrimaryVertices.position_rho", Vertex_position_rho);
  stream.select("recoVertex_offlinePrimaryVertices.x", Vertex_x);
  stream.select("recoVertex_offlinePrimaryVertices.y", Vertex_y);
  stream.select("recoVertex_offlinePrimaryVertices.z", Vertex_z);
  stream.select("edmEventHelper_info.bunchCrossing", edmEventHelper_bunchCrossing);
  stream.select("edmEventHelper_info.event", edmEventHelper_event);
  stream.select("edmEventHelper_info.isRealData", edmEventHelper_isRealData);
  stream.select("edmEventHelper_info.luminosityBlock", edmEventHelper_luminosityBlock);
  stream.select("edmEventHelper_info.orbitNumber", edmEventHelper_orbitNumber);
  stream.select("edmEventHelper_info.run", edmEventHelper_run);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v21", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v21);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v22", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v22);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v23", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v23);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v24", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v24);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v25", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v25);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9", edmTriggerResultsHelper_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v1", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v10", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v11", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v12", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v13", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v14", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v15", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v16", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v17", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v18", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v19", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v2", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v20", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v3", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v4", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v5", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v6", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v7", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v8", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Ele27_WP80_v9", edmTriggerResultsHelper_prescale_HLT_Ele27_WP80_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v1", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v10", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v11", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v12", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v13", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v14", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v15", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v16", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v17", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v18", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v19", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v2", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v20", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v3", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v4", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v5", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v6", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v7", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v8", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_IsoMu24_eta2p1_v9", edmTriggerResultsHelper_prescale_HLT_IsoMu24_eta2p1_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v1", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v10", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v11", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v12", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v13", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v14", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v15", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v16", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v17", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v18", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v19", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v2", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v20", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v3", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v4", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v5", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v6", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v7", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v8", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu13_Mu8_v9", edmTriggerResultsHelper_prescale_HLT_Mu13_Mu8_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v1", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v10", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v11", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v12", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v13", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v14", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v15", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v16", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v17", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v18", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v19", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v2", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v20", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v3", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v4", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v5", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v6", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v7", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v8", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescale_HLT_Mu17_Mu8_v9", edmTriggerResultsHelper_prescale_HLT_Mu17_Mu8_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v21", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v21);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v22", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v22);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v23", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v23);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v24", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v24);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v25", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v25);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9", edmTriggerResultsHelper_value_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v1", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v10", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v11", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v12", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v13", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v14", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v15", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v16", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v17", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v18", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v19", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v2", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v20", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v3", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v4", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v5", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v6", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v7", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v8", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Ele27_WP80_v9", edmTriggerResultsHelper_value_HLT_Ele27_WP80_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v1", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v10", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v11", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v12", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v13", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v14", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v15", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v16", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v17", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v18", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v19", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v2", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v20", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v3", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v4", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v5", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v6", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v7", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v8", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_IsoMu24_eta2p1_v9", edmTriggerResultsHelper_value_HLT_IsoMu24_eta2p1_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v1", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v10", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v11", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v12", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v13", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v14", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v15", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v16", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v17", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v18", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v19", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v2", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v20", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v3", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v4", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v5", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v6", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v7", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v8", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu13_Mu8_v9", edmTriggerResultsHelper_value_HLT_Mu13_Mu8_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v1", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v10", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v11", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v12", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v13", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v13);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v14", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v14);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v15", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v15);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v16", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v16);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v17", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v17);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v18", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v18);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v19", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v19);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v2", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v20", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v20);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v3", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v4", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v5", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v6", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v7", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v8", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.value_HLT_Mu17_Mu8_v9", edmTriggerResultsHelper_value_HLT_Mu17_Mu8_v9);
  stream.select("npatElectronHelper_selectedPatElectrons", nElectron);
  stream.select("npatElectronHelper_patElectronsLoosePFlow", nElectronPFlow);
  stream.select("nrecoGenParticle_genParticlesReduced", nGenParticle);
  stream.select("npatJet_selectedPatJetsPFlow", nJet);
  stream.select("npatMuon_selectedPatMuonsLoose", nMuon);
  stream.select("npatMuon_selectedPatMuonsLoosePFlow", nMuonPFlow);
  stream.select("nPileupSummaryInfo_addPileupInfo", nPileupSummaryInfo);
  stream.select("npatTau_selectedPatTaus", nTau);
  stream.select("nrecoTrackHelper_generalTracksReduced", nTrack);
  stream.select("nrecoVertex_offlinePrimaryVertices", nVertex);
  stream.select("sdouble_kt6PFJets_rho.value", sdoublePF_value);
  stream.select("sdouble_kt6CaloJets_rho.value", sdouble_value);

}
}; // end namespace evt
#endif
