//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "config.h"
#include "analyzer.h"
#include "functions.h"
#include "selection.h"
#include "hitInformation.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TStyle.h"
#include <time.h>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace evt;

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{

  // Get file list and histogram filename from command line
  clock_t t;
  t = clock();
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples
  bool isSignal = false;
  vector<string> filenames = getFilenames(cmdline.filelist);
  if(filenames[0].find("pMSSM") != std::string::npos || filenames[0].find("RECO_RAW2DIGI_L1Reco_RECO") != std::string::npos) isSignal = true;
  string outputfilename = cmdline.outputfilename;

  double TargetLifetime  = 0;
  double CurrentLifetime = 0;
  if(isSignal){
    TargetLifetime = getLifetime(outputfilename);
    cout<<endl<<"TargetLifetime = "<<TargetLifetime<<endl<<endl;
  }
  bool isData = false;
  if(filenames[0].find("MET_Run2012*_22Jan2013") != std::string::npos) isData = true;
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  //cout << "Number of events: " << nevents << endl;

  // Select variables to be read

  selectVariables(stream);


  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the
  // following line

  // TApplication app("analyzer", &argc, argv);

  /**
     Notes 1
     -------
     1. Use
     ofile = outputFile(cmdline.outputfile, stream)

     to skim events to output file in addition to writing out histograms.

     2. Use
     ofile.addEvent(event-weight)

     to specify that the current event is to be added to the output file.
     If omitted, the event-weight is defaulted to 1.

     3. Use
     ofile.count(cut-name, event-weight)

     to keep track, in the count histogram, of the number of events
     passing a given cut. If omitted, the event-weight is taken to be 1.
     If you want the counts in the count histogram to appear in a given
     order, specify the order, before entering the event loop, as in
     the example below

     ofile.count("NoCuts", 0)
     ofile.count("GoodEvent", 0)
     ofile.count("Vertex", 0)
     ofile.count("MET", 0)

     Notes 2
     -------
     By default all variables are saved. Before the event loop, you can use
  
     select(objectname)
	  
     e.g.,
	
     select("GenParticle")
  
     to declare that you intend to select objects of this type. The
     selection is done using

     select(objectname, index)
	  
     e.g.,
	  
     select("GenParticle", 3),
  
     which is called within the event loop. Call saveSelectedObjects()
     before a call to addEvent if you wish to save the selected objects.
     All other objects are saved by default.
	 
     NB: If you declare your intention to select objects of a given type
     by calling select(objectname), but subsequently fail to select
     them using select(objectname, index) then none will be saved!
  */
  

  outputFile ofile(cmdline.outputfilename);

  //-------------------------------------------------------------------------------------
  // Declare pileup histograms
  //-------------------------------------------------------------------------------------
  gStyle->SetOptStat(111111);

  TH1D *nVertices_0            = new TH1D("nVertices_0","nVertices_0",100,0,100);
  TH1D *hPU_NumInteractions_0  = new TH1D("hPU_NumInteractions_0","hPU_NumInteractions_0",100,0,100);
  TH1D *hTrueNumInteractions_0 = new TH1D("hTrueNumInteractions_0","hTrueNumInteractions_0",100,0,100);

  //-------------------------------------------------------------------------------------
  // Declaration of Variables
  //-------------------------------------------------------------------------------------
  cout<<endl<<endl<<"------------ Is it signal : "<<isSignal<<" --------------"<<endl<<endl;
  cout<<endl<<endl<<"------------ Is it data   : "<<isData<<" --------------"<<endl<<endl;
  //class Event noSelection("noSelection",ofile);
  //class Event triggerRequirements("triggerRequirements",ofile);
  //triggerRequirements.triggerRequirements = true;
  //class Event preselection("preselection",ofile);
  //preselection.triggerRequirements        = true;
  //preselection.trackPreselection          = true;
  //preselection.qcdSupression              = true;
  //class Event fullSelection("fullSelection",ofile);
  //fullSelection.triggerRequirements       = true;
  //fullSelection.trackPreselection         = true;
  //fullSelection.qcdSupression             = true;
  //fullSelection.trackCandidateCutFinal    = true;
  //fullSelection.TrackPtRequirement        = true;
  //fullSelection.NumOfLostOuterRequirement = true;
  //fullSelection.CaloIsolationRequirement  = true;
  //fullSelection.DeDxRequirement           = true;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Only Chi
  //class Event chiTracksnoSelection("chiTracksnoSelection",ofile);
  //if(isSignal) chiTracksnoSelection.onlyChi = true;

  //class Event noSelectionTightMuons("noSelectionTightMuons",ofile);
  //noSelectionTightMuons.tightMuonCut                =true;
  //class Event noSelectionTightElectrons("noSelectionTightElectrons",ofile);
  //noSelectionTightElectrons.tightElectronCut            =true;
  
  class Event chiTracksnoSelectionTightMuons("chiTracksnoSelectionTightMuons",ofile);
  if(isSignal) chiTracksnoSelectionTightMuons.onlyChi = true;
  chiTracksnoSelectionTightMuons.tightMuonCut         = true;

  class Event chiTracksnoSelectionTightElectrons("chiTracksnoSelectionTightElectrons",ofile);
  if(isSignal) chiTracksnoSelectionTightElectrons.onlyChi = true;
  chiTracksnoSelectionTightElectrons.tightElectronCut     = true;

  class Event chiTrackspreselectionTightMuons("chiTrackspreselectionTightMuons",ofile);
  if(isSignal) chiTrackspreselectionTightMuons.onlyChi = true;
  chiTrackspreselectionTightMuons.tightMuonCut         = true;
  chiTrackspreselectionTightMuons.trackPreselection    = true;

  class Event chiTrackspreselectionTightElectrons("chiTrackspreselectionTightElectrons",ofile);
  if(isSignal) chiTrackspreselectionTightElectrons.onlyChi = true;
  chiTrackspreselectionTightElectrons.tightElectronCut     = true;
  chiTrackspreselectionTightElectrons.trackPreselection    = true;

  class Event chiTrackspreselectionTagAndProbeElectron("chiTrackspreselectionTagAndProbeElectron",ofile);
  if(isSignal) chiTrackspreselectionTagAndProbeElectron.onlyChi    = true;
  chiTrackspreselectionTagAndProbeElectron.TagAndProbeElectronCut  = true;
  chiTrackspreselectionTagAndProbeElectron.trackPreselection       = true;

  class Event chiTrackspreselectionTagAndProbeMuon("chiTrackspreselectionTagAndProbeMuon",ofile);
  if(isSignal) chiTrackspreselectionTagAndProbeMuon.onlyChi    = true;
  chiTrackspreselectionTagAndProbeMuon.TagAndProbeMuonCut      = true;
  chiTrackspreselectionTagAndProbeMuon.trackPreselection       = true;

  class Event chiTrackspreselectionTagAndProbeTau("chiTrackspreselectionTagAndProbeTau",ofile);
  if(isSignal) chiTrackspreselectionTagAndProbeTau.onlyChi    = true;
  chiTrackspreselectionTagAndProbeTau.TagAndProbeTauCut       = true;
  chiTrackspreselectionTagAndProbeTau.trackPreselection       = true;

  class Event ReconstructionEfficiencyTagAndProbe("ReconstructionEfficiencyTagAndProbe",ofile);
  if(isSignal) ReconstructionEfficiencyTagAndProbe.onlyChi    = true;
  ReconstructionEfficiencyTagAndProbe.OneTightMuonCut         = true;


  //class Event chiTrackstriggerRequirements("chiTrackstriggerRequirements",ofile);
  //if(isSignal) chiTrackstriggerRequirements.onlyChi = true;
  //chiTrackstriggerRequirements.triggerRequirements  = true;

  //class Event chiTrackspreselection("chiTrackspreselection",ofile);
  //if(isSignal) chiTrackspreselection.onlyChi= true;
  //chiTrackspreselection.triggerRequirements = true;
  //chiTrackspreselection.trackPreselection   = true;
  //chiTrackspreselection.qcdSupression       = true;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // SM model control regions
  //class Event chiTracksSMControlCalo("chiTracksSMControlCalo",ofile);
  //if(isSignal) chiTracksSMControlCalo.onlyChi   = true;
  //chiTracksSMControlCalo.triggerRequirements    = true;
  //chiTracksSMControlCalo.trackPreselection      = true;
  //chiTracksSMControlCalo.qcdSupression          = true;
  //chiTracksSMControlCalo.trackCandidateCutFinal = true;
  //chiTracksSMControlCalo.TrackPtRequirement             = false;
  //chiTracksSMControlCalo.NumOfLostOuterRequirement      = false;
  //chiTracksSMControlCalo.CaloIsolationRequirement       = true;
  //chiTracksSMControlCalo.DeDxRequirement                = false;
  //chiTracksSMControlCalo.invertCaloIsolationRequirement = true;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //-------------------------------------------------------------------------------------
  // Read file for bad ecal cells and csc chambers 
  //-------------------------------------------------------------------------------------
  ifstream inputFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/BadCSCChambers.txt");
  int lines;
  int i=0;
  double n,m;
  while(inputFile>>n>>m){
    etaCSC.push_back(n);
    phiCSC.push_back(m);
    i++;
  }
  lines=i;
  inputFile.close();

  inputFile.open("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/DeadEcalChannelsNew.txt");
  i=0;
  while(inputFile>>n>>m){
    etaEcal.push_back(n);
    phiEcal.push_back(m);
    i++;
  }
  lines=i;
  
  //-------------------------------------------------------------------------------------
  // Read discriminator templates
  //-------------------------------------------------------------------------------------
  /*  if(isData){
    template_pixel = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Discrim_template_pixel_data_2012.root");
    template_strip = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Data7TeV_Deco_SiStripDeDxMip_3D_Rcd.root");
  }
  else{
    template_pixel = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Discrim_template_pixel_mc_2012.root");
    template_strip = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Discrim_Templates_MC_2012.root");
  }
  */
  //-------------------------------------------------------------------------------------
  // Declare additional branch addresses for hit information
  //-------------------------------------------------------------------------------------
  /*
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsDeDx",&HitsDeDx);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsPathlength",&HitsPathlength);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsShapetest",&HitsShapetest);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsSubdetId",&HitsSubdetid);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsEta",&HitsEta);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsPhi",&HitsPhi);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsTransverse",&HitsTransverse);
  */
  //-------------------------------------------------------------------------------------
  // Loop over events
  //-------------------------------------------------------------------------------------

  clock_t start, stop;
  double time = 0.0;

  assert((start = clock())!=-1);
  // ***********************************************************************************************************************
  // ----------  Stuff for PU reweighing -----------------------------------------------------------------------------------
  TFile file_PUdata("/nfs/dust/cms/user/rathjd/VBF-LS-tau/PU/DataPUFile_22Jan2013ReReco_Run2012.root", "read");
  TFile file_PUmc("/afs/desy.de/user/t/tlenz/HSCPworkdir/PUhistos/TrueNumInteractions_0.root", "read");
  //TFile file_PUmc("/nfs/dust/cms/user/rathjd/VBF-LS-tau/PU/S10MC_PUFile.root", "read");
 
  //TH1F *PUweights = (TH1F*)file_PU.Get("ratio");
  TH1F *PUweights = (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
  PUweights->Scale(1/PUweights->Integral());
  //cout<<"PUweights = "<<PUweights->Integral()<<endl;
  //  TH1F *PUmc = (TH1F*)file_PUmc.Get("analyzeHiMassTau/NVertices_0");
  TH1F *PUmc = (TH1F*)file_PUmc.Get("hTrueNumInteractions_0");
  PUmc->Scale(1/PUmc->Integral());
  //cout<<"PUmc = "<<PUmc->Integral()<<endl;
  
  PUweights->Divide(PUmc);
  //cout<<"PUweights = "<<PUweights->Integral()<<endl;
 
  weight = 1.;
  // -----------------------------------------------------------------------------------------------------------------------
  // ***********************************************************************************************************************

  cout<<endl<<"Number Of Events = "<<nevents<<endl<<endl;
  
  for(int entry=0; entry <nevents; ++entry)
    {

      // Read event into memory
      stream.read(entry);
      stream._chain->GetEntry(entry);

      // NB: call to clear object selection map (indexmap)
      initialize();
	  
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file analyzer.h to find out what structs
      // are available. Alternatively, you can call individual fill functions.
      fillObjects();

     
      /******************************************************************************************************************************
       ******************************************************************************************************************************
       ******************************************************************************************************************************
       *****************************************************************************************************************************/
      weight=1.;
      findChiInGenParticleCollection();
      //findChiInSimTrackCollection();
      //findChiDecayVertex();
      

      // Chargino event reweighting:
      
      if(isSignal && TargetLifetime !=0){
	CurrentLifetime = getLifetime(stream.filename());
      	//cout<<"CurrentLifetime = "<<CurrentLifetime<<endl;

	// 1.) Get both charginos from genParticle collection and their proper lifetime

	if (ChiTrack.size() > 2) cout << "Too many charginos!: " << ChiTrack.size() << endl;
	if (ChiTrack.size() < 2) cout << "Too few charginos!: " << ChiTrack.size() << endl;
	for(unsigned int i=0; i<ChiTrack.size();i++){

	 
	  double beta  = ChiTrack[i].genp/ChiTrack[i].genenergy;
	  double gamma = 1./std::sqrt(1.-pow(beta,2));
 
	  double rho      = sqrt(pow(ChiTrack[i].SimVertexposition_x - Vertex[0].x,2)+pow(ChiTrack[i].SimVertexposition_y - Vertex[0].y,2));
	  double z        = abs(ChiTrack[i].SimVertexposition_z - Vertex[0].z);
	  double distance = sqrt( pow(rho,2) + pow(z,2) );
	  double ProperLifetime = distance/(beta*gamma);

	  double wtTarget       = (1. / TargetLifetime)  * TMath::Exp(-(ProperLifetime) /  TargetLifetime  );  
	  double wtCurrent      = (1. / CurrentLifetime) * TMath::Exp(-(ProperLifetime) /  CurrentLifetime );
	  double wt = wtTarget / wtCurrent;

	  if(ProperLifetime<0){
	    cout<<"Warning:  Found event with ctau<0."<<endl;
	    wt=0.;
	  }

	  weight = weight*wt;
	}
      }
      /******************************************************************************************************************************
       ******************************************************************************************************************************
       ******************************************************************************************************************************
       *****************************************************************************************************************************/

      if(!edmEventHelper_isRealData){
	weight*=PUweights->GetBinContent(PUweights->FindBin(PileupSummaryInfo_getTrueNumInteractions[0]));
      }
      //weight =1.;
      ofile.count("NoCuts", weight);


      //------------- Calculate track dependent variables from hits and save them in trk coll ----
      /*
      for(unsigned int i=0; i<evt::Track.size();i++){
	
	
	double ASmiOnTheFly            = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,0); 
	double ASmiNPOnTheFly          = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,0); 
	double ASmiOnTheFly_3          = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,3);
	double ASmiNPOnTheFly_3        = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,0,3);
	double ASmiOnTheFly_7          = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,7);
	double ASmiNPOnTheFly_7        = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,0,7);
	double ASmiOnTheFly_woLastHit  = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,-1);

	Track[i].ASmi           = ASmiOnTheFly;
       	Track[i].ASmiNP         = ASmiNPOnTheFly;
	Track[i].ASmi_3         = ASmiOnTheFly_3;
       	Track[i].ASmiNP_3       = ASmiNPOnTheFly_3;
	Track[i].ASmi_7         = ASmiOnTheFly_7;
       	Track[i].ASmiNP_7       = ASmiNPOnTheFly_7;
	Track[i].ASmi_woLastHit = ASmiOnTheFly_woLastHit;

	double test           = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,0,&Track[i].DeDx1,&Track[i].DeDx2,&Track[i].DeDx3,&Track[i].DeDx4,&Track[i].Dx1,&Track[i].Dx2,&Track[i].Dx3,&Track[i].Dx4,&Track[i].MeasSize);

      }
      */
      
      //-------------------------------------------------------------- Cuts ---------------------------------------------------------
      
      nVertices_0->Fill(nVertex);
      hTrueNumInteractions_0->Fill(PileupSummaryInfo[0].getTrueNumInteractions);
      hPU_NumInteractions_0->Fill(PileupSummaryInfo[0].getPU_NumInteractions);


      //noSelection.Selection();
      //triggerRequirements.Selection();
      //preselection.Selection();
      //fullSelection.Selection();
      
      //chiTracksnoSelection.Selection();
      chiTracksnoSelectionTightMuons.Selection();
      chiTracksnoSelectionTightElectrons.Selection();
      chiTrackspreselectionTightMuons.Selection();
      chiTrackspreselectionTightElectrons.Selection();
      chiTrackspreselectionTagAndProbeElectron.Selection();
      chiTrackspreselectionTagAndProbeMuon.Selection();
      chiTrackspreselectionTagAndProbeTau.Selection();
      ReconstructionEfficiencyTagAndProbe.Selection();
      
      
    }//end of loop over events
 
  stop = clock();
  time = (double) (stop-start)/CLOCKS_PER_SEC;
  cout<<endl<<endl<<"time = "<<time/60.<<endl;

  stream.close();
  ofile.close();
  
  cout<<endl;
  cout<<"nChiInSimVertex = "<<nChiInSimVertex<<endl;
  cout<<"nChiInSimTrack  = "<<nChiInSimTrack<<endl<<endl;

  return 0;
}
