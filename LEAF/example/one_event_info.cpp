#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"

#include "LEAF.hh"	
#include "HKManager.hh"	

#define ID_EVENT 1


void AnalyseEvent(WCSimRootEvent * tEvent, int iEventType);
void eventInfo(char * fname, int entryNum);

void AnalyseEvent(WCSimRootEvent * tEvent, int iEventType){
  WCSimRootTrigger * fRootTrigger;
  for(int iTrig = 0; iTrig < tEvent->GetNumberOfEvents(); iTrig++){
    fRootTrigger = tEvent->GetTrigger(iTrig);
		std::cout<<"Event Details for Debugging"<<std::endl;
		std::cout<<"Vtx info"<<std::endl;
		for(int i=0; i<3; i++){
			std::cout<<fRootTrigger->GetVtxs(0, i)<<std::endl;
		}
		std::cout<<"Dir info"<<std::endl;
		TObject * TrackDebug = (fRootTrigger->GetTracks())->At(0);
		WCSimRootTrack *wcTrackDebug = dynamic_cast<WCSimRootTrack*>(TrackDebug);
		for(int i=0; i<3; i++){
			std::cout<<wcTrackDebug->GetDir(i)<<std::endl;
		}
		std::cout<<"ParticleID: "<<wcTrackDebug->GetIpnu()<<std::endl;
		std::cout<<"Energy: "<<wcTrackDebug->GetE()<<std::endl;
		std::cout<<"Time: "<<wcTrackDebug->GetTime()<<std::endl;
		std::cout<<"Cherenkov Hits: "<<fRootTrigger->GetNcherenkovhittimes()<<std::endl;
		std::cout<<"N_Chr_digiHits: "<<fRootTrigger->GetNcherenkovdigihits()<<std::endl;
		std::cout<<"Trigger Type: "<<fRootTrigger->GetTriggerType()<<std::endl;
		std::cout<<"Trigger Info"<<std::endl;
		if(fRootTrigger->GetTriggerType() != -1){
			for(int i=0; i<3; i++){
				std::cout<<fRootTrigger->GetTriggerInfo()[i]<<std::endl;
			}
			std::cout<<"Event successfully worked"<<std::endl;
		}
	}
}

void eventInfo(char * fname, int entryNum){
  TFile * fInputFile = new TFile(fname, "READ");
  TTree * fInputTree = (TTree*) fInputFile->Get("wcsimT");
	WCSimRootEvent * fIDevent  = new WCSimRootEvent();

	fInputTree->SetBranchAddress("wcsimrootevent" , &fIDevent );
	fInputTree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
	
	int nPrimaryEvents = fInputTree->GetEntries();
	std::cout << "Event # = " << entryNum << " / " << nPrimaryEvents <<std::endl;
		
	fInputTree->GetEntry(entryNum);
	AnalyseEvent(fIDevent, ID_EVENT);
	std::cout<<"Informed successfully"<<std::endl;
}

int main(int argc, char **argv){
  // Optional arguments for infile, outfile: adapted from Analysis code for WCSim
  char * filename=NULL;
	int evt = 0; 
  int c = -1;

  while( (c = getopt(argc,argv,"f:e:a")) != -1 ){
    // Input in c the argument (-f etc...) and in optarg the next argument.
    // When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
    case 'f':
			std::cout<<optarg<<std::endl;
      filename = optarg;
      break;
		case 'e':
			if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
   			optarg = argv[optind++];
			}
			if (optarg == NULL) {
				std::cerr << "Error: Argument for -e option is missing." << std::endl;
				exit(EXIT_FAILURE);
			}
			evt = atoi(optarg);
      break;
    }
  }
	eventInfo(filename, evt);
  return 0;
}

