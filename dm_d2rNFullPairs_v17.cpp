/**********************************************************************
Created by Alex Cherlin, 20/04/2015
**********************************************************************/

#include <iostream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TSystemDirectory.h>
#include <fstream>
#include <vector>
#include <TRandom3.h>
#include <TStopwatch.h>

using namespace std;

const Int_t TempChannelBase = 62000;
const Int_t Cath1ChannelBase = 60000;
const Int_t Cath2ChannelBase = 61000;

TString hname;
ofstream out_file12;
ofstream out_file3;

ifstream in_file;
string line;
TString str, fname_cfg, fname;
TString dir2do = "";

#define BYTES_PER_LINE 13

std::vector<Int_t> buffAM;
std::vector<Int_t> buffGM;
std::vector<Int_t> buffTime;
std::vector<Int_t> buffPixel;
std::vector<Int_t> buffData;
std::vector<Int_t> buffDataNeg;
std::vector<Int_t> buffTrigFlag;
std::vector<Int_t> buffTimeDetect;
std::vector<Int_t> buffTimeDetectNeg;

//
// Running parameters
//

///// Common path for input data, can be either of two:
///// 1. Directory containing multiple files belonging to one measurement.
///// 2. Directory containing any number of subdirectories each one of them containing multiple files belonging to one measurement.
/////    The outputs of all directories are not consolidated! Each directory is treated separately.	
/////
TString pathRoot = ".";
const TString setup_file = "./analysis_setup_d2rNFullPairs_v17.txt";

const TString data_file_suff = ".dat";

////////////////////////////////////////////////////////////////////////////////////////
// Running parameters, set here some default values.
// These parameters are re-set with values read from the analysis setup file.
////////////////////////////////////////////////////////////////////////////////////////

Int_t firstPixel = 1;
Int_t lastPixel = 968;
Int_t nPixXY = 22;
Int_t nAMs = 2;
Int_t GM2do = 1;

Bool_t printOutSetupFileParameters = kFALSE;
TString pathSetup = "";
TString fname_pixel_mapping_file = "";

Int_t prntOutF = 1000000;

Bool_t selectAllPairsOnly = kFALSE;
Bool_t selectPairsWithCathodesOnly = kFALSE;

Bool_t makeSortedEcalNTuple = kFALSE;
TString EcalibFileName = "";

Bool_t includeNegativeSignalData = kFALSE;
TString BL_NEG_FileName = "";

Bool_t saveSelfCheckASCIIoutputFile = kFALSE;

////////////////////////////////////////////////////////////////////////////////////////
// End of running parameters
////////////////////////////////////////////////////////////////////////////////////////

Long_t event = -1;
Long_t eventAll = 0;
Double_t *ADCperKEV;
Double_t *ADC_scale_offset;
Double_t ADCperKEV_cath[1100] = {0};
Double_t ADC_scale_offset_cath[1100] = {0};
TRandom3 *rand3;
Float_t *bl_neg;

Int_t getTemperatureChannelNumber(const Int_t, const Int_t);
Int_t getCath1Channel(const Int_t, const Int_t, const Int_t);
Int_t getCath2Channel(const Int_t, const Int_t, const Int_t);
Bool_t readAnalysisSetupFile(TString);
TStopwatch localTimer;
Float_t totalTimeElapced = 0;
TFile *f, *f2;

Int_t main()
{
	if (!readAnalysisSetupFile(setup_file))
	{
		cerr << "ERROR reading analysis setup file " << pathRoot + "/" + setup_file << ". Exiting." << endl;
		return 0;
	}
	if (selectPairsWithCathodesOnly && selectAllPairsOnly)
	{
		cerr << "ERROR: both selectAllPairsOnly and selectPairsWithCathodesOnly cannot be TRUE simultaneously. Exiting." << endl;
		return 0;
	}
	
	bl_neg = new Float_t[lastPixel+1];
	memset(bl_neg, 0, sizeof(bl_neg));
	
	rand3 = new TRandom3();
	rand3->SetSeed(0);
	
	while (pathSetup.BeginsWith(" ")) pathSetup.Remove(0,1);
	if (!((TString) pathSetup[pathSetup.Length()]).IsAlnum()) pathSetup = pathSetup.Chop();

	if (makeSortedEcalNTuple)
	{
		ADCperKEV = new Double_t[lastPixel+1];
		ADC_scale_offset = new Double_t[lastPixel+1];
		for (Int_t i = 0; i <= lastPixel; i++)
		{
			ADCperKEV[i] = 0;
			ADC_scale_offset[i] = 0;
		}		
		while (EcalibFileName.BeginsWith(" ")) EcalibFileName.Remove(0,1);
		if (!((TString) EcalibFileName[EcalibFileName.Length()]).IsAlnum()) EcalibFileName = EcalibFileName.Chop();
		TString fname2 = Form("%s/%s.txt",pathSetup.Data(),EcalibFileName.Data());
		if (gSystem->AccessPathName(fname2))  // Strange convention - this function return 0;s 1 (true) if path name doesn't exist !!!!
		{
			cerr << "ERROR: Energy calibration data file " << fname2 << " doesn't exist. Exiting." << endl;
			return 0;
		}
		in_file.open(fname2, ios::in);
		while (!in_file.eof())
		{
			Int_t pix = -1;
			Float_t a1 = -1, a2 = -1, boo = 0;
			in_file >> pix >> a2 >> a1;
			if (pix < 0) continue;
			if (pix <= lastPixel)
			{
				ADCperKEV[pix] = a2;
				ADC_scale_offset[pix] = a1;
			}
			if (pix >= Cath1ChannelBase)
			{
				ADCperKEV_cath[pix-Cath1ChannelBase] = a2;
				ADC_scale_offset_cath[pix-Cath1ChannelBase] = a1;
			}
		}
		in_file.close();
		in_file.clear();

		if (includeNegativeSignalData)
		{
			while (BL_NEG_FileName.BeginsWith(" ")) BL_NEG_FileName.Remove(0,1);
			if (!((TString) BL_NEG_FileName[BL_NEG_FileName.Length()]).IsAlnum()) BL_NEG_FileName = BL_NEG_FileName.Chop();
			fname2 = Form("%s/%s.txt",pathSetup.Data(),BL_NEG_FileName.Data());
			if (gSystem->AccessPathName(fname2))  // Strange convention - this function return 0;s 1 (true) if path name doesn't exist !!!!
			{
				cerr << "ERROR: BL NEG calibration data file " << fname2 << " doesn't exist. Exiting." << endl;
				return 0;
			}
			in_file.open(fname2, ios::in);
			while (!in_file.eof())
			{
				Int_t pix = -1;
				Float_t a1 = -1, boo = 0;
				in_file >> boo >> boo >> boo >> pix >> a1;
				if (pix < firstPixel || pix > lastPixel) continue;
				bl_neg[pix] = a1;
			}
			in_file.close();
			in_file.clear();
		}
	}
	Float_t *nEventsWithCathodeInTree = new Float_t[nAMs];
	Float_t *nEventsInTree = new Float_t[nAMs];
	Float_t nPairs = 0;
	Float_t nPairsWithCathodes = 0;
	for (Int_t im = 0; im < nAMs; im++)
	{	
		nEventsWithCathodeInTree[im] = 0;
		nEventsInTree[im] = 0;
	}
	
	localTimer.Start();
	const Char_t *ext = ".";
	TSystemDirectory dir("./","./");
	pathRoot = "./";
	TList *files = dir.GetListOfFiles();
	Int_t Nfiles = 0;
	TSystemFile *file;
	if (files)
	{
		TString fname, fnameRoot;
		TIter next(files);
		cout << "Counting files in " << dir2do << " ... ";
		
		const Char_t *ext = data_file_suff;
		while (file=(TSystemFile*)next())
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) Nfiles++;
			else files->Remove(file);
		}
		cout <<  Nfiles << " found." << endl;
		if (Nfiles == 0)
		{
			cerr << "ERROR: Current target directory:" << endl; 
			cerr << pathRoot + "/" + dir2do << endl;
			cerr << "doesn't contain valid data. Going to next directory, if there is any." << endl;
			return kFALSE;
		}
		next.Reset();

		Int_t toInit = kTRUE;
		Int_t toSkip = kFALSE;
		Int_t cnt11 = 0;
		Int_t nTSloops = 0;
		Int_t lastTimeStamp = -1;
		Int_t cnt = 0;
		Int_t cntNeg = 0;
		Char_t str1[10];
		Double_t lineNumber = 0;
		Int_t nFile = 0;
		Double_t *data;
		Double_t *data2;
		TTree **tree = new TTree*[nAMs];
		TTree *tree_ev;
		TTree **tree2 = new TTree*[nAMs];
		TTree *tree2_ev;
		Bool_t *eventHasCathode = new Bool_t[nAMs];
		Int_t *cathodeAmp = new Int_t[nAMs];
		Int_t *cathodeAmpNeg = new Int_t[nAMs];
		Float_t *Temp = new Float_t[nAMs];
		Int_t *nTrigPixels = new Int_t[nAMs];
		Int_t *cntAM = new Int_t[nAMs];
		Int_t *Nvar = new Int_t[nAMs];
		Int_t *Nvar2 = new Int_t[nAMs];
		for (Int_t im = 0; im < nAMs; im++)
		{
			eventHasCathode[im] = kFALSE;
			cathodeAmp[im] = 0;
			cathodeAmpNeg[im] = 0;
			Temp[im] = 0;
			nTrigPixels[im] = 0;
			cntAM[im] = 0;
			Nvar[im] = 0;
			Nvar2[im] = 0;
		}

		Long_t *event_d = new Long_t[nAMs];
		Int_t *nTrigPixels_d = new Int_t[nAMs];
		Int_t *AM_d = new Int_t[nAMs];
		Int_t *GM_d = new Int_t[nAMs];
		Int_t *nAMsInEvent_d = new Int_t[nAMs];
		Int_t *pixel_d = new Int_t[nAMs];
		Int_t *ADC_d = new Int_t[nAMs];
		Int_t *ADC_neg_d = new Int_t[nAMs];
		Int_t *time_d = new Int_t[nAMs];
		Int_t *triggerFlag_d = new Int_t[nAMs];
		Int_t *nloop_d = new Int_t[nAMs];
		Int_t *timeDetect_d = new Int_t[nAMs];
		Int_t *timeDetectPos_d = new Int_t[nAMs];
		Int_t *Temp_d = new Int_t[nAMs];
		Int_t *cathodeADC_d = new Int_t[nAMs];
		Int_t *cathodeADCneg_d = new Int_t[nAMs];
		Long_t event_ev;
		Int_t *AM_flag_ev = new Int_t[nAMs];
		Int_t *cath_flag_ev = new Int_t[nAMs];
		Int_t nAMsInEvent_ev;
		Int_t *nTrigPixels_ev = new Int_t[nAMs];
		
		Long_t *event_e = new Long_t[nAMs];
		Int_t *nTrigPixels_e = new Int_t[nAMs];
		Int_t *AM_e = new Int_t[nAMs];
		Int_t *GM_e = new Int_t[nAMs];
		Int_t *nAMsInEvent_e = new Int_t[nAMs];
		Int_t *pixel_e = new Int_t[nAMs];
		Float_t *E = new Float_t[nAMs];
		Float_t *E_neg = new Float_t[nAMs];
		Int_t *time_e = new Int_t[nAMs];
		Int_t *triggerFlag_e = new Int_t[nAMs];
		Int_t *nloop_e = new Int_t[nAMs];
		Int_t *timeDetect_e = new Int_t[nAMs];
		Int_t *timeDetectPos_e = new Int_t[nAMs];
		Float_t *Temp_e = new Float_t[nAMs];
		Int_t *pos_e = new Int_t[nAMs];
		Long_t event_eev;
		Int_t *AM_flag_eev = new Int_t[nAMs];
		Int_t *cath_flag_eev = new Int_t[nAMs];
		Int_t nAMsInEvent_eev;
		Int_t *nTrigPixels_eev = new Int_t[nAMs];

		Int_t cntFile = 0;
		while (file=(TSystemFile*)next())
		{
			fname = file->GetName();
			cout << fname << endl;
			if (toInit)
			{
				fnameRoot = fname;
				fnameRoot.ReplaceAll(data_file_suff,".root");
				TString rfileOut = pathRoot + "/" + dir2do + "/sortedData";
				if (!selectAllPairsOnly && !selectPairsWithCathodesOnly) rfileOut += "PairAnalysis";
				if (selectAllPairsOnly) rfileOut += "_allPairs";
				if (selectPairsWithCathodesOnly) rfileOut += "_allPairsWithCathodes";
				f = new TFile(rfileOut + ".root","recreate");
				if (makeSortedEcalNTuple) f2 = new TFile(rfileOut + "_Ecal.root","recreate");
				toInit = kFALSE;
				
				f->cd();
				tree_ev = new TTree("tree_events","Event tree of ADC data");
				tree_ev->Branch("event",&event_ev,"event/L");
				tree_ev->Branch("nAMsInEvent",&nAMsInEvent_ev,"nAMsInEvent/I");
				tree_ev->Branch("nTrigPixels",nTrigPixels_ev,Form("nTrigPixels[%d]/I",nAMs));
				tree_ev->Branch("AM_flag",AM_flag_ev,Form("AM_flag[%d]/I",nAMs));
				tree_ev->Branch("cath_flag",cath_flag_ev,Form("cath_flag[%d]/I",nAMs));
				
				for (Int_t im = 0; im < nAMs; im++)
				{
					tree[im] = new TTree(Form("tree_AM%d",im), Form("ADC data tree AM%d",im));
					tree[im]->Branch("event",&event_d[im],"event/L");
					tree[im]->Branch("nTrigPixels",&nTrigPixels_d[im],"nTrigPixels/I");
					tree[im]->Branch("AM",&AM_d[im],"AM/I");
					tree[im]->Branch("GM",&GM_d[im],"GM/I");
					tree[im]->Branch("nAMsInEvent",&nAMsInEvent_d[im],"nAMsInEvent/I");
					tree[im]->Branch("pixel",&pixel_d[im],"pixel/I");
					tree[im]->Branch("ADC",&ADC_d[im],"ADC/I");
					tree[im]->Branch("ADC_neg",&ADC_neg_d[im],"ADC_neg/I");
					tree[im]->Branch("time",&time_d[im],"time/I");
					tree[im]->Branch("cathodeADC",&cathodeADC_d[im],"cathodeADC/I");
					tree[im]->Branch("cathodeADCneg",&cathodeADCneg_d[im],"cathodeADCneg/I");
					tree[im]->Branch("triggerFlag",&triggerFlag_d[im],"triggerFlag/I");
					tree[im]->Branch("nloop",&nloop_d[im],"nloop/I");
					tree[im]->Branch("timeDetect",&timeDetect_d[im],"timeDetect/I");
					tree[im]->Branch("timeDetectPos",&timeDetectPos_d[im],"timeDetectPos/I");
					tree[im]->Branch("Temp",&Temp_d[im],"Temp/I");
					Nvar[im] = tree[im]->GetNbranches();
				}
				if (makeSortedEcalNTuple)
				{
					f2->cd();
					tree2_ev = new TTree("tree2_events","Event tree of E calib data");
					tree2_ev->Branch("event",&event_eev,"event/L");
					tree2_ev->Branch("nAMsInEvent",&nAMsInEvent_eev,"nAMsInEvent/I");
					tree2_ev->Branch("nTrigPixels",nTrigPixels_eev,Form("nTrigPixels[%d]/I",nAMs));
					tree2_ev->Branch("AM_flag",AM_flag_eev,Form("AM_flag[%d]/I",nAMs));
					tree2_ev->Branch("cath_flag",cath_flag_eev,Form("cath_flag[%d]/I",nAMs));
					
					for (Int_t im = 0; im < nAMs; im++)
					{
						tree2[im] = new TTree(Form("tree2_AM%d",im), Form("E calib data tree AM%d",im));
						tree2[im]->Branch("event",&event_e[im],"event/L");
						tree2[im]->Branch("nTrigPixels",&nTrigPixels_e[im],"nTrigPixels/I");
						tree2[im]->Branch("AM",&AM_e[im],"AM/I");
						tree2[im]->Branch("GM",&GM_e[im],"GM/I");
						tree2[im]->Branch("nAMsInEvent",&nAMsInEvent_e[im],"nAMsInEvent/I");
						tree2[im]->Branch("pixel",&pixel_e[im],"pixel/I");
						tree2[im]->Branch("E",&E[im],"E/F");
						tree2[im]->Branch("E_neg",&E_neg[im],"E_neg/F");
						tree2[im]->Branch("pos",&pos_e[im],"pos/I");
						tree2[im]->Branch("time",&time_e[im],"time/I");
						tree2[im]->Branch("triggerFlag",&triggerFlag_e[im],"triggerFlag/I");
						tree2[im]->Branch("nloop",&nloop_e[im],"nloop/I");
						tree2[im]->Branch("timeDetect",&timeDetect_e[im],"timeDetect/I");
						tree2[im]->Branch("timeDetectPos",&timeDetectPos_e[im],"timeDetectPos/I");
						tree2[im]->Branch("Temp",&Temp_e[im],"Temp/F");
						Nvar2[im] = tree2[im]->GetNbranches();
					}
				}
			}
			
			cout << "Processing dir " << dir2do << ", file " << fname << endl;
			
			ofstream out_file;
			TString fname2 = fname;
			fname2.ReplaceAll(data_file_suff,".output.txt");
			if (saveSelfCheckASCIIoutputFile) out_file.open(pathRoot+ "/" + dir2do + "/" + fname2, ios::out);		
			ofstream out_file2;
			TString fname22 = fname;
			fname22.ReplaceAll(data_file_suff,".Esorted.txt");
			if (saveSelfCheckASCIIoutputFile) out_file2.open(pathRoot+ "/" + dir2do + "/" + fname22, ios::out);		
			in_file.open(pathRoot+ "/" + dir2do + "/" + fname, ios::in | ios::binary);
			TH1F **baseline = new TH1F*[lastPixel+1];
			TH1F **baselineSel = new TH1F*[lastPixel+1];
			
			nTSloops = 0;
			lastTimeStamp = -1;
			cnt = 0;
			cntNeg = 0;
			for (Int_t im = 0; im < nAMs; im++)
			{
				eventHasCathode[im] = kFALSE;
				cathodeAmp[im] = 0;
				cathodeAmpNeg[im] = 0;
				Temp[im] = 0;
				nTrigPixels[im] = 0;
				nEventsWithCathodeInTree[im] = 0;
				cntAM[im] = 0;
			}
			
			buffAM.clear();
			buffGM.clear();
			buffTime.clear();
			buffPixel.clear();
			buffData.clear();
			buffDataNeg.clear();
			buffTrigFlag.clear();
			buffTimeDetect.clear();
			buffTimeDetectNeg.clear();			
			lineNumber = 0;
			cntFile++;
			while (!in_file.eof())
			{
				unsigned char amID=0;
				unsigned char gmID=0;
				unsigned short tstamp = 0;
				unsigned short pixelNum = 0;
				unsigned short energy = 0;
				unsigned char enPosEvent_flag=0;
				unsigned char threshold_flag = 0;
				unsigned short timeDetect = 0;
				unsigned char tdPosEvent_flag=0;
				char *buffer = new char [BYTES_PER_LINE];
				unsigned char Tl=0;
				unsigned char Th=0;
				
				in_file.read(buffer, BYTES_PER_LINE);
				lineNumber++;
				amID = (unsigned char)buffer[0];
				gmID = (unsigned char)buffer[1]; 
				tstamp = (unsigned char)buffer[2];
				tstamp = tstamp << 8;
				tstamp |= (unsigned char)buffer[3];
				pixelNum = (unsigned char)buffer[4];
				pixelNum = pixelNum << 8;
				pixelNum |= (unsigned char)buffer[5];
				energy = (unsigned char)buffer[6];
				energy = energy << 8;
				energy |= (unsigned char)buffer[7];
				enPosEvent_flag = (unsigned char)buffer[8];
				threshold_flag = (unsigned char)buffer[9];
				timeDetect = (unsigned char)buffer[10];
				timeDetect = timeDetect << 8;
				timeDetect |= (unsigned char)buffer[11];
				tdPosEvent_flag = (unsigned char)buffer[12];
				
				if (energy < 1) continue;
				if (GM2do != Int_t(gmID+0.5)) continue;
				if (pixelNum <= 0 || pixelNum >= 65000) continue;
				if (lastTimeStamp < 0) lastTimeStamp = tstamp;				
				if (pixelNum >= TempChannelBase)
				{
					if (pixelNum == getTemperatureChannelNumber(Int_t(amID+0.5),Int_t(gmID+0.5))) Temp[Int_t(amID+0.5)] = energy;
					continue;
				}
				if (tstamp != lastTimeStamp)
				{
					eventAll++;
					if (tstamp < lastTimeStamp) nTSloops++;
					lastTimeStamp = tstamp;
					Int_t isPair = 0;
					Int_t nCathodes = 0;
					for (Int_t im = 0; im < nAMs; im++)
					{
						if (eventHasCathode[im])
						{
							nAMsInEvent_d[im]++;
							nAMsInEvent_e[im]++;
							nEventsWithCathodeInTree[im]++;
							cath_flag_ev[im] = 1;
							cath_flag_eev[im] = 1;
							nCathodes++;
						}
						if (cntAM[im] > 0)
						{	
							nEventsInTree[im]++;
							isPair++;
							AM_flag_ev[im] = 1;
							AM_flag_eev[im] = 1;
							nAMsInEvent_ev++;
							nAMsInEvent_eev++;
						}
						nTrigPixels_ev[im] = nTrigPixels[im];
						nTrigPixels_eev[im] = nTrigPixels[im];
					}
					if (nCathodes == nAMs) nPairsWithCathodes++;
					if (isPair == nAMs) nPairs++;
					if (selectAllPairsOnly && isPair != nAMs) toSkip = kTRUE;
					if (selectPairsWithCathodesOnly && nCathodes != nAMs) toSkip = kTRUE;
					if (!toSkip)
					{
						event++;
						if (event%prntOutF == 0)
						{
							TString event_str = Form("%ld",event);
							for (Int_t is = event_str.Length(); is > 1; is--)
							{
								if ((is-event_str.Length()+2)%4 == 0) event_str.Insert(is-1,",",1);
							}
							cout << "Event " << event_str << endl;
						}
						f->cd();
						if (saveSelfCheckASCIIoutputFile) out_file << "event " << event << ":" << endl;
						for (Int_t im = 0; im < nAMs; im++)
						{
							if (saveSelfCheckASCIIoutputFile) out_file << "AM " << im << ":" << endl;
							for (Int_t m = 0; m < cnt; m++)
							{
								if (buffAM[m] != im) continue;
								if (buffData[m] < 1) continue;
								event_d[im] = event;
								nTrigPixels_d[im] = nTrigPixels[im];
								AM_d[im] = buffAM[m];
								GM_d[im] = buffGM[m];
								pixel_d[im] = buffPixel[m];
								ADC_d[im] = buffData[m];
								ADC_neg_d[im] = buffDataNeg[m];
								time_d[im] = buffTime[m];
								triggerFlag_d[im] = buffTrigFlag[m];
								nloop_d[im] = nTSloops;
								timeDetect_d[im] = buffTimeDetect[m];
								timeDetectPos_d[im] = buffTimeDetectNeg[m];
								Temp_d[im] = Temp[buffAM[m]];
								cathodeADC_d[im] = cathodeAmp[buffAM[m]];
								cathodeADCneg_d[im] = cathodeAmpNeg[buffAM[m]];
								tree[im]->Fill();
								if (saveSelfCheckASCIIoutputFile)
								{
									out_file << Form("%ld %d %d %d %d %d %d %d %d %d %d %d %d",event_d[im],nTrigPixels_d[im],AM_d[im],GM_d[im],pixel_d[im],ADC_d[im],ADC_neg_d[im],time_d[im],triggerFlag_d[im],nloop_d[im],timeDetect_d[im],timeDetectPos_d[im],Temp_d[im]) << endl;
								}
							}
						}
						out_file << "=====================================================" << endl;
						event_ev = event;
						tree_ev->Fill();
						
						if (makeSortedEcalNTuple)
						{
							f2->cd();
							if (saveSelfCheckASCIIoutputFile) out_file2 << "event " << event << ":" << endl;
							for (Int_t im = 0; im < nAMs; im++)
							{
								if (saveSelfCheckASCIIoutputFile) out_file2 << "AM " << im << ":" << endl;
								Int_t *idx = new Int_t[cnt];
								Float_t *arr = new Float_t[cnt];
								Int_t cnt_loc = 0;
								for (Int_t m = 0; m < cnt; m++)
								{
									if (buffAM[m] == im && buffPixel[m] <= lastPixel && ADCperKEV[buffPixel[m]] > 0.01) arr[m] = (buffData[m]+rand3->Uniform(1.)-ADC_scale_offset[buffPixel[m]])/ADCperKEV[buffPixel[m]];
									else arr[m] = -99999;
									cnt_loc++;
								}
								TMath::Sort(cnt,arr,idx,1);
								Int_t cntLoc = 0;
								for (Int_t m = 0; m < cnt; m++)
								{
									if (buffAM[idx[m]] != im) continue;
									if (buffData[idx[m]] < 1) continue;
									Float_t en = buffData[idx[m]];
									if (buffPixel[idx[m]] <= lastPixel) en = arr[idx[m]];
									if (buffPixel[idx[m]] <= lastPixel && en < -9000) continue;
									if (buffPixel[idx[m]] >= Cath1ChannelBase) en = (buffData[idx[m]]+rand3->Uniform(1.)-ADC_scale_offset_cath[buffPixel[idx[m]]-Cath1ChannelBase])/ADCperKEV_cath[buffPixel[idx[m]]-Cath1ChannelBase];
									Float_t enNEG = buffDataNeg[idx[m]];
									if (buffPixel[idx[m]] <= lastPixel) enNEG -= bl_neg[buffPixel[idx[m]]];
									event_e[im] = event;
									nTrigPixels_e[im] = nTrigPixels[im];
									AM_e[im] = buffAM[idx[m]];
									GM_e[im] = buffGM[idx[m]];
									pixel_e[im] = buffPixel[idx[m]];
									E[im] = en;
									E_neg[im] = enNEG;
									pos_e[im] = cntLoc;
									time_e[im] = buffTime[idx[m]];
									triggerFlag_e[im] = buffTrigFlag[idx[m]];
									nloop_e[im] = nTSloops;
									timeDetect_e[im] = buffTimeDetect[idx[m]];
									timeDetectPos_e[im] = buffTimeDetectNeg[idx[m]];
									Temp_e[im] = (Temp[buffAM[idx[m]]]*0.0005 - 1.525)/0.00567;
									tree2[im]->Fill();
									cntLoc++;
									if (saveSelfCheckASCIIoutputFile)
									{
										out_file2 << Form("%ld %d %d %d %d %.2f %.2f %d %d %d %d %d %d %.1f",event_e[im],nTrigPixels_e[im],AM_e[im],GM_e[im],pixel_e[im],E[im],E_neg[im],pos_e[im],time_e[im],triggerFlag_e[im],nloop_e[im],timeDetect_e[im],timeDetectPos_e[im],Temp_e[im]) << endl;
									}
								}
							}
							out_file2 << "=====================================================" << endl;
							event_eev = event;
							tree2_ev->Fill();
						}
					}
					buffAM.clear();
					buffGM.clear();
					buffTime.clear();
					buffPixel.clear();
					buffData.clear();
					buffDataNeg.clear();
					buffTrigFlag.clear();
					buffTimeDetect.clear();
					buffTimeDetectNeg.clear();					
					cnt = 0;
					cntNeg = 0;
					toSkip = kFALSE;
					nAMsInEvent_ev = 0;
					nAMsInEvent_eev = 0;
					for (Int_t im = 0; im < nAMs; im++)
					{
						eventHasCathode[im] = kFALSE;
						cathodeAmp[im] = 0;
						cathodeAmpNeg[im] = 0;
						Temp[im] = 0;
						nAMsInEvent_d[im] = 0;
						nTrigPixels[im] = 0;
						cntAM[im] = 0;
						AM_flag_ev[im] = 0;
						nTrigPixels_ev[im] = 0;
						AM_flag_eev[im] = 0;
						nTrigPixels_eev[im] = 0;
						cath_flag_ev[im] = 0;
						cath_flag_eev[im] = 0;
					}
				}
				if (Int_t(enPosEvent_flag+0.5) == 1)
				{
					buffAM.push_back(Int_t(amID+0.5));
					buffGM.push_back(Int_t(gmID+0.5));
					buffTime.push_back(tstamp);
					buffPixel.push_back(pixelNum);
					buffData.push_back(energy);
					buffTrigFlag.push_back(Int_t(threshold_flag+0.5));
					buffTimeDetect.push_back(timeDetect);
					if (pixelNum <= lastPixel && Int_t(threshold_flag+0.5) == 1) nTrigPixels[Int_t(amID+0.5)]++;
					buffDataNeg.push_back(0);
					buffTimeDetectNeg.push_back(0);
					cnt++;
					cntAM[Int_t(amID+0.5)]++;
				}
				if (pixelNum == getCath1Channel(Int_t(amID+0.5),Int_t(gmID+0.5),0))
				{
					eventHasCathode[Int_t(amID+0.5)] = kTRUE;
					cathodeAmp[Int_t(amID+0.5)] = energy;
				}
				if (Int_t(enPosEvent_flag+0.5) == 0)
				{
					for (Int_t p = 0; p < cnt; p++)
					{
						if (buffPixel[p] == pixelNum)
						{
							buffDataNeg[p] = energy;
							buffTimeDetectNeg[p] = timeDetect;
							cntNeg++;
							break;
						}
					}
				}
				delete buffer;
			}
			in_file.close();
			in_file.clear();
			
			if (saveSelfCheckASCIIoutputFile) out_file.close();
			if (saveSelfCheckASCIIoutputFile) out_file2.close();
			f->cd();
			tree_ev->Write();
			if (makeSortedEcalNTuple)
			{
				f2->cd();
				tree2_ev->Write();
			}
			for (Int_t im = 0; im < nAMs; im++)
			{
				f->cd();
				tree[im]->Write();
				if (makeSortedEcalNTuple)
				{
					f2->cd();
					tree2[im]->Write();
				}
			}
			cnt11++;
			nFile++;
			
			cout << Form("RealTime: %.2fs, CPUTime: %.2fs",localTimer.RealTime(), localTimer.CpuTime()) << endl;
			totalTimeElapced += localTimer.RealTime();
			localTimer.ResetRealTime();
			localTimer.ResetCpuTime();
			localTimer.Start();
		}
		f->Close();
		if (makeSortedEcalNTuple) f2->Close();
	}
	else
	{
		cerr << "ERROR: Current target directory:" << endl; 
		cerr << pathRoot + "/" + dir2do << endl;
		cerr << "doesn't contain valid data. Going to next directory, if there is any." << endl;
		return kFALSE;
	}
	
	if (totalTimeElapced < 60) cout << Form("Total running time is %.2f sec",totalTimeElapced) << endl;
	if (totalTimeElapced >= 60 && totalTimeElapced < 3600) cout << Form("Total running time is %dm:%.0f sec",Int_t(totalTimeElapced/60),totalTimeElapced - Int_t(totalTimeElapced/60)*60) << endl;
	if (totalTimeElapced >= 3600) cout << Form("Total running time is %dh:%dm:%.0fs",Int_t(totalTimeElapced/3600),Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60),
	totalTimeElapced - Int_t(totalTimeElapced/3600)*3600 - Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60)*60) << endl;

	cout << endl;
	out_file3.open(pathRoot+ "/" + dir2do + "/" + "logfile.txt", ios::out);
	TString num_str, num_str2;
	num_str = Form("%ld",event+1);
	for (Int_t is = num_str.Length(); is > 1; is--)
		if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
	cout << "Number of events in the tree = " << num_str << endl;
	out_file3 << "Number of events in the tree = " << num_str << endl;
	num_str = Form("%ld",eventAll);
	for (Int_t is = num_str.Length(); is > 1; is--)
		if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
	cout << "Total number of events = " << num_str << endl;
	out_file3 << "Total number of events = " << num_str << endl;
	for (Int_t im = 0; im < nAMs; im++)
	{
		num_str = Form("%.0f",nEventsInTree[im]);
		for (Int_t is = num_str.Length(); is > 1; is--)
			if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
		num_str2 = Form("%.0f",nEventsWithCathodeInTree[im]);
		for (Int_t is = num_str2.Length(); is > 1; is--)
			if ((is-num_str2.Length()+2)%4 == 0) num_str2.Insert(is-1,",",1);
		cout << Form("Number of events in the AM%d tree: %s, among them with cathode signal: %s or %.1f%%",im,num_str.Data(),num_str2.Data(),nEventsWithCathodeInTree[im]/nEventsInTree[im]*100) << endl;
		out_file3 << Form("Number of events in the AM%d tree: %s, among them with cathode signal: %s or %.1f%%",im,num_str.Data(),num_str2.Data(),nEventsWithCathodeInTree[im]/nEventsInTree[im]*100) << endl;
	}
	
	num_str = Form("%.0f",nPairs);
	for (Int_t is = num_str.Length(); is > 1; is--)
		if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
	cout << "Number of pairs: " << num_str << endl; 
	out_file3 << "Number of pairs: " << num_str << endl; 
	num_str = Form("%.0f",nPairsWithCathodes);
	for (Int_t is = num_str.Length(); is > 1; is--)
		if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
	cout << "Number of pairs with cathode signals: " << num_str << endl; 
	out_file3 << "Number of pairs with cathode signals: " << num_str << endl; 
	out_file3.close();
}

Int_t getTemperatureChannelNumber(const Int_t AM_number, const Int_t ASIC_number)
{
	return TempChannelBase + AM_number*4 + ASIC_number;
}

Int_t getCath1Channel(const Int_t AM_number, const Int_t ASIC_number, const Int_t serial)
{
	return Cath1ChannelBase + AM_number*16 + ASIC_number*4 + serial;
}

Int_t getCath2Channel(const Int_t AM_number, const Int_t ASIC_number, const Int_t serial)
{
	return Cath2ChannelBase + AM_number*16 + ASIC_number*4 + serial;
}

Bool_t readAnalysisSetupFile(TString fname)
{
	cout << "Reading analysis setup file:" << endl;
	cout << fname << endl;
	if (gSystem->AccessPathName(fname))  // Strange convention - this function return 0;s 1 (true) if path name doesn't exist !!!!
	{
		cerr << "ERROR: Analysis setup file " << fname << " doesn't exist. Exiting." << endl;
		return kFALSE;
	}
	
	Int_t buf = 0;
	in_file.open(fname, ios::in);
	while (!in_file.eof())
	{
		getline(in_file,line);
		TString sline = TString(line);
		if (sline.BeginsWith("/")) continue;
		
		if (sline.Contains("printOutSetupFileParameters ="))
		{
			sline.ReplaceAll("printOutSetupFileParameters =","");
			buf = sline.Atoi();
			printOutSetupFileParameters = kTRUE;
			if (buf == 0) printOutSetupFileParameters = kFALSE;
			if (printOutSetupFileParameters) cout << "printOutSetupFileParameters = " << printOutSetupFileParameters << endl;
		}

		if (sline.Contains("selectAllPairsOnly ="))
		{
			sline.ReplaceAll("selectAllPairsOnly =","");
			buf = sline.Atoi();
			selectAllPairsOnly = kTRUE;
			if (buf == 0) selectAllPairsOnly = kFALSE;
			if (printOutSetupFileParameters) cout << "selectAllPairsOnly = " << selectAllPairsOnly << endl;
		}
		
		if (sline.Contains("selectPairsWithCathodesOnly ="))
		{
			sline.ReplaceAll("selectPairsWithCathodesOnly =","");
			buf = sline.Atoi();
			selectPairsWithCathodesOnly = kTRUE;
			if (buf == 0) selectPairsWithCathodesOnly = kFALSE;
			if (printOutSetupFileParameters) cout << "selectPairsWithCathodesOnly = " << selectPairsWithCathodesOnly << endl;
		}
		
		if (sline.Contains("firstPixel ="))
		{
			sline.ReplaceAll("firstPixel =","");
			firstPixel = sline.Atoi();
			if (printOutSetupFileParameters) cout << "firstPixel = " << firstPixel << endl;
		}
		if (sline.Contains("lastPixel ="))
		{
			sline.ReplaceAll("lastPixel =","");
			lastPixel = sline.Atoi();
			if (printOutSetupFileParameters) cout << "lastPixel = " << lastPixel << endl;
		}
		if (sline.Contains("nPixXY ="))
		{
			sline.ReplaceAll("nPixXY =","");
			nPixXY = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nPixXY = " << nPixXY << endl;
		}
		
		if (sline.Contains("nAMs ="))
		{
			sline.ReplaceAll("nAMs =","");
			nAMs = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nAMs = " << nAMs << endl;
		}
		
		if (sline.Contains("GM2do ="))
		{
			sline.ReplaceAll("GM2do =","");
			GM2do = sline.Atoi();
			if (printOutSetupFileParameters) cout << "GM2do = " << GM2do << endl;
		}
		
		if (sline.Contains("pathSetup ="))
		{
			sline.ReplaceAll("pathSetup =","");
			pathSetup = sline;
			if (printOutSetupFileParameters) cout << "pathSetup = " << pathSetup << endl;
		}

		if (sline.Contains("fname_pixel_mapping_file ="))
		{
			sline.ReplaceAll("fname_pixel_mapping_file =","");
			fname_pixel_mapping_file = sline;
			if (printOutSetupFileParameters) cout << "fname_pixel_mapping_file = " << fname_pixel_mapping_file << endl;
		}
		
		if (sline.Contains("prntOutF ="))
		{
			sline.ReplaceAll("prntOutF =","");
			prntOutF = sline.Atoi();
			if (printOutSetupFileParameters) cout << "prntOutF = " << prntOutF << endl;
		}

		if (sline.Contains("makeSortedEcalNTuple ="))
		{
			sline.ReplaceAll("makeSortedEcalNTuple =","");
			buf = sline.Atoi();
			makeSortedEcalNTuple = kTRUE;
			if (buf == 0) makeSortedEcalNTuple = kFALSE;
			if (printOutSetupFileParameters) cout << "makeSortedEcalNTuple = " << makeSortedEcalNTuple << endl;
		}

		if (sline.Contains("EcalibFileName ="))
		{
			sline.ReplaceAll("EcalibFileName =","");
			EcalibFileName = sline;
			if (printOutSetupFileParameters) cout << "EcalibFileName = " << EcalibFileName << endl;
		}
		
		if (sline.Contains("includeNegativeSignalData ="))
		{
			sline.ReplaceAll("includeNegativeSignalData =","");
			buf = sline.Atoi();
			includeNegativeSignalData = kTRUE;
			if (buf == 0) includeNegativeSignalData = kFALSE;
			if (printOutSetupFileParameters) cout << "includeNegativeSignalData = " << includeNegativeSignalData << endl;
		}

		if (sline.Contains("BL_NEG_FileName ="))
		{
			sline.ReplaceAll("BL_NEG_FileName =","");
			BL_NEG_FileName = sline;
			if (printOutSetupFileParameters) cout << "BL_NEG_FileName = " << BL_NEG_FileName << endl;
		}
		
		if (sline.Contains("saveSelfCheckASCIIoutputFile ="))
		{
			sline.ReplaceAll("saveSelfCheckASCIIoutputFile =","");
			buf = sline.Atoi();
			saveSelfCheckASCIIoutputFile = kTRUE;
			if (buf == 0) saveSelfCheckASCIIoutputFile = kFALSE;
			if (printOutSetupFileParameters) cout << "saveSelfCheckASCIIoutputFile = " << saveSelfCheckASCIIoutputFile << endl;
		}
	}
	in_file.close();
	in_file.clear();
	return kTRUE;
}

