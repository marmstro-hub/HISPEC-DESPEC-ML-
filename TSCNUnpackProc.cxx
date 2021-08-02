
#include "TSCNUnpackProc.h"
#include "TGo4UserException.h"
#include "Riostream.h"

// ROOT Includes //
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TTree.h"
#include "TArc.h"

#include <time.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "Go4StatusBase/TGo4Picture.h"
#include "TGo4MbsEvent.h"
#include "TGo4Analysis.h"
#include "TGo4MbsSubEvent.h"
#include "TSCNUnpackEvent.h"
#include "QDC.h"
#include "TGo4Log.h"


#define RAW_DATA_LENGTH 400
#define RAW_BUFFER_LEN 2048

#define bPLASTIC_TAMEX_MODULES 9
#define bPLASTIC_TAMEX_CHANNELS 16

#include "Detector_System.cxx"
#include "PLASTIC_TWINPEAKS_Detector_System.h"

using namespace std;

//-----------------------------------------------------------
TSCNUnpackProc::TSCNUnpackProc():
  TGo4EventProcessor() {
    ffill = kTRUE;
    fshift = 6;
    load_bPlasticTamex_Allocationfile();
  }
//-----------------------------------------------------------
TSCNUnpackProc::TSCNUnpackProc(const char * name):
  TGo4EventProcessor(name) {
    TGo4Log::Info("TSCNUnpackProc: Create");

    // get_used_systems();

    Detector_Systems = new Detector_System * [7];
    Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_TWINPEAKS_Detector_System();

    RAW = new Raw_Event();
    load_bPlasticTamex_Allocationfile();

    //--- Histograms
    for (int i = 0; i < SCN_NUM_CHAN; i++) {
      fCr1Ch[i] = MakeTH1('D', Form("LYSO/QDC1Raw/QDC1channel%2d", i), Form("QDC1 channel %2d", i + 1), 5000, 0., 5000.);
      fCr2Ch[i] = MakeTH1('D', Form("LYSO/QDC2Raw/QDC2channel%2d", i), Form("QDC2 channel %2d", i + 1), 5000, 0., 5000.); ////2QDC SETUP
    }
  }
//-----------------------------------------------------------
TSCNUnpackProc::~TSCNUnpackProc() {
  TGo4Log::Info("TSCNUnpackProc: Delete");
}

//-----------------------------------------------------------
Bool_t TSCNUnpackProc::BuildEvent(TGo4EventElement * dest) {
  Bool_t isValid = kFALSE; // validity of output event

  TGo4MbsEvent * inp_evt = (TGo4MbsEvent * ) GetInputEvent(); // from this
  TSCNUnpackEvent * out_evt = (TSCNUnpackEvent * ) dest;

  if (inp_evt == 0) {
    cout << "SCNUnpackProc: no input event !" << endl;
    out_evt -> SetValid(isValid); // to store or not to store
    // default calling Fill method will set validity of out_evt to return value!
    return isValid;
  }
  isValid = kTRUE;
  #ifdef EXAMPLE_CODE
  if (inp_evt -> GetTrigger() > 11) {
    cout << "**** Tsis3316Proc: Skip trigger event" << endl;
    return kFALSE;
  }
  #endif
  //cout << "found an event!!!! " << inp_evt -> GetTrigger() << endl; //$$
  inp_evt -> ResetIterator();
  TGo4MbsSubEvent * psubevt(0);

  Int_t value = 0;
  Int_t chNo = 0;

  while ((psubevt = inp_evt -> NextSubEvent()) != 0) // subevent loop
  {
    Int_t * pdata = psubevt -> GetDataField();
    //Int_t lwords = psubevt -> GetIntLen();
    Int_t PrcID = psubevt -> GetProcid();
    //cout << "Proc ID " << PrcID << endl;

    switch (PrcID) {
    case 20: //It's the QDC. 
      //Int_t ModNo=floor(lwords/SCN_NUM_CHAN); //To know if I have more that one module
      Int_t EvCounter;
      Int_t geo_ad;
      Int_t Channels;
      //We only have 2 QDC, one later in the real experiment
      for (Int_t j = 0; j < 4; j++){
        Header * h = (Header * ) pdata;
        //cout << " header  " <<   h->check2 << endl;
        if (h -> check2 != 2) {
          //cout << "There's a problem reading the header" << endl;
          break;
        }
        geo_ad = h -> geo;
        // 	      cout<<"The geographcal addres is "<<geo_ad<<endl;
        Channels = h -> channels;
        //cout<<"1:The number of channels is "<<Channels<<endl;
        //cout<<"1 pdata "<<pdata<<"  "<<j<<endl;
        pdata++; //13march2020
        //     cout<<"2 pdata "<<pdata<<endl;
        // 	    cout << " !!!!!!!!!!  GEO "  << geo_ad << endl;
        for (Int_t i = 0; i < Channels; i++) {
          Data * d = (Data * ) pdata;

          if (d -> check != 0) {
            //cout<<"There's a problem reading the data"<<endl;
            break;
          }
          chNo = d -> channel;
          //cout<<"Reading channel "<<chNo<<endl;
          value = d -> value;
          if (geo_ad == 10) {
            fCr1Ch[chNo] -> Fill(value);
            out_evt -> fiCrate1[chNo] = value;
            //cout << "!!!!!!!!!!!!!!!!!!!!!!QDC 1 VAL : " << value << endl;
          } else if (geo_ad == 6) {
            fCr2Ch[chNo] -> Fill(value);
            //cout << "!!!!!!!!!!!!!!!!!!!!!!QDC 2 VAL : " << value << endl;
            out_evt -> fiCrate2[chNo] = value;
          }
          pdata++;
        }
        End * e = (End * ) pdata;

        if (e -> check != 4) {
          cout << "There was a problem reading the ending" << endl;
          //cout << "check was at "<<e->check<<endl ; 
          break;
        }

        EvCounter = e -> evCounter;
        // cout<<"Events: "<<EvCounter<<endl;
        //cout<<"5 pdata "<<pdata<<endl;
        pdata++;
        // cout<<"6 pdata "<<pdata<<endl;
      }
      //pdata++;// to skip the end 
      // cout<<"7 pdata "<<pdata<<endl;

      /*End of QDC case ------------------------------------*/
      break;

    case 80:
      Detector_Systems[2] -> Process_MBS(psubevt);
      Detector_Systems[2] -> Process_MBS(pdata);

      //get mbs stream data from unpacker (pointer copy solution)
      pdata = Detector_Systems[2] -> get_pdata();

      //get data from subevent and send to RAW
      Detector_Systems[2] -> get_Event_data(RAW);

      ///--------------------------------------------------------------------------------------------///
      /**Output bPlast Twin Peaks TAMEX **/
      ///--------------------------------------------------------------------------------------------///
      int bPlasfired[9];
      int Phys_Channel_Lead_Fast_bPlast[4][256];
      int Phys_Channel_Trail_Fast_bPlast[4][256];
      int Phys_Channel_Lead_Slow_bPlast[4][256];
      int Phys_Channel_Trail_Slow_bPlast[4][256];
      int bPlasdetnum_fast = -1;
      int bPlasdetnum_slow = -1;
      ///----------------------------------------------------------///
      /**Output bPlast Twin peaks TAMEX **/
      ///---------------------------------------------------------///
      //cout << "get_bPLAST_TWINPEAKS_tamex_hits" << RAW -> get_bPLAST_TWINPEAKS_tamex_hits() << endl;

      for (int i = 0; i < RAW -> get_bPLAST_TWINPEAKS_tamex_hits(); i++) { ///Loop over tamex ID's

        bPlasfired[i] = RAW -> get_bPLAST_TWINPEAKS_am_Fired(i);

        for (int j = 0; j < bPlasfired[i]; j++) { ///Loop over hits per board
          //cout<<"Input UNPACK RAW->get_bPLAST_TWINPEAKS_CH_ID(i,j) " <<RAW->get_bPLAST_TWINPEAKS_CH_ID(i,j) <<" RAW->get_bPLAST_TWINPEAKS_lead_T(i,j) " <<RAW->get_bPLAST_TWINPEAKS_lead_T(i,j) <<  " RAW->get_bPLAST_TWINPEAKS_trail_T(i,j) " <<RAW->get_bPLAST_TWINPEAKS_trail_T(i,j) << " i " << i << " j " << j <<  endl;

          ////NOW DEFINE FAST (ODD CHANNELS) AND SLOW  (EVEN)     
          if (j % 2 == 0) { //Lead even hits
            //cout<<"RAW->get_bPLAST_TWINPEAKS_CH_ID(i,j)  " <<RAW->get_bPLAST_TWINPEAKS_CH_ID(i,j)  << endl;
            ///Fast lead channels odd
            if (RAW -> get_bPLAST_TWINPEAKS_CH_ID(i, j) % 2 == 1) {
              if (RAW -> get_bPLAST_TWINPEAKS_lead_T(i, j) > 0) {
                Phys_Channel_Lead_Fast_bPlast[i][j] = TAMEX_bPlast_Chan[i][((RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) + 1) / 2) - 1];
                //cout<<"1 UNPACK Phys_Channel_Lead_Fast_bPlast[i][j] " <<Phys_Channel_Lead_Fast_bPlast[i][j] << " i " << i << " (RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1)/2 " << (RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1)/2 << " RAW->get_bPLAST_TWINPEAKS_lead_T(i,j) " <<RAW->get_bPLAST_TWINPEAKS_lead_T(i,j) << endl;
                //cout<<" Phys_Channel_Lead_Fast_bPlast[i][j] " <<  Phys_Channel_Lead_Fast_bPlast[i][j] << " (RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1)/2 " <<RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1/2 <<" i " << i << " j " << j << endl;
                bPlasdetnum_fast = TAMEX_bPlast_Det[i][((RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) + 1) / 2) - 1];
                out_evt -> fbPlasDetNum_Fast = bPlasdetnum_fast;
                //cout<<"bPlasdetnum_fast " <<bPlasdetnum_fast << " (RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j))/2 " <<(RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j))/2 << " i " << i << " RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) " <<RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) << endl;
                int chan_bPlast_fast_lead = Phys_Channel_Lead_Fast_bPlast[i][j];
                //cout<<"LEAD FAST bPlasdetnum_fast " <<bPlasdetnum_fast << " chan_bPlast_fast_lead " <<chan_bPlast_fast_lead << " i " << i << " j " << j << endl;
                out_evt -> fbPlas_FastChan[bPlasdetnum_fast] = chan_bPlast_fast_lead;
                //cout<<"Phys_Channel_Lead_Fast_bPlast[i][j] " <<Phys_Channel_Lead_Fast_bPlast[i][j] << " i " << i << " j " << j <<" RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) " <<RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) <<  endl;
                
		if (chan_bPlast_fast_lead > -1 && chan_bPlast_fast_lead < bPLASTIC_CHAN_PER_DET) {
                  //cout<<"chan_bPlast_fast_lead " << chan_bPlast_fast_lead <<" i " << i << " j " << j <<" (RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1)/2 " <<(RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1)/2<< endl;
                  int N1_fast = out_evt -> fbPlast_Fast_Lead_N[bPlasdetnum_fast][chan_bPlast_fast_lead]++;
                  out_evt -> fbPlast_Fast_Lead[bPlasdetnum_fast][chan_bPlast_fast_lead][N1_fast] = RAW -> get_bPLAST_TWINPEAKS_lead_T(i, j);
                  //cout<<"2 UNPACK EVENT " << event_number << " FAST bPlasdetnum_fast " << bPlasdetnum_fast << " chan_bPlast_fast_lead " << chan_bPlast_fast_lead << " N1_fast " <<N1_fast <<" out_evt->fbPlast_Fast_Lead[bPlasdetnum_fast][chan_bPlast_fast_lead][N1_fast] " <<out_evt->fbPlast_Fast_Lead[bPlasdetnum_fast][chan_bPlast_fast_lead][N1_fast] << " i " << i << " j " << j <<endl;       
                  //cout<<"FAST LEAD RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) " << RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) << " chan_bPlast_fast_lead " <<chan_bPlast_fast_lead << " N1_fast " <<N1_fast << " out_evt->fbPlast_Lead_Fast[chan_bPlast_fast_lead][N1_fast]  " <<out_evt->fbPlast_Lead_Fast[chan_bPlast_fast_lead][N1_fast]  << " i " << i << " j " << j << endl;
                }
              }
            }
            ///Slow lead channels, even 
            if (RAW -> get_bPLAST_TWINPEAKS_CH_ID(i, j) % 2 == 0) {
              Phys_Channel_Lead_Slow_bPlast[i][j] = TAMEX_bPlast_Chan[i][(RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) / 2) - 1];
              int chan_bPlast_slow_lead = Phys_Channel_Lead_Slow_bPlast[i][j];
              bPlasdetnum_slow = TAMEX_bPlast_Det[i][((RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) + 1) / 2) - 1];
              out_evt -> fbPlasDetNum_Slow = bPlasdetnum_slow;
              out_evt -> fbPlas_SlowChan[bPlasdetnum_slow] = chan_bPlast_slow_lead;
              //cout<<"LEAD SLOW bPlasdetnum_slow " <<bPlasdetnum_slow << " chan_bPlast_slow_lead " <<chan_bPlast_slow_lead << endl;

              if (chan_bPlast_slow_lead > -1 && chan_bPlast_slow_lead < bPLASTIC_CHAN_PER_DET) {
                int N1_slow = out_evt -> fbPlast_Slow_Lead_N[bPlasdetnum_fast][chan_bPlast_slow_lead]++;
                out_evt -> fbPlast_Slow_Lead[bPlasdetnum_fast][chan_bPlast_slow_lead][N1_slow] = RAW -> get_bPLAST_TWINPEAKS_lead_T(i, j);
                //cout<<"FAST bPlasdetnum_slow " << bPlasdetnum_slow << " chan_bPlast_fast_lead " << chan_bPlast_slow_lead << " N1_fast " <<N1_fast <<" out_evt->fbPlast_Lead_Fast[bPlasdetnum_fast][chan_bPlast_fast_lead][N1_fast] " <<out_evt->fbPlast_Lead_Fast[bPlasdetnum_fast][chan_bPlast_fast_lead][N1_fast] << " i " << i << " j " << j <<endl;
                //cout<<"SLOW LEAD RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) " << RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) << " chan_bPlast_slow_lead " <<chan_bPlast_slow_lead << " N1_slow " <<N1_slow << " out_evt->fbPlast_Lead_Slow[chan_bPlast_slow_lead][N1_slow]  " <<out_evt->fbPlast_Lead_Slow[chan_bPlast_slow_lead][N1_slow]  << " i " << i << " j " << j << endl;
              }
            }
          } ///End of lead hits

          if (j % 2 == 1) { //TRAIL 
            ///Fast trail channels even
            if (RAW -> get_bPLAST_TWINPEAKS_CH_ID(i, j) % 2 == 0 && (RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) + 1) / 2 < 256) {
              Phys_Channel_Trail_Fast_bPlast[i][j] = TAMEX_bPlast_Chan[i][(RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) / 2)];
              int chan_bPlast_fast_trail = Phys_Channel_Trail_Fast_bPlast[i][j];
              bPlasdetnum_fast = TAMEX_bPlast_Det[i][((RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) + 1) / 2) - 1];
              //cout<<"FAST CHAN " << Phys_Channel_Trail_Fast_bPlast[i][j]  << " bPlasdetnum_fast " <<bPlasdetnum_fast << endl;
              //cout<<"TRAIL FAST bPlasdetnum_fast " <<bPlasdetnum_fast << " chan_bPlast_fast_trail " <<chan_bPlast_fast_trail << endl;

              if (chan_bPlast_fast_trail > -1 && chan_bPlast_fast_trail < bPLASTIC_CHAN_PER_DET) {
                int N1_fast = out_evt -> fbPlast_Fast_Trail_N[bPlasdetnum_fast][chan_bPlast_fast_trail]++;
                out_evt -> fbPlast_Fast_Trail[bPlasdetnum_fast][chan_bPlast_fast_trail][N1_fast] = RAW -> get_bPLAST_TWINPEAKS_trail_T(i, j);
                //cout<<"FAST TRAIL RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) " << RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) << " chan_bPlast_fast_trail " <<chan_bPlast_fast_trail << " N1_fast " <<N1_fast << " out_evt->fbPlast_Trail_Fast[chan_bPlast_fast_trail][N1_fast]  " <<out_evt->fbPlast_Trail_Fast[chan_bPlast_fast_trail][N1_fast]  << " i " << i << " j " << j << endl;
              }
            }
            ///Slow trail channels even
            if (RAW -> get_bPLAST_TWINPEAKS_CH_ID(i, j) % 2 == 1 && RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) < 256) {
              Phys_Channel_Trail_Slow_bPlast[i][j] = TAMEX_bPlast_Chan[i][(RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) / 2) - 1];
              bPlasdetnum_slow = TAMEX_bPlast_Det[i][((RAW -> get_bPLAST_TWINPEAKS_physical_channel(i, j) + 1) / 2) - 1];
              //cout<<"Phys_Channel_Trail_Slow_bPlast[i][j] " <<Phys_Channel_Trail_Slow_bPlast[i][j] << " RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)+1/2 " <<RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j)/2-1 << " RAW " << RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) << endl; 
              //cout<<"SLOW CHAN " << Phys_Channel_Trail_Slow_bPlast[i][j]  << " bPlasdetnum_Slow " <<bPlasdetnum_slow << endl;
              int chan_bPlast_slow_trail = Phys_Channel_Trail_Slow_bPlast[i][j];

              if (chan_bPlast_slow_trail > -1 && chan_bPlast_slow_trail < bPLASTIC_CHAN_PER_DET) {
                int N1_slow = out_evt -> fbPlast_Slow_Trail_N[bPlasdetnum_slow][chan_bPlast_slow_trail]++;
                out_evt -> fbPlast_Slow_Trail[bPlasdetnum_slow][chan_bPlast_slow_trail][N1_slow] = RAW -> get_bPLAST_TWINPEAKS_trail_T(i, j);
                //cout<<"TRAIL SLOW bPlasdetnum_slow " <<bPlasdetnum_slow << " chan_bPlast_slow_trail " <<chan_bPlast_slow_trail << endl;
                //cout<<"SLOW TRAIL RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) " << RAW->get_bPLAST_TWINPEAKS_physical_channel(i, j) << " chan_bPlast_slow_trail " <<chan_bPlast_slow_trail << " N1_slow " <<N1_slow << " out_evt->fbPlast_Trail_Slow[chan_bPlast_slow_trail][N1_slow]  " <<out_evt->fbPlast_Trail_Slow[chan_bPlast_slow_trail][N1_slow]  << " i " << i << " j " << j << endl;

              }
            }
          }
        }
      }

    } //Switch
  } //while
  out_evt -> SetValid(isValid);
  return isValid;
}

//-----------------------------------------------------------
void TSCNUnpackProc::load_bPlasticTamex_Allocationfile() { //this works

  const char * format = "%d %d %d %d";
  ifstream data("Configuration_Files/bPlast_TAMEX_allocation.txt");
  if (data.fail()) {
    cerr << "Could not find bPlast_TAMEX_allocation file!" << endl;
    exit(0);
  }

  for (int i = 0; i < bPLASTIC_TAMEX_MODULES; i++) {
    for (int j = 0; j < bPLASTIC_TAMEX_CHANNELS; j++) {
      TAMEX_bPlast_Chan[i][j] = 0;
      TAMEX_bPlast_Det[i][j] = 0;
    }
  }
  int bPlastTamID = 0;
  int bPlastTamCh = 0;
  int bPlast_det = 0;
  int bPlast_ch = 0;
  string line;
  while (data.good()) {

    getline(data, line, '\n');
    if (line[0] == '#') continue;
    sscanf(line.c_str(), format, & bPlastTamID, & bPlastTamCh, & bPlast_det, & bPlast_ch);

    TAMEX_bPlast_Det[bPlastTamID][bPlastTamCh] = bPlast_det;
    TAMEX_bPlast_Chan[bPlastTamID][bPlastTamCh] = bPlast_ch;
    //cout << "TAMEX_bPlast_Det " << TAMEX_bPlast_Det[bPlastTamID] << " TAMEX_bPlast_Chan " << TAMEX_bPlast_Chan[bPlastTamID][bPlastTamCh] << " bPlast_det " << bPlast_det << " bPlastTamID " << bPlastTamID << " bPlastTamCh " << bPlastTamCh << " bPlast_ch" << bPlast_ch << endl;
  }
}
