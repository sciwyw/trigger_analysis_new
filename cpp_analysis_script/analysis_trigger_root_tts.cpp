//
// Created by ineffablord on 23-11-28.
//
#include "fstream"
#include "iostream"
#include "sstream"
#include "vector"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "random"
#include "analysis_trigger_root_tts.h"
#include "set"
#include "cmath"

Double_t GetTotalTime(TTree* aintree){
    Long64_t hitnum;
    Double_t TotalTime;
    hitnum = aintree->GetEntries();
    TBranch* hittimebranch = aintree->GetBranch("hittime");
    hittimebranch->SetAddress(&TotalTime);
    aintree->GetEntry(hitnum - 1);
    return TotalTime;
}

void ApplyQE(TTree* aintree, const std::vector<QEdata>& qe, TTree* anewtree, std::mt19937 agenerator){
    //calculate wavelength
    Double_t energy;
    TBranch* energybranch = aintree->GetBranch("energy");
    energybranch->SetAddress(&energy);
    Double_t wavelength;
    TBranch* wavelengthbranch = aintree->Branch("wavelength", &wavelength, "wavelength/D");
    for (Long64_t i = 0; i < aintree->GetEntries(); ++i){
        aintree->GetEntry(i);
        wavelength = 1240 / energy;
        wavelengthbranch->Fill();
    }
    aintree->Write("",TObject::kOverwrite);

    anewtree->AddFriend(aintree);
    TBranch* eventid = aintree->GetBranch("eventid");
    Double_t eid;
    eventid->SetAddress(&eid);
    anewtree->Branch("eventid",&eid,"eventid/I");
    TBranch* hittime = aintree->GetBranch("hittime");
    Double_t hitt;
    hittime->SetAddress(&hitt);
    anewtree->Branch("hittime",&hitt,"hittime/D");
    TBranch* PMTid = aintree->GetBranch("PMTid");
    Double_t pmtid;
    PMTid->SetAddress(&pmtid);
    anewtree->Branch("PMTid",&pmtid,"PMTid/D");
    std::uniform_real_distribution<double> uniform(0,1);
    for (Long64_t i = 0; i < aintree->GetEntries(); ++i){
        aintree->GetEntry(i);
        Double_t QE = interpolate(qe,wavelength) / 100;
        Double_t probability = uniform(agenerator);
        if (QE >= probability){
            anewtree->Fill();
        }
    }
    anewtree->Write("",TObject::kOverwrite);
}

std::vector<QEdata> Readcsv(const std::string& csvfile) {
    std::vector<QEdata> qedata;
    std::ifstream file(csvfile);
    std::string line;
    //discard the first line
    if (std::getline(file, line)){}
    while (std::getline(file, line)){
        std::istringstream ss(line);
        std::string token;
        QEdata point{};
        std::getline(ss, token, ',');
        std::getline(ss, token, ',');
        point.wl = std::stod(token);
        std::getline(ss, token);
        point.qe = std::stod(token);
        qedata.push_back(point);
    }
    return qedata;
}

Double_t interpolate(const std::vector<QEdata>& qe, Double_t x){
    for (size_t i = 0; i < qe.size() - 1; i++){
        if (x >= qe[i].wl && x <= qe[i +1].wl){
            Double_t x0 = qe[i].wl;
            Double_t x1 = qe[i+1].wl;
            Double_t y0 = qe[i].qe;
            Double_t y1 = qe[i+1].qe;
            Double_t slope = (y1 - y0)/(x1 - x0);
            Double_t intery = y0 + slope*(x - x0);
            return intery;
        }
    }
    if (x <= qe.front().wl){
        return qe.front().qe;
    }
    else if (x >= qe.back().wl){
        return qe.back().qe;
    }
    else return 0;
}

void Addtts(const Double_t tts, TTree* anewtree, std::mt19937 agenerator) {
    std::normal_distribution<double> Normal(0,tts/2.3548);
    Double_t hitt;
    TBranch* hittime = anewtree->GetBranch("hittime");
    hittime->SetAddress(&hitt);
    Double_t atts;
    TBranch* addtts = anewtree->Branch("addtts", &atts, "addtts/D");
    for (Long64_t i = 0; i < hittime->GetEntries(); i++){
        hittime->GetEntry(i);
        atts = hitt + Normal(agenerator);
        addtts->Fill();
    }
    anewtree->Write("",TObject::kOverwrite);
}

void Trigger(TTree* thenewtree, TTree* outputtree, Int_t coin, Double_t window, std::ofstream* csvfile){
    Int_t id;
    TBranch* idbranch = thenewtree->GetBranch("eventid");
    idbranch->SetAddress(&id);
    Double_t time, timeflag;
    TBranch* timebranch = thenewtree->GetBranch("addtts");
    timebranch->SetAddress(&time);
    Double_t pmtid;
    TBranch* pmtidbranch = thenewtree->GetBranch("PMTid");
    pmtidbranch->SetAddress(&pmtid);
    Int_t oid;
    Double_t t;
    Bool_t tri = kFALSE;
    Double_t pmt;
    outputtree->Branch("eventid",&oid,"eventid/I");
    outputtree->Branch("time", &t, "time/D");
    outputtree->Branch("PMTid",&pmt,"PMTid/D");
    outputtree->Branch("state",&tri,"trigger_state/O");
    const Long64_t entries = thenewtree->GetEntries();
    Long64_t x = 0;
    for ( ; x < entries - coin + 1; x++){
        idbranch->GetEntry(x);
        pmtidbranch->GetEntry(x);
        timebranch->GetEntry(x);
        outputtree->GetEntry(x);
        oid = id;
        pmt = pmtid;
        t = timeflag = time;
        timebranch->GetEntry(x+coin-1);
        if (std::abs(time-timeflag)<=window){
            std::set<Double_t> pmts;
            Int_t number_of_pmts = 0;
            Long64_t I = 0;
            while (std::abs(time-timeflag)<=window and x+I < entries){
                if (pmts.count(pmt) == 0){
                    number_of_pmts++;
                }
                pmts.insert(pmt);
                I++;
                timebranch->GetEntry(x+I);
                pmtidbranch->GetEntry(x+I);
                pmt = pmtid;
            }
            if (number_of_pmts>=coin){
                outputtree->GetEntry(x);
                idbranch->GetEntry(x);
                timebranch->GetEntry(x);
                pmtidbranch->GetEntry(x);
                oid=id;t=time;pmt=pmtid;
                tri=kTRUE;
                outputtree->Fill();
                for (Long64_t i = 1; i < I; i++){
                    x++;
                    outputtree->GetEntry(x);
                    idbranch->GetEntry(x);
                    timebranch->GetEntry(x);
                    pmtidbranch->GetEntry(x);
                    oid=id;t=time;pmt=pmtid;
                    tri=kFALSE;
                    outputtree->Fill();
                }
            }else{
                outputtree->GetEntry(x);
                idbranch->GetEntry(x);
                timebranch->GetEntry(x);
                pmtidbranch->GetEntry(x);
                oid=id;t=time;pmt=pmtid;
                tri=kFALSE;
                outputtree->Fill();
            }
        }
        else{
            outputtree->GetEntry(x);
            tri = kFALSE;
            outputtree->Fill();
        }
    }
    for ( ; x < entries; x++){
        idbranch->GetEntry(x);
        timebranch->GetEntry(x);
        pmtidbranch->GetEntry(x);
        outputtree->GetEntry(x);
        oid = id;
        pmt = pmtid;
        t = time;
        tri = kFALSE;
        outputtree->Fill();
    }

    Long64_t trigger_count;
    trigger_count = outputtree->Draw("","state == 1","goff");
    *csvfile << int(trigger_count) << ",";
    outputtree->Write("",TObject::kOverwrite);
}

int main(int argc, char** argv){
    if (argc != 2) {
        std::cerr << "Wrong input arguments! The right way is: ./main [path/to/root file]";
        return -1;
    }

    //const Int_t totalevent = 20000000;

    const Char_t* ROOTPath = argv[1];
    const Double_t TTS = 1.3;
    std::random_device rd;
    std::mt19937 generator(rd());

    //open ROOT file. Tree:Hit0 Branch:eventid hittime PMTid energy x y z decaytime
    auto *infile= new TFile(ROOTPath, "UPDATE");
    TTree *intree;
    intree = dynamic_cast<TTree *>(infile->Get("Hit0"));
    TTree* newtree;
    newtree = new TTree("analysis", "fortrigger");
    newtree->SetDirectory(infile);
    //---Get total simulation time
    Double_t sim_time;
    sim_time = GetTotalTime(intree);
    std::cout << sim_time << std::endl;
    //---Get qe
    std::string qepath = "../../data/qe.csv";
    std::vector<QEdata> QE = Readcsv(qepath);
    //Apply qe in a new tree
    ApplyQE(intree,QE,newtree,generator);
    //Now the ROOT file should be: Tree:Hit0 Branch:eventid hittime PMTid energy x y z decaytime wavelength
    //                             Tree:analysis Branch:eventid hittime PMTid
    //---Add tts
    Addtts(TTS,newtree,generator);
    //Now the ROOT file should be: Tree:Hit0 Branch:eventid hittime PMTid energy x y z decaytime wavelength
    //                             Tree:analysis Branch:eventid hittime PMTid addtts
    //---Trigger
    std::ofstream outcsv;
    outcsv.open("trigger_num.csv", std::ios_base::app);
    const int coin_num[] = {1,2,3,4,5,6,7};
    const double time_win[] = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
    for (int i : coin_num){
        for (double j : time_win){
            std::string outpath = std::to_string(i) + "_" + std::to_string(int(j)) + "_" + "copy.root";
            auto* opath = const_cast<char *>(outpath.c_str());
            auto *outfile= new TFile(opath, "RECREATE");
            auto* trigger = new TTree("trigger", "batch_trigger");
            trigger->SetDirectory(outfile);
            Trigger(newtree,trigger,i,j,&outcsv);
            outfile->Write();
            outfile->Close();
            delete outfile;
        }
    }

    outcsv<<std::endl;
    outcsv.close();
    infile->Write();
    infile->Close();
    delete infile;
    return 0;
}