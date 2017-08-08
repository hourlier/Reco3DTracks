#include "fstream"
#include "iostream"
#include <string>
#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TFile.h"

double LArDensity = 1.396; // g/cm^3
void MakeMuonSplines();
void MakeProtonSplines();


int main(){
    TFile *fout = TFile::Open("Proton_Muon_Range_dEdx_LAr_TSplines.root","RECREATE");
    std::cout << "  Muons" << std::endl;
    MakeMuonSplines();
    std::cout << "=========" << std::endl;
    std::cout << " Protons" << std::endl;
    MakeProtonSplines();
    fout->Close();
}


void MakeProtonSplines(){
    std::ifstream file("../data/ProtonArgonStoppingRangePSTAR.txt");
    std::string line;
    std::vector<double> lineInfo(7);

    std::vector< std::vector<double> > info;
    for(int i = 0;i<8;i++){
        getline(file,line);
        //std::cout << line << std::endl;
    }
    int iteration = 0;
    while(!file.eof() && iteration < 1000){
        iteration++;
        for(int item = 0;item<7;item++){
            file >> lineInfo[item];
        }
        if(lineInfo[0] <= 1e3){
            info.push_back(lineInfo);
        }
        if(file.eof()){iteration+=1e9;}
    }

    const int npx = info.size();
    double dEdx[npx],T[npx],Range[npx];
    for(int i = 0;i<info.size();i++){
        T[i]     = info[i][0];
        dEdx[i]  = info[i][3]/LArDensity;
        Range[i] = info[i][4]/LArDensity;
    }


    //TFile *fout = TFile::Open("ProtonTSplines.root","RECREATE");

    TSpline3 *sRange2T = new TSpline3("sProtonRange2T;CSDA Range [cm];T [MeV]",Range,T,npx);
    sRange2T->SetName("sProtonRange2T");
    sRange2T->Write();

    TSpline3 *sT2dEdx = new TSpline3("sProtonT2dEdx;T [MeV];dEdX [Mev/cm]",T,dEdx,npx);
    sT2dEdx->SetName("sProtonT2dEdx");
    sT2dEdx->Write();

    //fout->Close();
}
void MakeMuonSplines(){
    std::ifstream file("../data/MuonStoppingPower.txt");
    std::string line;
    std::vector<double> lineInfo(12);
    std::vector< std::vector<double> > info;
    for(int i = 0;i<10;i++){
        getline(file,line);
        //std::cout << line << std::endl;
    }
    int iteration = 0;
    while(!file.eof() && iteration < 1000){
        iteration++;
        for(int item = 0;item<12;item++){
            file >> lineInfo[item];
        }
        if(lineInfo[0] <= 1e3){
            info.push_back(lineInfo);
        }
        if(file.eof()){iteration+=1e9;}
    }

    const int npx = info.size();
    double dEdx[npx],T[npx],p[npx],Range[npx];
    for(int i = 0;i<info.size();i++){
        T[i]     = info[i][0];
        p[i]     = info[i][1];
        dEdx[i]  = info[i][7]/LArDensity;
        Range[i] = info[i][8]/LArDensity;
    }


    //TFile *fout = TFile::Open("MuonTSplines.root","RECREATE");

    TSpline3 *sRange2T = new TSpline3("sMuonRange2T;CSDA Range [cm];T [MeV]",Range,T,npx);
    sRange2T->SetName("sMuonRange2T");
    sRange2T->Write();

    TSpline3 *sP2dEdx = new TSpline3("sMuonT2dEdx;T [MeV];dEdX [Mev/cm]",T,dEdx,npx);
    sP2dEdx->SetName("sMuonT2dEdx");
    sP2dEdx->Write();
    
    //fout->Close();
}
