#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// config/storage: from LArCV
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h" // example of data product

// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Vertex.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/ImageMeta.h"

// larlitecv
#include "GapChs/GapChProcessor.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "ThruMu/AStar3DAlgo.h"
#include "Reco3DTracks/Reco3D.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TColor.h"

// my includes
#include "Reco3DTracks/DrawImage2D.h"

int main( int nargs, char** argv ) {
    std::string codename = "\t \033[97m[run_3DTrackReco]\033[00m   ";
    std::cout << codename << std::endl;
    Reco3D *myReco = new Reco3D();
    myReco->initialize();

    std::cout << codename <<"Set larlitecv" << std::endl;

    larlitecv::DataCoordinator dataco;
    dataco.add_inputfile("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/app/Cincinnati/data/numu_8000.root","larcv");
    dataco.configure("config.cfg","StorageManager","IOManager","3DRecoConfigurationFile");

    dataco.initialize();
    std::cout << dataco.get_nentries("larcv") << " entries in larcv input file" << std::endl;
    //return 0;

    std::cout << codename << "Read root file" << std::endl;

    //TFile *fIN = TFile::Open("test_adrien.root");
    //TFile *fIN = TFile::Open("aaa.root");
    //TFile *fIN = TFile::Open("test.root");
    TFile *fIN = TFile::Open("a.root");
    if(!fIN){std::cout << "ERROR, could no open input file" << std::endl;return 1;}
    std::ofstream myFlux;
    bool write3DPoint = false;
    if(write3DPoint){myFlux.open("Res.txt",std::ofstream::out | std::ofstream::trunc);}

    TTree *EventTree = (TTree*)fIN->Get("EventTree");

    UInt_t run,subrun,event,entry;
    int previousRun(-1);
    bool DrawPictures = true;
    std::vector<std::vector<larcv::EventImage2D> > *particle_vv = 0;
    std::vector<std::vector<larcv::Vertex> > *particle_start_vv = 0;
    std::vector<std::vector<larcv::Vertex> > *particle_end_vv = 0;
    std::vector<std::vector<std::vector<larcv::Vertex> > > *particle_start2d_vvv = 0;
    std::vector<std::vector<std::vector<larcv::Vertex> > > *particle_end2d_vvv = 0;

    larcv::EventImage2D Track_img_v;


    EventTree->SetBranchAddress("run",&run);
    EventTree->SetBranchAddress("subrun",&subrun);
    EventTree->SetBranchAddress("event",&event);
    EventTree->SetBranchAddress("entry",&entry);
    EventTree->SetBranchAddress("particle_vv",&particle_vv);
    EventTree->SetBranchAddress("particle_start_vv",&particle_start_vv);
    EventTree->SetBranchAddress("particle_end_vv",&particle_end_vv);
    EventTree->SetBranchAddress("particle_start2d_vvv",&particle_start2d_vvv);
    EventTree->SetBranchAddress("particle_end2d_vvv",&particle_end2d_vvv);

    TCanvas *cViews = new TCanvas("cViews","cViews",1200,400);
    std::string histName = "";
    TH2D *hFullVertex[3];
    //TH2D *hOneTrack[3];
    TH2D *hBadChannels[3];
    TH2D *hMasked[3];
    TGraph *gStart[3];
    TGraph *gEnd[3];
    TGraph *gReco[3];
    for(int iPlane = 0;iPlane<3;iPlane++){
        gStart[iPlane] = new TGraph();
        gEnd[iPlane]   = new TGraph();
        gReco[iPlane]  = new TGraph();
    }


    for(int i = 0;i<EventTree->GetEntries();i++){
        EventTree->GetEntry(i);
        std::cout << run << "\t" << subrun << "\t" << event << "\t" << entry << "\t"  << std::endl;

        dataco.goto_entry(entry,"larcv");
        myReco->initOneEvt(dataco);

        for(int iVertex = 0; iVertex < particle_vv->size(); iVertex++){
            //_______________________________________________
            // define histograms and graphs
            //-----------------------------------------------
            if(DrawPictures){
                for(int iPlane = 0;iPlane<3;iPlane++){
                    histName = Form("hFullVertex_%d%d%d%d%d",run,subrun,event,iVertex,iPlane);
                    hFullVertex[iPlane]=DrawImages(myReco->GetEvImg2D()->Image2DArray()[iPlane],histName);
                }
            }

            //_______________________________________________
            // Loop over the tracks
            //-----------------------------------------------

            for(int iTrack = 0; iTrack< particle_vv->at(iVertex).size();iTrack++){
                Track_img_v = particle_vv->at(iVertex)[iTrack]; // Get Image2D vector from Root file
                myReco->SetTrackImg(Track_img_v);
                //_______________________________________________
                // Check that I have all three image valid
                //-----------------------------------------------
                if(myReco->Are3planesOK() == false){continue;}
                std::cout << "vertex #" << Form("%3d",iVertex) << "\t particle #" << iTrack << std::endl;

                //_______________________________________________
                // Get TGraphs with start and end points
                //-----------------------------------------------

                for(int iPlane = 0;iPlane<3;iPlane++){
                    while(gStart[iPlane]->GetN() > 0){gStart[iPlane]->RemovePoint(0);}
                    while(gEnd[iPlane]->GetN() > 0){gEnd[iPlane]->RemovePoint(0);}
                }

                int start_row = particle_start2d_vvv->at(iVertex)[iTrack][0].X();
                int end_row = particle_end2d_vvv->at(iVertex)[iTrack][0].X();
                std::vector<int> start_cols(3);
                std::vector<int> end_cols(3);

                for(int iPlane = 0;iPlane<3;iPlane++){
                    start_cols[iPlane] = particle_start2d_vvv->at(iVertex)[iTrack][iPlane].Y();
                    end_cols[iPlane] = particle_end2d_vvv->at(iVertex)[iTrack][iPlane].Y();

                    gStart[iPlane]->SetPoint(0,start_cols[iPlane],start_row);
                    gEnd[iPlane]->SetPoint(0,end_cols[iPlane],end_row);
                }

                //_______________________________________________
                // If new run => get bad channels
                //-----------------------------------------------
                if(myReco->GetRun() != previousRun){//new run, evaluate the new dead channels
                    std::cout << "MakeBadChImage" << std::endl;
                    myReco->SetRunBadCh();
                    previousRun = myReco->GetRun();
                }
                else{
                    std::cout << "same run, no need to evaluate the new dead channels" << std::endl;
                }

                //_______________________________________________
                // Crop bad channels images
                //-----------------------------------------------
                myReco->CropBadCh2TrackImg();

                //_______________________________________________
                // Label gap channels
                //-----------------------------------------------
                /*int maxgap = 200;
                std::vector< larcv::Image2D> gapchimgs_v = emptyalgo.findMissingBadChs( track_img_v.Image2DArray(), track_badchs_v, 5, maxgap );
                std::cout << "original image : (w) " << ev_Image2D->Image2DArray()[0].meta().width() << "\t(h) " << ev_Image2D->Image2DArray()[0].meta().height() << std::endl;
                std::cout << "bad ch image   : (w) " << track_badchs_v.at(0).meta().width() << "\t(h) " << track_badchs_v.at(0).meta().height() << std::endl;
                std::cout << "label image    : (w) " << gapchimgs_v.at(0).meta().width()    << "\t(h) " << gapchimgs_v.at(0).meta().height() << std::endl << std::endl;
                //for ( size_t p=0; p<track_badchs_v.size(); p++ ) {
                //    larcv::Image2D& gapchimg = gapchimgs_v.at(p);
                //    gapchimg += track_badchs_v.at(p);
                //}
                std::cout << "label image += : (w) " << gapchimgs_v.at(0).meta().width()    << "\t(h) " << gapchimgs_v.at(0).meta().height() << std::endl << std::endl;*/

                //_______________________________________________
                // Create image with pixel values
                // at 0 for dead channels
                //-----------------------------------------------
                myReco->MaskDeadWires();

                //_______________________________________________
                // run A*
                //-----------------------------------------------

                myReco->SetAstarEndPoints(start_row, end_row, start_cols, end_cols);
                myReco->RunAStarAlgo();

                //_______________________________________________
                // Draw outputs
                //-----------------------------------------------

                for(int iPlane = 0;iPlane<3;iPlane++){
                    while(gReco[iPlane]->GetN() > 0){gReco[iPlane]->RemovePoint(0);}
                }

                int nNodes[3] = {0,0,0};
                for(int iNode = 0;iNode<myReco->GetRecoedPath().size();iNode++){
                    for(int iPlane=0;iPlane<3;iPlane++){
                        gReco[iPlane]->SetPoint(nNodes[iPlane],myReco->GetRecoedPath().at(iNode).cols[iPlane],myReco->GetRecoedPath().at(iNode).row);
                        nNodes[iPlane]++;
                    }
                }

                if(DrawPictures){
                    for(int iPlane = 0;iPlane<3;iPlane++){
                        histName = Form("hMasked_%d%d%d%d%d%d",myReco->GetRun(),myReco->GetSubRun(),myReco->GetEvent(),iVertex,iTrack,iPlane);
                        hMasked[iPlane] = DrawImages(myReco->GetMaskedImg().at(iPlane),histName);
                        histName = Form("hBadChannels_%d%d%d%d%d%d",myReco->GetRun(),myReco->GetSubRun(),myReco->GetEvent(),iVertex,iTrack,iPlane);
                        hBadChannels[iPlane] = DrawImages(myReco->GetTrckBdCh().at(iPlane),histName);
                    }

                    cViews->Clear();
                    cViews->Divide(3,1);
                    for(int iPlane = 0;iPlane<3;iPlane++){
                        cViews->cd(iPlane+1);

                        hFullVertex[iPlane]->Rebin2D();
                        hFullVertex[iPlane]->Draw("colz");

                        hBadChannels[iPlane]->SetMarkerColor(kGray);
                        hBadChannels[iPlane]->Draw("colz");

                        hMasked[iPlane]->SetMarkerStyle(1);
                        hMasked[iPlane]->SetMarkerSize(0.5);
                        hMasked[iPlane]->Draw("same colz");

                        gStart[iPlane]->SetMarkerStyle(20);
                        gStart[iPlane]->SetMarkerColor(2);
                        gStart[iPlane]->Draw("P same");

                        gEnd[iPlane]->SetMarkerStyle(20);
                        gEnd[iPlane]->SetMarkerColor(8);
                        gEnd[iPlane]->Draw("P same");

                        gReco[iPlane]->SetMarkerColor(4);
                        gReco[iPlane]->SetMarkerStyle(4);
                        gReco[iPlane]->SetMarkerSize(0.25);
                        gReco[iPlane]->SetLineColor(4);
                        gReco[iPlane]->SetLineStyle(1);
                        gReco[iPlane]->Draw("LP same");
                    }
                    cViews->Modified();
                    cViews->Update();
                    cViews->SaveAs(Form("views_%d_%d_%d_%d_%d.png",myReco->GetRun(),myReco->GetSubRun(),myReco->GetEvent(),iVertex,iTrack));
                }

                std::cout << std::endl;
            }
        }
    }
    if(write3DPoint)myFlux.close();
    return 0;
}

