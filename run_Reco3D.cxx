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
#include "Reco3D.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TColor.h"

// my includes
#include "DrawImage2D.h"
#include "Reco3D.h"

int main( int nargs, char** argv ) {
    std::string codename = "\t \033[97m[run_3DTrackReco]\033[00m   ";
    std::cout << codename << std::endl;
    Reco3D *myReco = new Reco3D();
    myReco->SayHi();

    std::cout << codename <<"Set larlitecv" << std::endl;

    larlitecv::DataCoordinator dataco;
    dataco.add_inputfile("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/app/Cincinnati/data/numu_8000.root","larcv");
    dataco.configure("config.cfg","StorageManager","IOManager","3DRecoConfigurationFile");

    dataco.initialize();
    std::cout << dataco.get_nentries("larcv") << "entries in larcv input file" << std::endl;
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
    int Run, SubRun,Event, previousRun(-1);
    bool DrawPictures = true;
    bool Valid3planes = true;
    std::vector<std::vector<larcv::EventImage2D> > *particle_vv = 0;
    std::vector<std::vector<larcv::Vertex> > *particle_start_vv = 0;
    std::vector<std::vector<larcv::Vertex> > *particle_end_vv = 0;
    std::vector<std::vector<std::vector<larcv::Vertex> > > *particle_start2d_vvv = 0;
    std::vector<std::vector<std::vector<larcv::Vertex> > > *particle_end2d_vvv = 0;

    larcv::EventImage2D track_img_v;
    std::vector<larcv::Image2D> track_badchs_v(3);
    std::vector<larcv::Image2D> run_badchs_v;
    std::vector<larcv::Image2D> full_vtx_img_v(3);


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
        larcv::EventChStatus* ev_status = (larcv::EventChStatus*)(dataco.get_larcv_data(larcv::kProductChStatus, "tpc"));
        larcv::EventImage2D* ev_Image2D = (larcv::EventImage2D*)(dataco.get_larcv_data(larcv::kProductImage2D, "tpc_hires_crop"));

        Run    = ev_status->run();
        SubRun = ev_status->subrun();
        Event  = ev_status->event();
        if(!(run == Run && subrun == SubRun && event == Event)){
            std::cout << "ERROR, accessed wrong larcv entry" << std::endl;
            std::cout << Run << "\t" << SubRun << "\t" << Event << std::endl;
            continue;
        }

        for(int iVertex = 0; iVertex < particle_vv->size(); iVertex++){
            //_______________________________________________
            // define histograms and graphs
            //-----------------------------------------------
            if(DrawPictures){
                for(int iPlane = 0;iPlane<3;iPlane++){
                    histName = Form("hFullVertex_%d%d%d%d%d",Run,SubRun,Event,iVertex,iPlane);
                    hFullVertex[iPlane]=DrawImages(ev_Image2D->Image2DArray()[iPlane],histName);
                }
            }

            //_______________________________________________
            // Loop over the tracks
            //-----------------------------------------------

            for(int iTrack = 0; iTrack< particle_vv->at(iVertex).size();iTrack++){
                Valid3planes = true;
                track_img_v = particle_vv->at(iVertex)[iTrack];
                for(int iPlane = 0;iPlane<3;iPlane++){
                    if(track_img_v.Image2DArray()[iPlane].meta().rows() == 0 || track_img_v.Image2DArray()[iPlane].meta().cols() == 0){
                        std::cout << "ERROR NO 3 VALID PLANES => ABORT" << std::endl;
                        Valid3planes = false;
                    }
                }
                if(!Valid3planes){continue;}

                std::cout << "vertex #" << Form("%3d",iVertex) << "\t particle #" << iTrack << std::endl;
                std::cout << "\t\t\t starts at \t\t\tends  at" << std::endl;
                for(int iPlane = 0;iPlane <3;iPlane++){
                    std::cout << "\t\t\t(" << particle_start2d_vvv->at(iVertex)[iTrack][iPlane].X() << "," << particle_start2d_vvv->at(iVertex)[iTrack][iPlane].Y() <<")\t\t\t(" << particle_end2d_vvv->at(iVertex)[iTrack][iPlane].X() << "," << particle_end2d_vvv->at(iVertex)[iTrack][iPlane].Y() <<")" << std::endl;
                }

                //_______________________________________________
                // Get bad channels
                //-----------------------------------------------
                larlitecv::EmptyChannelAlgo emptyalgo;
                larcv::ImageMeta track_meta = ev_Image2D->Image2DArray()[0].meta();
                for(int iPlane = 0;iPlane<3;iPlane++){
                    while(gStart[iPlane]->GetN() > 0){gStart[iPlane]->RemovePoint(0);}
                    while(gEnd[iPlane]->GetN() > 0){gEnd[iPlane]->RemovePoint(0);}
                    gStart[iPlane]->SetPoint(0,particle_start2d_vvv->at(iVertex)[iTrack][iPlane].Y(),particle_start2d_vvv->at(iVertex)[iTrack][iPlane].X());
                    gEnd[iPlane]->SetPoint(0,particle_end2d_vvv->at(iVertex)[iTrack][iPlane].Y(),particle_end2d_vvv->at(iVertex)[iTrack][iPlane].X());
                }

                if(Run != previousRun){//new run, evaluate the new dead channels
                    int minstatus = larcv::chstatus::kGOOD; //=4
                    int nplanes = 3;
                    int start_tick = 0;//track_meta.max_y();
                    int nticks = 120000;//track_meta.height();
                    int nchannels = 20000;
                    int time_downsample_factor = (int)track_meta.pixel_height();
                    int wire_downsample_factor = (int)track_meta.pixel_width();

                    std::cout << "MakeBadChImage" << std::endl;
                    run_badchs_v = emptyalgo.makeBadChImage(minstatus,nplanes,start_tick, nticks, nchannels,time_downsample_factor,wire_downsample_factor,*ev_status);

                    previousRun = Run;
                }
                else{
                    std::cout << "same run, no need to evaluate the new dead channels" << std::endl;
                }

                /*for(int iPlane = 0;iPlane<3;iPlane++){
                    std::cout << "original image : " << ev_Image2D->Image2DArray()[iPlane].meta().dump() << std::endl;
                    std::cout << "bad ch   image : " << run_badchs_v.at(iPlane).meta().dump() << std::endl;
                }
                std::cout << std::endl;*/

                //_______________________________________________
                // Crop bad channels images
                //-----------------------------------------------
                std::cout << "Cropping bad ch images" << std::endl;
                for(size_t p = 0;p<run_badchs_v.size();p++){
                    track_badchs_v.at(p) = run_badchs_v.at(p).crop(ev_Image2D->Image2DArray()[p].meta().overlap(run_badchs_v.at(p).meta()));
                    full_vtx_img_v.at(p) = ev_Image2D->Image2DArray()[p].crop(track_badchs_v.at(p).meta());
                    //ev_Image2D->Image2DArray()[p] = ev_Image2D->Image2DArray()[p].crop(track_badchs_v.at(p).meta());
                }

                /*for(int iPlane = 0;iPlane<3;iPlane++){
                    std::cout << "original image : " << full_vtx_img_v.at(iPlane).meta().dump() << std::endl;
                    std::cout << "bad ch   image : " << track_badchs_v.at(iPlane).meta().dump() << std::endl;
                }
                std::cout << std::endl;
                std::cout << "original image : (w) " << full_vtx_img_v.at(0).meta().width() << "\t(h) " << full_vtx_img_v.at(0).meta().height() << std::endl;
                std::cout << "bad ch   image : (w) " << track_badchs_v.at(0).meta().width() << "\t(h) " << track_badchs_v.at(0).meta().height() << std::endl << std::endl;*/

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
                std::cout << "mask pixels corresponding to bad channels" << std::endl;
                std::vector<larcv::Image2D> MaskedImg_v(3);
                std::vector<larcv::Image2D> tagged_img_v(3);
                for(int iPlane = 0;iPlane<3;iPlane++){
                    MaskedImg_v.at(iPlane) = MaskImage2D(track_img_v.Image2DArray()[iPlane],track_badchs_v.at(iPlane),full_vtx_img_v.at(iPlane));
                    larcv::Image2D tagged_img(MaskedImg_v.at(iPlane).meta());
                    //std::cout << "painting tagged image with 0" << std::endl;
                    tagged_img.paint(0);
                    tagged_img_v.at(iPlane) = tagged_img;
                    //std::cout << "painting done and vector filled" << std::endl;
                }


                //_______________________________________________
                // Configure A* algo
                //-----------------------------------------------
                std::cout << "starting A*" << std::endl;
                larlitecv::AStar3DAlgoConfig config;
                config.accept_badch_nodes = true;
                config.astar_threshold.resize(3,0.001);
                //config.astar_threshold[2] = 10.0;
                config.astar_neighborhood.resize(3,10); //can jump over n empty pixels
                config.astar_start_padding = 5;// allowed region around the start point
                config.astar_end_padding = 5;  // allowed region around the end point
                config.lattice_padding = 10; // margin around the edges
                config.min_nplanes_w_hitpixel = 2;
                config.restrict_path = false; // do I want to restrict to a cylinder around the strainght line of radius defined bellow
                config.path_restriction_radius = 30.0;

                //_______________________________________________
                // Define A* algo and start/end points
                //-----------------------------------------------
                larlitecv::AStar3DAlgo algo( config );
                algo.setVerbose(0);
                int goal_reached(0);
                int start_row = particle_start2d_vvv->at(iVertex)[iTrack][0].X();
                int end_row = particle_end2d_vvv->at(iVertex)[iTrack][0].X();
                std::vector<int> start_cols(3);
                std::vector<int> end_cols(3);

                for(int iPlane = 0;iPlane<3;iPlane++){
                    start_cols[iPlane] = particle_start2d_vvv->at(iVertex)[iTrack][iPlane].Y();
                    end_cols[iPlane] = particle_end2d_vvv->at(iVertex)[iTrack][iPlane].Y();
                }


                //_______________________________________________
                // run A*
                //-----------------------------------------------

                std::vector<larlitecv::AStar3DNode> path = algo.findpath( MaskedImg_v, track_badchs_v, tagged_img_v, start_row, end_row, start_cols, end_cols, goal_reached );
                std::cout << "pathsize=" << path.size() << std::endl;

                //_______________________________________________
                // Draw outputs
                //-----------------------------------------------
                std::cout << std::endl;
                std::cout << "start_row  : " << start_row << std::endl;
                std::cout << "start_cols : " << start_cols[0] << "\t" << start_cols[1] << "\t" << start_cols[2] << std::endl;
                std::cout << std::endl;

                for(int iPlane = 0;iPlane<3;iPlane++){
                    while(gReco[iPlane]->GetN() > 0){gReco[iPlane]->RemovePoint(0);}
                }
                int nNodes[3] = {0,0,0};
                for(int iNode = 0;iNode<path.size();iNode++){
                    //std::cout << iNode << "\trow : " << path.at(iNode).row << "\t";
                    if(write3DPoint)myFlux << Run << " " << SubRun << " " << Event << " " << iVertex << " " << iTrack << " " << iNode << " " << path.at(iNode).u << " " << path.at(iNode).v << " " << path.at(iNode).w << std::endl;
                    for(int iPlane=0;iPlane<3;iPlane++){
                        gReco[iPlane]->SetPoint(nNodes[iPlane],path.at(iNode).cols[iPlane],path.at(iNode).row);
                        //std::cout << "\tcol " << iPlane << " : " << path.at(iNode).cols[iPlane] << "\t";
                        nNodes[iPlane]++;
                    }
                    //std::cout << std::endl;
                }
                std::cout << std::endl;
                //std::cout << "end_row    : " << end_row << std::endl;
                //std::cout << "end_cols   : " << end_cols[0] << "\t" << end_cols[1] << "\t" << end_cols[2] << std::endl;
                //std::cout << "Goal Reached? : " << goal_reached << std::endl;
                //std::cout << std::endl;
                if(DrawPictures){
                    for(int iPlane = 0;iPlane<3;iPlane++){
                        histName = Form("hMasked_%d%d%d%d%d%d",Run,SubRun,Event,iVertex,iTrack,iPlane);
                        hMasked[iPlane] = DrawImages(MaskedImg_v.at(iPlane),histName);
                        histName = Form("hBadChannels_%d%d%d%d%d%d",Run,SubRun,Event,iVertex,iTrack,iPlane);
                        hBadChannels[iPlane] = DrawImages(track_badchs_v.at(iPlane),histName);
                    }
                    cViews->Clear();
                    cViews->Divide(3,1);
                    for(int iPlane = 0;iPlane<3;iPlane++){
                        cViews->cd(iPlane+1);
                        hBadChannels[iPlane]->SetMarkerColor(kGray);
                        hFullVertex[iPlane]->Rebin2D();
                        hFullVertex[iPlane]->Draw("colz");
                        hBadChannels[iPlane]->Draw("same");
                        hMasked[iPlane]->SetMarkerStyle(1);
                        hMasked[iPlane]->SetMarkerSize(0.5);
                        hMasked[iPlane]->Draw("same colz");
                        gStart[iPlane]->SetMarkerStyle(20);
                        gStart[iPlane]->SetMarkerColor(2);
                        gEnd[iPlane]->SetMarkerStyle(20);
                        gEnd[iPlane]->SetMarkerColor(8);
                        gStart[iPlane]->Draw("P same");
                        gEnd[iPlane]->Draw("P same");
                        gReco[iPlane]->SetMarkerColor(4);
                        gReco[iPlane]->SetMarkerStyle(4);
                        gReco[iPlane]->SetMarkerSize(0.25);
                        gReco[iPlane]->SetLineColor(4);
                        gReco[iPlane]->SetLineStyle(2);
                        gReco[iPlane]->Draw("LP same");
                    }
                    cViews->Modified();
                    cViews->Update();
                    cViews->SaveAs(Form("views_%d_%d_%d_%d_%d.png",Run, SubRun,Event,iVertex,iTrack));
                }

                std::cout << std::endl;
            }
        }
    }
    if(write3DPoint)myFlux.close();
    return 0;
}

