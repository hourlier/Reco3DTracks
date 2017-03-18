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

//my includes
#include "DrawImage2D.h"
#include "Reco3D.h"


/*****************/
/*****************/
void Reco3D::initialize(){
    MaskedImg_v.resize(3);
    track_badchs_v.resize(3);
    tagged_img_v.resize(3);
    run_badchs_v.resize(3);
    full_vtx_img_v.resize(3);
    minstatus = larcv::chstatus::kGOOD;
    nplanes = 3;
    start_tick = 0;//track_meta.max_y();
    nticks = 120000;//track_meta.height();
    nchannels = 20000;
}
void Reco3D::initOneEvt(larlitecv::DataCoordinator &dataco){// starting from a larcv file
    ev_status = (larcv::EventChStatus*)(dataco.get_larcv_data(larcv::kProductChStatus, "tpc"));
    ev_Image2D = (larcv::EventImage2D*)(dataco.get_larcv_data(larcv::kProductImage2D, "tpc_hires_crop"));
    
    Run    = ev_status->run();
    SubRun = ev_status->subrun();
    Event  = ev_status->event();
    track_meta = ev_Image2D->Image2DArray()[0].meta();
    std::cout << Run << "\t" << SubRun << "\t" << Event << std::endl;

}
bool Reco3D::Are3planesOK(){
    bool Valid3planes = true;
    for(int iPlane = 0;iPlane<3;iPlane++){
        if(track_img_v.Image2DArray()[iPlane].meta().rows() == 0 || track_img_v.Image2DArray()[iPlane].meta().cols() == 0){
            std::cout << "ERROR NO 3 VALID PLANES => ABORT" << std::endl;
            Valid3planes = false;
        }
    }
    return Valid3planes;
}


void Reco3D::SetRunBadCh(){
    SetBadChParam();
    run_badchs_v = emptyalgo.makeBadChImage(BadChParams[1],BadChParams[1],BadChParams[2],BadChParams[3],BadChParams[4],BadChParams[5],BadChParams[6],*ev_status);
}

void Reco3D::SetBadChParam(){
    time_downsample_factor = (int)track_meta.pixel_height();
    wire_downsample_factor = (int)track_meta.pixel_width();

    BadChParams.clear();
    BadChParams.push_back(minstatus);//0
    BadChParams.push_back(nplanes);//1
    BadChParams.push_back(start_tick);//2
    BadChParams.push_back(nticks);//3    ==> here ultimately use the Geo package to ask for the maximum number of ticks and wires
    BadChParams.push_back(nchannels);//4
    BadChParams.push_back(time_downsample_factor);//5
    BadChParams.push_back(wire_downsample_factor);//6
}

void Reco3D::CropBadCh2TrackImg(){
    for(size_t iPlane = 0;iPlane<run_badchs_v.size();iPlane++){
        track_badchs_v.at(iPlane) = run_badchs_v.at(iPlane).crop(ev_Image2D->Image2DArray()[iPlane].meta().overlap(run_badchs_v.at(iPlane).meta()));
        full_vtx_img_v.at(iPlane) = ev_Image2D->Image2DArray()[iPlane].crop(track_badchs_v.at(iPlane).meta());
    }
}

void Reco3D::MaskDeadWires(){
    for(int iPlane = 0;iPlane<3;iPlane++){
        MaskedImg_v.at(iPlane) = MaskImage2D(track_img_v.Image2DArray()[iPlane],track_badchs_v.at(iPlane),full_vtx_img_v.at(iPlane));
        larcv::Image2D tagged_img(MaskedImg_v.at(iPlane).meta());
        tagged_img.paint(0);
        tagged_img_v.at(iPlane) = tagged_img;
    }
}

larcv::Image2D Reco3D::MaskImage2D(larcv::Image2D img_orig,larcv::Image2D img_bdch, larcv::Image2D img_true){
    larcv::Image2D img_masked(img_bdch.meta());
    img_masked.paint(0);
    for(int row = 0;row<img_masked.meta().rows();row++){
        for(int col = 0;col<img_masked.meta().cols();col++){
            if(img_bdch.pixel(row, col) == 0 && img_orig.pixel(row,col) !=0){
                img_masked.set_pixel(row,col,img_orig.pixel(row,col));
            }
            else{
                img_masked.set_pixel(row,col,0);
            }
        }
    }
    return img_masked;
}

void Reco3D::SetAstarEndPoints(int strt_rw, int end_rw, std::vector<int> strt_v, std::vector<int> end_v){
    start_row = strt_rw;
    end_row = end_rw;
    start_cols = strt_v;
    end_cols = end_v;
}

void Reco3D::RunAStarAlgo(){
    //_______________________________________________
    // Configure A* algo
    //-----------------------------------------------
    std::cout << "starting A*" << std::endl;
    larlitecv::AStar3DAlgoConfig config;
    config.accept_badch_nodes = true;
    config.astar_threshold.resize(3,1);
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

    //_______________________________________________
    // run A*
    //-----------------------------------------------
    RecoedPath = algo.findpath( MaskedImg_v, track_badchs_v, tagged_img_v, start_row, end_row, start_cols, end_cols, goal_reached );
    std::cout << "pathsize=" << RecoedPath.size() << std::endl;
}

/*****************/
/*****************/
void Reco3D::RecoEvent(){}
/*****************/
/*****************/
void Reco3D::finalize(){}
/*****************/
/*****************/
