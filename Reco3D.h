#ifndef RECO3D_H
#define RECO3D_H
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

class Reco3D{
public:
    Reco3D(){}
    virtual ~Reco3D(){}
    void initialize();
    void initOneEvt(larlitecv::DataCoordinator &dataco);

    void SetTrackImg(larcv::EventImage2D img_v){track_img_v = img_v;}
    /*
     this function will need to be updated when reading hits from larlite files
     actually, track_img_v will need to be set directly from the data coordinator from larlite in initOneEvt
     */
    bool Are3planesOK();

    larcv::EventImage2D*        GetEvImg2D(){return ev_Image2D;}
    larcv::EventChStatus*       GetEvStat(){return ev_status;}
    larcv::ImageMeta            GetTrackMeta(){return track_meta;}
    larlitecv::EmptyChannelAlgo GetEmptyAlgo(){return emptyalgo;}
    std::vector<larcv::Image2D> GetRunBadCh(){return run_badchs_v;}
    std::vector<larcv::Image2D> GetFullVtxImg(){return full_vtx_img_v;}
    larcv::EventImage2D         GetTrackImg(){return track_img_v;}
    std::vector<larcv::Image2D> GetTrckBdCh(){return track_badchs_v;}
    std::vector<larcv::Image2D> GetMaskedImg(){return MaskedImg_v;}

    void SetTrackMeta(larcv::ImageMeta meta){track_meta = meta;}
    void SetRunBadCh();
    void SetMinStatus(int minstat){minstatus = minstat;}
    void SetBadChParam();
    void RecoEvent();
    void finalize();
    void CropBadCh2TrackImg();
    void MaskDeadWires();

    void SetAstarEndPoints(int strt_rw, int end_rw, std::vector<int> strt_v, std::vector<int> end_v);
    void RunAStarAlgo();

    int GetRun(){return Run;}
    int GetSubRun(){return SubRun;}
    int GetEvent(){return Event;}
    int GetMinStatus(){return minstatus;}
    int GetNplanes(){return nplanes;}
    std::vector<int> GetBadChParams(){return BadChParams;}

    larcv::Image2D MaskImage2D(larcv::Image2D img_orig,larcv::Image2D img_bdch, larcv::Image2D img_true);
    std::vector<larlitecv::AStar3DNode> GetRecoedPath(){return RecoedPath;}


private:
    larlitecv::EmptyChannelAlgo emptyalgo;
    std::vector<larcv::Image2D> MaskedImg_v;
    std::vector<larcv::Image2D> track_badchs_v;
    std::vector<larcv::Image2D> tagged_img_v;
    std::vector<larcv::Image2D> run_badchs_v;
    std::vector<larcv::Image2D> full_vtx_img_v;
    larcv::EventImage2D         track_img_v;
    larcv::EventChStatus*       ev_status;
    larcv::EventImage2D*        ev_Image2D;
    larcv::ImageMeta            track_meta;
    std::vector<int>          BadChParams;
    int                       Run;
    int                       SubRun;
    int                       Event;
    int                       minstatus;
    int                       nplanes;
    int                       start_tick;
    int                       nticks;
    int                       nchannels;
    int                       time_downsample_factor;
    int                       wire_downsample_factor;

    int start_row;
    int end_row;
    std::vector<int> start_cols;
    std::vector<int> end_cols;

    std::vector<larlitecv::AStar3DNode> RecoedPath;
};

#endif
