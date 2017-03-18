#ifndef RECO3D_H
#define RECO3D_H
//larlite
#include "Base/DataFormatConstants.h"
#include "DataFormat/track.h"

//larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Vertex.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/ImageMeta.h"

//larlitecv
#include "GapChs/GapChProcessor.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "ThruMu/AStar3DAlgo.h"

class Reco3D{
public:
    Reco3D(){}
    virtual ~Reco3D(){}
    void initialize();
    void RecoEvent();
    void finalize();
    void SayHi();
private:
    std::vector<larcv::Image2D> MaskedImg_v;
};

#endif
