#include "Base/DataFormatConstants.h"
// larcv data
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/Vertex.h"
#include "DataFormat/EventChStatus.h"
#include "DataFormat/ImageMeta.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "Base/DataCoordinator.h"

// larlitecv
#include "GapChs/GapChProcessor.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "ThruMu/AStar3DAlgo.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"


TH2D* DrawImages(larcv::Image2D img, std::string histName);
