//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace larlitecv;
#pragma link C++ std::vector<std::vector<larcv::EventImage2D> >+;
#pragma link C++ namespace Reco3D+;
//ADD_NEW_CLASS ... do not change this line
#endif
