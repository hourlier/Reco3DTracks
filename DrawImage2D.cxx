#include "DrawImage2D.h"

TH2D* DrawImages(larcv::Image2D img, std::string histName){
    int Npx_x = img.meta().cols();
    int Npx_y = img.meta().rows();
    int start_x = 0;
    int start_y = 0;
    int end_x = start_x+Npx_x;
    int end_y = start_y+Npx_y;
    TH2D *himg = new TH2D(Form("%s",histName.c_str()),Form("%s;wire/col;time/row",histName.c_str()),Npx_x,start_x,end_x,Npx_y,start_y,end_y);
    for(int col = 0;col<Npx_x;col++){
        for(int row=0; row < Npx_y;row++){
            himg->SetBinContent(col,row,img.pixel(row,col));
        }
    }
    return himg;
}
