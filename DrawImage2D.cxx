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

larcv::Image2D MaskImage2D(larcv::Image2D img_orig,larcv::Image2D img_bdch, larcv::Image2D img_true){
    //std::cout << "entering MaskImage2D" << std::endl;
    larcv::Image2D img_masked(img_bdch.meta());
    img_masked.paint(0);
    //std::cout << "img_masked => " << img_masked.meta().dump() << std::endl;
    //std::cout << "img_orig   => " << img_orig.meta().dump() << std::endl;
    //std::cout << "img_bdch   => " << img_bdch.meta().dump() << std::endl;

    //std::cout << img_masked.meta().rows() << "rows and " << img_masked.meta().cols() << " cols" << std::endl;
    //std::cout << "img_masked declared" << std::endl;
    for(int row = 0;row<img_masked.meta().rows();row++){
        for(int col = 0;col<img_masked.meta().cols();col++){
            //std::cout << "MaskImage2D: " << row << "," << col << std::endl;
            if(img_bdch.pixel(row, col) == 0 && img_orig.pixel(row,col) !=0){
                //std::cout << "MaskImage2D: setting pixel " <<row<<","<<col << " with value " << img_orig.pixel(row,col)<< std::endl;
                img_masked.set_pixel(row,col,img_orig.pixel(row,col));
            }
            else{
                img_masked.set_pixel(row,col,0);
            }
        }
    }
    return img_masked;
}
