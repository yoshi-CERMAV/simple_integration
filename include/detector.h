//
//  detector.h
//  
//
//  Created by Yoshiharu Nishiyama on 10/12/2023.
//

#define NOFILTER
#ifndef detector_h
#define detector_h

#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>
#include <iostream>
#include <fstream>
#include <memory>
//#ifndef NOFILTER
//#include "Savitzky_Golay_2d.h"
//#endif
//#include "check_mask.h"


using namespace H5;
using namespace std;
class apply_poni;
typedef std::shared_ptr<H5::H5File> H5FilePtr;

class detector
{
    void copy_data(double *dat, int *idat){
        for(int i = 0; i != image_size; i++){
            dat[i] = idat[i];//*flat[i];
        }
    }
    template<class T>
    void flat_data(T *dat){
        for(int i = 0; i != image_size; i++){
            dat[i] *= response[i];
        }
    }
    template<class T1, class T2>
    void flat_data(T1 *dat, T2 *out_dat){
        for(int i = 0; i != image_size; i++){
            out_dat = dat[i] * response[i];
        }
    }

    int cook(double *dat, double *clean_dat){
        flat_data(dat);
  //      cout << "flattened"<<endl;
  //      if(!check) return 1;
  //      check->get_from_neighbors(dat);
  //      cout <<"filled"<<endl;
#ifndef NOFILTER
        filter->apply(dat, clean_dat, image_size);
#endif
  //      cout << "filtered"<<endl;
        return 0;
    }
public:
    detector(){
        mask = NULL;
        response=NULL;
   //     check = NULL;
    }
    void prepare(int *data){
   //     cout <<"preparing"<<endl;
        copy_data(temp, data);
//        cout << "copied"<<endl;
        cook(temp, clean_data);
    }
    void copy_to_bg()
    {
        for(int i = 0; i < image_size; i++){
            clean_bg[i] = clean_data[i];
        }
    }
    
    void init(const char H5[]){
        H5File* file = 0;
        try {file = new H5File(H5, H5F_ACC_RDONLY);
        } catch(const FileIException& f){
            f.printErrorStack();
        }
        DataSet dataset = file->openDataSet("/entry_0000/pyFAI/Detector/pixel_corners");
        DataSpace dataspace = dataset.getSpace();
        
        int rank = dataspace.getSimpleExtentNdims() ;
        hsize_t dims_out[rank];
        const int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
        pixel_corners_dim = 1;
        for(int i = 0; i < rank; i++) {
            cout << dims_out[i]<<endl;
            pixel_corners_dim *= dims_out[i];
        }
        dim0 = dims_out[1];
        dim1 = dims_out[0];
        image_size = dims_out[0] * dims_out[1];
        pixel_corners = new float [pixel_corners_dim];
        center_pos = new float[image_size*3];
        candidates = new int[image_size];
        candidates1 = new int[image_size];

        mask = new char [image_size];
        response = new float[image_size];
        cout<< "no filter"<<endl;
#ifndef NOFILTER
        filter = new Filter2d (9, 9, 0, 2, dim0, dim1);
#endif
        filtered = new double[image_size];
        temp = new double[image_size];
        clean_data = new double[image_size];
        clean_bg = new double[image_size];
        clean_subtracted = new double[image_size];
        
        
        for(int i = 0; i < image_size; i++){
            candidates[i] = i;
            mask[i]= 0;
            response[i] = 1;
        }

        candidate_size = image_size;
        dataset.read(pixel_corners, PredType::IEEE_F32LE, dataspace, dataspace);
        
        float *ptr1 = center_pos;
        float *ptr = pixel_corners;
        for(int i = 0; i < image_size; i++, ptr1+=3, ptr +=12){
            for(int j = 0; j < 3; j++){
                ptr1[j] = ptr[j] + ptr[j+3] + ptr[j+6] + ptr[j+9];
                ptr1[j] *= 0.25;
            }
        }
        get_min_max(min_x, max_x ,min_y, max_y);
        float shape = (max_y-min_y)/(max_x-min_x);
    }
    detector(const char H5[])
    {
        mask = NULL;
        response=NULL;
    //    check = NULL;
        init(H5);
    }
    
    double sqr(double x){return x*x;}

    void add_mask(double x0, double y0, double r){
        double center_x = x0+r;
        double center_y = y0+r;
        int ix = x0;
        int iy = y0;
        for(int j = 0; j < 2*r + 2; j++){
            for(int i = 0; i < 2*r + 2; i++){
                if (sqr(i-r)+sqr(j-r)< r*r) mask[ix + i + 578* (iy+j)] = 1;
            }
        }
    }
    
    void load_mask(const char file[])
    {
        ifstream fi(file);
        fi.read(mask, image_size);
        int count = 0;
        for(int i = 0; i < image_size; i++){
            if(! mask[i]){
                candidates[i] = count;
                count++;
            }
        }
//        if(! check) check = new check_mask(mask, dim0, dim1);
//        else check->reset(mask);
        candidate_size =count;
        cout << "count = "<< count<<endl;
        get_positions();
        masked = true;
    }
    void load_negative_mask(const char file[])
    {
        ifstream fi(file);
        fi.read(mask, image_size);
        int count = 0;
        for(int i = 0; i < image_size; i++){
            mask[i] = (mask[i])? 0:1;
            
            if(!mask[i]){
                candidates[i] = count;
                count++;
            }
        }
 //       if(! check) check = new check_mask(mask, dim0, dim1);
//        else check->reset(mask);
        cout << "load negative"<<endl;
        candidate_size =count;
        cout << "count = "<< count<<endl;
        get_positions();
        masked = true;
    }
    
    void no_mask()
    {
         for(int i = 0; i < image_size; i++){
            if(! mask[i]){
                candidates[i] = i;
            }
        }

        candidate_size =image_size;
 //       masked = true;
    }

    void load_response(const char file[])
    {
        ifstream fi(file);
        fi.read(reinterpret_cast<char *>(response), image_size*sizeof(float));
        flat = true;
    }
    
    void get_min_max(float &min_x, float &max_x, float &min_y, float &max_y)
    {
        min_x = 1e5;
        max_x = -1e5;
        min_y = 1e5;
        max_y = -1e5;
        float *ptr = center_pos;
        for(int i = 0; i < image_size; i++, ptr+=3){
            if(ptr[2] < min_x ) min_x = ptr[2];
            if(ptr[2] > max_x ) max_x = ptr[2];
            if(ptr[1] < min_y ) min_y = ptr[1];
            if(ptr[1] > max_y ) max_y = ptr[1];
        }
    }
    
    int find_pos(float x, float y){
        float *ptr = center_pos;
        float min = 1e5;
        int pos;
        
        for(int i = 0; i < image_size; i++, ptr+=3){
            float dx, dy;
            if((dx=(x-ptr[2])) > step) continue;
            if((dy=(y-ptr[1])) > step) continue;
            float dist2 = dx*dx + dy*dy;
            if (dist2 < min){
                min = dist2;
                pos = i;
            }
        }
        if (min < step) return pos;
        else return -1;
    }
    
    void zoom(float in_step){
        step = in_step;
        get_positions();
    }
    void print_xy(int i){
        float x = center_pos[i*3+2];
        float y = center_pos[i*3+1];
        cout << "position "<<x <<" "<< y << endl;
    }
    
    void get_xy(int i, float *xy){
        xy[0] = center_pos[i*3+2];
        xy[1] = center_pos[i*3+1];
    }
    
    void get_positions()
    {
        for(int i = 0; i < view_size; i++) pixel_pos[i] = -1;
        for(int i = 0; i < view_size; i++) dist[i] = 1e5;

        int count = 0;
        float rstep = 1./step;
        cout << "get_position "<<endl;
 //       candidate_size = image_size;
 //       for(int i = 0; i < image_size; i++)
//        {
//            candidates[i] = i;
//        }
        cout << " candidate_size "<<candidate_size<<endl;
        for(int j = 0; j < candidate_size; j++){
            int i = candidates[j];
            //            if (!(i%1000)) cout << i << endl;
            float x = center_pos[i*3+2];
            float y = center_pos[i*3+1];
            float x1 = x*rstep;
            float y1 = y*rstep;
            int xi = int(x1+0.5);
            int yi = int(y1+0.5);
            float xs = x1-xi;
            float ys = y1-yi;
            
            if (xi < 0 ) continue; if (xi > (view_width2)) continue;
            if (yi < 0 ) continue; if (yi > (view_height2)) continue;

            candidates1[count] = i;
            count++;
            
            int pos = xi + view_width * yi;
            float a = xs*xs + ys*ys ; //chose the closest.
            
            if(a < dist[pos]) {
                dist[pos] = a; //currently the closest
                pixel_pos[pos] = i;
            }
        }
   //     for(int i = 0; i < view_size; i++) cout << pixel_pos[i]<<" ";
   //     cout << endl;
        cout <<"end get position"<<endl;
        cout << "candidate_size "<<candidate_size<<endl;
        candidate_size = count;
        cout << "image_size "<<image_size<<endl;
        swap();
    }
    
    void init_quick_view(int width, int height){
        view_width = width;
        view_height = height;
        view_width2 = view_width-2;
        view_height2 =view_height-2;
        step = ( max_x - min_x) / width;
        cout << " step = "<<step<<endl;
        view_size = width * height;
        pixel_pos = new int[view_size];
        dist = new float[view_size];
        cout <<" allocated pixel_pos"<<endl;
        get_positions();
        cout <<"initialized"<<endl;
    }
    
    template<class T>
    void make_image1(unsigned char *image, T *data, float x0, float x1){
        double scale = 255./(x1-x0);
  //      make_image(image, data, x0, scale);
        make_flattened_image(image, data, x0, scale);
    }
    void make_image1(unsigned char *image, int *data, int *bg, float cof, float x0, float x1){
        double scale = 255./(x1-x0);
  //      make_image(image, data, x0, scale);
        make_flattened_image(image, data,bg, cof, x0, scale);
    }
    template<class T>
    void make_flattened_image(unsigned char *image, T *data, int *bg, float cof, float x0, float scale)
    {
        if(!image)cerr<<"make_image error: no image space";
        if(!data)cerr<<"make_image error: no data space";
        for(int i = 0; i < view_size; i++){
            int ipos = pixel_pos[i];
            if(ipos<0) {
                image[i] = 0;}
            else{
                float a = ((data[ipos]*cof - bg[ipos])*response[ipos] - x0) * scale;
                if (a<0) a= 0;
                if (a>255)a=255;
                image[i] = (char)a;
            }
        }
    }
    template<class T>
    void make_flattened_image(unsigned char *image, T *data, float x0, float scale)
    {
        if(!image)cerr<<"make_image error: no image space";
        if(!data)cerr<<"make_image error: no data space";
        for(int i = 0; i < view_size; i++){
            int ipos = pixel_pos[i];
            if(ipos<0) {
                image[i] = 0;}
            else{
                float a = (data[ipos]*response[ipos] - x0) * scale;
                if (a<0) a= 0;
                if (a>255)a=255;
                image[i] = (char)a;
            }
        }
    }
    void make_image(unsigned char *image, int *data, float x0, float scale){
        //   cout <<"view size "<<view_size<<endl;
        if(!image)cerr<<"make_image error: no image space";
        if(!data)cerr<<"make_image error: no data space";
        for(int i = 0; i < view_size; i++){
            int ipos = pixel_pos[i];
            if(ipos<0) {
                image[i] = 0;}
            else{
                float a = (data[ipos] - x0) * scale;
                if (a<0) a= 0;
                if (a>255)a=255;
                image[i] = (char)a;
            }
        }
    }
    
    void make_image_zoom(unsigned char *image, int *data, float x0, float scale, int ix0, int iy0, int zoom){
        //   cout <<"view size "<<view_size<<endl;
        if(!image)cerr<<"make_image error: no image space";
        if(!data)cerr<<"make_image error: no data space";
        if(! zoom>0) cerr << "zoom "<< zoom << "is wrong";
        
        for(int i = 0; i < view_size; i++){
            int ipos = pixel_pos[i];
            if(ipos<0) {
                image[i] = 0;}
            else{
                float a = (data[ipos] - x0) * scale;
                if (a<0) a= 0;
                if (a>255)a=255;
                image[i] = (char)a;
            }
        }
    }
    
    float *get_flat(){return response;}
    int get_pos(int i){return pixel_pos[i];}
    void dump_pixel_corners(){
        float *ptr = pixel_corners;
        for(int i = 0; i < 4; i++, ptr+=3){
            cout << ptr[0] << " "<<ptr[1]<<" "<<ptr[2]<<endl;
        }
        cout <<endl;
        ptr = pixel_corners;
        for(int i = 0; i < 5; i++, ptr+=12){
            cout << ptr[0] << " "<<ptr[1]<<" "<<ptr[2]<<endl;
        }
        cout << endl;
        ptr = pixel_corners+12*dim0*120;
        for(int i = 0; i < 5; i++, ptr+=12){
            cout << ptr[0] << " "<<ptr[1]<<" "<<ptr[2]<<endl;
        }
    }
    int hide(int i){return mask[i];}
    char *get_mask(){return mask;}
    int size(){return image_size;}
    float *positions(){return center_pos;}
    bool is_masked(){return masked;}
    int ld(){return dim0;}
protected:
    friend class apply_poni;
    
    size_t dim0, dim1;
    size_t image_size;
    size_t pixel_corners_dim;
    float *pixel_corners = NULL;
    float *center_pos = NULL;
    int *pixel_pos;
    float *dist ;
    char *mask;
//    check_mask *check;
    bool flat;
    bool masked;
    char mask_filename[256];
    float *response;
    double *filtered;
    double *temp;
    double *clean_data;
    double *clean_bg;
    double *clean_subtracted;
#ifndef NOFILTER
    Filter2d *filter;
#endif
    int candidate_size;

    char response_filename[256];
    float min_x, max_x, min_y, max_y;
    float shape;
    float step;
    
    size_t view_size;
    int view_width;
    int view_height;
    int view_width2;
    int view_height2;

    int *candidates;
    int *candidates1;
    void swap(){int *temp = candidates1; candidates1 = candidates; candidates = temp;}
};

#endif /* detector_h */
