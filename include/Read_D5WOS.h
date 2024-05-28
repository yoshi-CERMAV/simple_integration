//
//  Read_D5WOS.h
//  D2AM2023
//
//  Created by yoshi on 13/09/2023.
//

#ifndef Read_D5WOS_h
#define Read_D5WOS_h
#define YN_MAX_LEN 256
#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>
#include <cstring>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <vector>
#include <numeric>
using namespace std;
using namespace H5;

class read_h5
{
public:
    int get_num_frame(int i, char fmt[]){
        char name[256];
        snprintf(name, 256, fmt,i);
        hsize_t dims_out[3];
        DataSpace dataspace = file_->openDataSet(name).getSpace();
        dataspace.getSimpleExtentDims(dims_out, NULL);
        return dims_out[0];
    }
    
    DataSet open(int i, char fmt []){
        char name[256];
        snprintf(name, 256, fmt, i );
        return file_->openDataSet(name);
    }
    
    float accumulate_pm(int i){
        DataSet dataset = open(i, PM1fmt_);
        DataSpace dataspace = dataset.getSpace();
        int n= dataspace.getSimpleExtentNpoints();
        float *pm = new float[n];
        dataset.read(pm, PredType::IEEE_F32LE, dataspace, dataspace);
        return accumulate(pm, pm+n, 0);
    }
    
    float pm1(int i){
        char name[256];
        snprintf(name, 256, PM1fmt_,i);
        DataSet dataset = file_->openDataSet(name);
        DataSpace dataspace = dataset.getSpace();
        float pm;
        dataset.read(&pm, PredType::IEEE_F32LE, dataspace, dataspace);
        return pm;
    }
    static herr_t
    scan_count(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
    {
        read_h5 *reader = static_cast<read_h5 *>(opdata);
        reader->num_scans++;
        return 0;
    }

    read_h5(vector<string> file_names){
        strncpy(filename_, file_names[0].c_str(), YN_MAX_LEN);
        file_ = new H5File(filename_, H5F_ACC_RDONLY);
        num_scans = 0;
        herr_t idx = H5Literate(file_->getId(), H5_INDEX_NAME, H5_ITER_INC, NULL,  scan_count, this);
    }
    
    void open(const char filename[])
    {
        file_->close();
        strncpy(filename_, filename, YN_MAX_LEN);
        file_->openFile(filename_, H5F_ACC_RDONLY);
        num_scans = 0;
        herr_t idx = H5Literate(file_->getId(), H5_INDEX_NAME, H5_ITER_INC, NULL,  scan_count, this);

    }
    read_h5(const char filename[]){
        strncpy(filename_,filename, YN_MAX_LEN);
        file_ = new H5File(filename, H5F_ACC_RDONLY);
        num_scans = 0;
        cout << "going to count"<<endl;
        herr_t idx = H5Literate(file_->getId(), H5_INDEX_NAME, H5_ITER_INC, NULL,  scan_count, this);
        cout << num_scans<<endl;
    }

    int set_path(const char path[])
    {
        struct stat sb;
        if(stat(path, &sb)) {cerr<<"directory does not exist"; return 0;}
        else strncpy(path_, path, YN_MAX_LEN);
        return 0;
    }
    int set_filename(const char filename[])
    {
        strncpy(filename_,filename, YN_MAX_LEN);
        struct stat buf;
        if(stat(filename, &buf)) {cerr<<"directory does not exist"; return 1;}
        if(!stat(filename, &buf)){
            if(file_) delete file_;
            file_ = new H5File(filename, H5F_ACC_RDONLY);
       }
        cout <<"fmt"<<D5fmt_<<endl;
        return 0;
    }
    
    int read_D5(int *data){
        DataSet dataset;
        dataset = file_->openDataSet("/entry_0000/ESRF-BM02/D5/data");
        DataSpace dataspace = dataset.getSpace();
        dataset.read(data, PredType::STD_I32LE, dataspace, dataspace);
    }
    
    int read_D5(int i, int *data){
        return read_data(i, data, D5fmt_, memspaceD5);
    }
    int read_WOS(int i, int *data){
        return read_data(i, data, WOSfmt_, memspaceWOS);
    }
    double epochtime(int i){
        char format[255];
        double data;
        DataSet dataset;
        snprintf(format, 255, fmtdat_,i);
        dataset = file_->openDataSet(format);
        DataSpace dataspace = dataset.getSpace();
        dataset.read(&data, PredType::IEEE_F64LE, dataspace, dataspace);
        return data;
    }
    
    DataSet open_dataset(char name[])
    {
        DataSet dataset;
        try {dataset = file_->openDataSet(name);}
        catch(FileIException &E){cerr <<"could not open dataspace";exit(0);}
        return dataset;
    }
    
    void add(int *data, int *temp_data, int image_size){
        for(int i = 0; i < image_size; i++){
            data[i] += temp_data[i];
        }
    }
    void reset(int *data, int image_size){
        for(int i = 0; i < image_size; i++) data[i] =0;
    }
    
    int accumulate_scan(DataSet &dataset, DataSpace &dataspace, hsize_t dims_out[], int *data)
    {
        hsize_t vec_size[2]={dims_out[1],dims_out[2]};
        DataSpace memspace (2, vec_size);
        int image_size = vec_size[0]*vec_size[1];
        int *temp_data = new int[image_size];
        reset(data, image_size);
        hsize_t dataCount[3] = {1,dims_out[1],dims_out[2]};
        hsize_t dataOffset[3] = {0, 0, 0};

        for(int i = 0; i < dims_out[0]; i++){
            dataOffset[0] = i;
            dataspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
            dataset.read(temp_data, PredType::STD_I32LE, memspace, dataspace);
            add(data, temp_data, image_size);
        }
    }
    
    int accumulate_scan_D5(int i, int *data){return accumulate_scan(i, data, D5fmt_);}
    int accumulate_scan_WOS(int i, int *data){return accumulate_scan(i, data, WOSfmt_);}

    int accumulate_scan(int i, int *data, const char fmt[])
    {
        if (!file_)set_filename (filename_);
        snprintf(name, YN_MAX_LEN, fmt, i);
        DataSet dataset = open_dataset(name);
        DataSpace dataspace = dataset.getSpace();
        hsize_t dims_out[3];
        dataspace.getSimpleExtentDims(dims_out, NULL);
        cout << dims_out[0] <<" "<<dims_out[1]<<" "<< dims_out[2]<<endl;
        if(dims_out[0]==1) dataset.read(data, PredType::STD_I32LE, dataspace, dataspace);
        else{ 
            accumulate_scan(dataset, dataspace, dims_out, data);
        }
    }
    
    int read_data(int i, int *data, const char fmt[], DataSpace &memspace){

        snprintf(name, YN_MAX_LEN, fmt, i);
        DataSet dataset;
        struct stat buf;

        if (!file_)set_filename (filename_);
        try {dataset = file_->openDataSet(name);}
        catch(FileIException &E){cerr <<"could not open dataspace";return 2;}
        if( !dataset.getId() )
            {
                cerr<<"ReportReaderHDF5: "
                     <<"Dataset " << name << " not found "
                     <<"in file: " << file_->getFileName();
                return -1;
            }
        DataSpace dataspace = dataset.getSpace();
//        size_t required_size = dataset.getStorageSize();
//        int rank = dataspace.getSimpleExtentNdims();
//        cout <<"rank "<< rank;
//        required_size =dataspace.getSimpleExtentNpoints();
//        cout << "size "<<required_size<<endl;;
//        cout << "required_size "<< required_size<<" "<< max_size<<endl;
//        cout << required_size / max_size <<endl;
//       if(required_size > max_size){cerr <<"not enough space "<<required_size  ; return 0;}
//        if(!data) {cout<< "no memory "<<endl; return 1;}
        dataset.read(data, PredType::STD_I32LE, memspace, dataspace);
        return 0;
    }
    ~read_h5(){file_->close();}
    int number_of_frames(){return num_frame;}
    int number_of_scans(){return num_scans;}

    char * filename(){return filename_;}
private:
    char path_[YN_MAX_LEN];
    char filename_[YN_MAX_LEN];
    char name[YN_MAX_LEN];
    static char D5fmt_[];// = "/%d.1/measurement/D5";
    static char WOSfmt_[];// = "/%d.1/measurement/WOS";
    static char PM1fmt_[];// = "/%d.1/measurement/pm1";
    static char fmtdat_[];// = "/%d.1/instrument/epoch/value";
    static DataSpace memspaceD5 ;
    static DataSpace memspaceWOS ;
    bool file_set = false;
    H5File *file_;
    int num_frame;
    int num_scans;
};

#endif /* Read_D5WOS_h */
