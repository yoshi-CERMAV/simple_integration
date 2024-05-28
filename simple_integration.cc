#include <unistd.h>
#include <detector.h>
#include "apply_poni.h"
#include <Read_D5WOS.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <functional>
////// you can modify here to change nuber of sampling and sampling mode
#define LOG_SAMPLING
int sampling = 300;
//////////////////////////////////////////////---------------------------------------------
int fit(vector<double> data, vector<double> &p, vector<double> &p_error);
char read_h5::D5fmt_[] = "/%d.1/measurement/D5";
char read_h5::WOSfmt_[] = "/%d.1/measurement/WOS";
char read_h5::PM1fmt_[] = "/%d.1/measurement/pm1";
char read_h5::fmtdat_[] = "/%d.1/instrument/epoch/value";

hsize_t vec_sizeD5[2]={578, 960};
hsize_t vec_sizeWOS[2]={1156, 600};

DataSpace read_h5::memspaceD5 (2, vec_sizeD5) ;
DataSpace read_h5::memspaceWOS (2, vec_sizeWOS) ;
enum D5_WOS{D5, WOS};
using namespace std;
detector *Detector;
int *sample_data = 0;
int *mt_data = 0;
apply_poni *poni;
int detector_type;
double scale = 1;
read_h5 *reader ;

FILE *gp = NULL;
int fit_list(char filename[]);

void plot(vector<float> &dat, int n){
    if(!gp)gp = popen("/opt/homebrew/bin/gnuplot --persist", "w");
    fputs("plot \"-\" binary format = \"%float%float%float\" record = ", gp);
    fprintf(gp, "%d using 1:2:3 palette\n", n);
    fwrite(reinterpret_cast<char *>(&(dat[0])), sizeof(float),3*n, gp);
    fflush(gp);
}

void plot(vector<double> &dat, int n){
    if(!gp)gp = popen("/opt/homebrew/bin/gnuplot --persist", "w");
    fputs("plot \"-\" binary format = \"%double%double%double\" record = ", gp);
    fprintf(gp, "%d using 1:2:3 palette\n", n);
    fwrite(reinterpret_cast<char *>(&(dat[0])), sizeof(double),3*n, gp);
    fflush(gp);
}

void plot(vector<float> &dat, int i, int n)
{
    if(!gp)gp = popen("/opt/homebrew/bin/gnuplot --persist", "w");
//    fputs("set cbrange [: 1e4] \n",gp);
    fputs("plot \"-\" binary format = \"%float\" record = ", gp);
    fprintf(gp, "%d  u 1", n);
    for(int k = 1; k < i; k++){
        fputs(", \"-\" binary format = \"%float\" record = ", gp);
        fprintf(gp, "%d u 1", n);
    }
    fputs("\n", gp);
    fwrite(reinterpret_cast<char *>(&(dat[0])), sizeof(float),n*i, gp);
    fflush(gp);
}

template<class T>
void read(const char fmt[], const char filename[], T *ptr, int size)
{
    char infile[255];
    snprintf(infile, 255, fmt, filename);
    ifstream fi(infile);
    fi.read(reinterpret_cast<char *>(ptr), sizeof(T)*size);
}

void read_paths(const char filename[], vector<string> &name)
{
   ifstream fi(filename);
   string s;
   while(fi >>s) name.push_back(s); 
}

void dump_data(const char filename[], void *ptr, int image_size)
{
    ofstream fo(filename);
    fo.write(reinterpret_cast<char *>(ptr), image_size*4);
}

void finalize_flat(float *flat, char *mask, int n, float min, float  max)
{
    int count0 = 0;
    int count1 = 0;
    for(int i = 0; i != n; i++){
        if(flat[i]< min) {flat[i] = 0; mask[i] = 1;count0++;}
        if(flat[i]> max) {flat[i] = 0; mask[i] = 1;count1++;}
    }
    cout << "pixels under "<<min << ": " << count0<< endl;
    cout << "pixels above "<<max << ": " << count1<< endl;
}

template<class T>
void readline(ifstream &fi, T &val)
{
    char buffer[256];
    fi.getline(buffer,256);
    istringstream istr(buffer);
    istr >> val;
}

int robust_average(vector<double> &data, double &avg, double &stdev, int &count)
{
    double *ptr = &data[0];
    int quarter = data.size()/4;
    double *ptr1 = ptr+quarter;
    double *ptr2 = ptr+data.size();
    nth_element(ptr, ptr1, ptr2);
    nth_element(ptr1, ptr2-quarter, ptr2);
    double sum=0;
    double sum2 = 0;
    int n = data.size()-2*quarter;
    if(!n) return 0;
    for(int i = 0; i < n; i++){
        sum+= ptr1[i];
        sum2+= ptr1[i] * ptr1[i];
    }
    sum/=n; sum2/=n;
    stdev = sqrt(sum2- sum*sum);
    avg = sum;
    double sigma3 = stdev*3;
    count = 0;
    sum = 0;
    sum2 = 0;
    for(vector<double>::iterator itr = data.begin(); itr!= data.end(); itr++){
        if(fabs((*itr)-avg) > sigma3) continue;
        sum += (*itr);
        sum2 += (*itr)*(*itr);
        count++;
    };
    if(!count) return 0;
    sum/=count;
    sum2/=count;
    avg = sum;
    stdev = sqrt(sum2-sum*sum);
    return 1;
}


int average_frame(int frame, int n, vector<double> &q, vector<double> &avg, vector<double> &stdev, vector<int> &count){
    q.resize(n);
    avg.resize(n);
    stdev.resize(n);
    count.resize(n);
    if(detector_type==D5_WOS::WOS){
        reader->accumulate_scan_WOS(frame, sample_data);
    }else{
        reader->accumulate_scan_D5(frame, sample_data);
    }
    double q0 = poni->qmin();
    double q1 = poni->qmax();
#ifdef LOG_SAMPLING
    double logq0 = log(q0);
    double logq1 = log(q1);
    double logqstep = (logq1-logq0)/n;
#else
    double qstep = (q1-q0)/n;
#endif
//    double qstep = (q1-q0)/n;
    vector<double> data;
    int pos = 0;
    for(int i = 0;i < n; i++){
        data.clear();
#ifdef LOG_SAMPLING
        q1 = exp(logq0 + logqstep);
        q0 = exp(logq0);
        logq0 = log(q1);
#else
        q1 = q0+qstep;
#endif
        cout << q0 << " "<<q1<<endl;
        q[i] = poni->get_data(data, sample_data, pos, q0,q1);
#ifndef LOG_SAMPLING
        q0 = q1;
#endif
        if(q[i]){
            double avg_temp, stdev_temp;
            int n;
            if(!data.size()) continue;
            robust_average(data, avg_temp, stdev_temp, n);
            //      q[i] = exp(log(q0)+0.5*logqstep;
            avg[i] = avg_temp;
            stdev[i] = stdev_temp;
            count[i] = n;
        }
    }
}
inline
void write(ofstream &fo, double &a, double &b, double &c, int &d)
{
    fo.write(reinterpret_cast<char *>(&a), sizeof(double));
    fo.write(reinterpret_cast<char *>(&b), sizeof(double));
    fo.write(reinterpret_cast<char *>(&c), sizeof(double));
    fo.write(reinterpret_cast<char *>(&d), sizeof(int));
}

int average_list(char filename[])
{
    ifstream fi(filename);
    int i;
    string outfile;
    vector<double> q;
    vector<double> avg;
    vector<double> stdev;
    vector<int> count;
    char fname[256];
    while (fi>> i ){
//        cout << i << " "<<outfile<<endl;
        if(detector_type==D5_WOS::WOS){
            snprintf(fname, 256, "%dWOS.dat", i);
        }else{
            snprintf(fname, 256, "%dD5.dat", i);
        }
        average_frame(i, sampling, q, avg, stdev, count);
        cout << "averaged"<<endl;
        ofstream fo(fname);
//        fo.write(reinterpret_cast<char *>)
        for(int i = 0; i < q.size(); i++){
            if(q[i])
            write(fo, q[i], avg[i], stdev[i], count[i] );
//            fo << q[i]<<" "<<avg[i]<<" "<<stdev[i]<<" "<<count[i]<<endl;
        }
    }
}


int main(int argc, char *argv[])
{
    char detector_name[256];
    char file_name[255];
    char poni_name[256];
    char flat_fname[256];
    char mask_fname[256];
    char profile_fname[256];
    int snap_h, snap_v;
    char snap_fname[256];
    int sample_scan;
    int mt_scan;
    double scale2;
    double mask_max;
    double mask_min;
    double image_max;
    double image_min;
    ifstream fi(argv[1]);
    readline(fi, detector_name);
    regex wos("WOS", regex_constants::icase);
    if(regex_search(detector_name, wos)) detector_type = D5_WOS::WOS;
    else detector_type=D5_WOS::D5;
    
    readline(fi, file_name);
    readline(fi, poni_name);
    readline(fi, flat_fname);
    readline(fi, mask_fname);
//    readline(fi, sample_scan);
    cout << "detector_name : "<<detector_name<<endl;
    Detector = new detector(detector_name);
    Detector->load_mask(mask_fname);
    Detector->load_response(flat_fname);
    
    reader = new read_h5(file_name);
    char buffer[256];
    
    int image_size;
    if(detector_type==D5_WOS::WOS){
        cout <<"WOS"<<endl;        
        image_size = vec_sizeWOS[0] * vec_sizeWOS[1];
        sample_data = new int[image_size];
        mt_data  = new int[image_size];
//        reader->accumulate_scan_WOS(sample_scan, sample_data);
    }else{
        cout <<"D5"<<endl;
        image_size = vec_sizeD5[0] * vec_sizeD5[1];
        sample_data = new int[image_size];
        mt_data  = new int[image_size];
 //       reader->accumulate_scan_D5(sample_scan, sample_data);
    }

    poni = new apply_poni(poni_name, Detector);
    poni->sort_index_by_q();
    average_list(argv[2]);
 }
