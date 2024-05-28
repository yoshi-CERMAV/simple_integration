//
//  apply_poni2003.h
//  
//
//  Created by Yoshiharu Nishiyama on 1/5/2024.
//

#ifndef apply_poni_h
#define apply_poni_h
#include <detector.h>
#include <map>
#include <rotator.h>
#include <Accelerate/Accelerate.h>
#include <vector>
#include <tile_grid.h>

const double four_pi = M_PI * 4;
const double two_pi = M_PI*2;
#ifndef PANEL
#define PANEL
class panel:public vector<int>
{
public:
    panel(){}
protected:
    bool is_empty;
};
#endif

class Comparator
{
public:
    Comparator(float *d){dat = d;}
    bool operator() (size_t i, size_t j){return dat[i*2] < dat[j*2];}
private:
    float *dat;
};

class Comparator_by_b
{
public:
    Comparator_by_b(float *d){dat = d;}
    bool operator() (size_t i, size_t j){return dat[i*2+1] < dat[j*2+1];}
private:
    float *dat;
};

class Comparator_sector
{
public:
    Comparator_sector(float *d, int sector, float offset){
        dat = d;
        sector_width = two_pi/sector;
        offset_=offset;
    }
    bool operator() (size_t i, size_t j){
        int s1 = sector(i);
        int s2 = sector(j);
        if(s1 < s2) return true;
        else if(s1 > s2) return false;
        else return (dat[i*2] < dat[j*2]);
    }
    inline int sector(int i){
        float t = dat[i*2+1];
        float s = (t + offset_);
        while (s < 0) s += two_pi;
        while (s > two_pi) s-= two_pi;
        return s/sector_width;
    }
    float *dat;
    float sector_width;
    float offset_;
    
};

class sector{
public:
    sector(){set_null();}
    void set_null(){i0=0; i1=0; num_points=0;}
    void print(){cout << i0<<" "<< i1 << " "<< num_points<<endl;}
    int i0, i1, num_points;
};

class apply_poni
{
public:
    apply_poni(detector *in_detector){
        Detector = in_detector;
        allocate(Detector->image_size);
        mask =new char[Detector->image_size];
        for(int i = 0; i <Detector->image_size; i++ ) mask[i] = 0;
    }
    
    apply_poni(const char poni_filename[], detector *in_detector){
        Detector = in_detector;
        allocate(Detector->image_size);
        read_poni(poni_filename);
        dim0=Detector->dim0;
        init_map();
        check_flat();
     }
    
    void read_poni(const char poni_filename[])
    {
        char key[255];
        double value;
        char *buffer = new char [1024];
        
        ifstream ponifile(poni_filename);
        
        while(ponifile.getline(buffer, 1024)){
            sscanf(buffer,"%[^:]:%lf",key, &value);
            poni[key]=value;
            cout << key <<": "<< value <<endl;
        }
        cout <<": "<<  poni["Poni1"]<<endl;
        r.init();
        r.roty(-poni["Rot1"]);
        r.rotx(-poni["Rot2"]);
        r.swap(); // swapping x and z;
        center[0] = poni["Distance"];
        center[1] = -poni["Poni1"];
        center[2] = -poni["Poni2"];
        cout << "center "<<center[0] <<" "<<center[1]<<" "<<center[2]<<endl;
        wavelength = poni["Wavelength"];
    }
 
    void read_poni_and_init(const char poni_filename[])
    {
        read_poni(poni_filename);
        init_map();
        sort_index_by_q();
    }
    
    void check_flat(){
        if(Detector->flat){
            response = Detector->response; // just copying the pointer
            for(int i = 0; i < Detector->image_size; i++ ){
                correct_response[i] = response[i] * center_corr[i];
            }
        }
        else {
            for(int i = 0; i < Detector->image_size; i++ ){
                correct_response[i] =  center_corr[i];
            }
            response = NULL;
        }
    }
 
    void print_q(int pos){
        cout << "q = "<<center_qb[pos*2]<<endl;
    }
    

   
    void allocate(int image_size){
        pixel_index = new size_t[image_size];
        center_qb = new float[image_size * 2];
        center_corr = new float[image_size];
        correct_response = new float[image_size];
        work = new double [image_size];
    }
    
    double get_q(int x, int y){
        size_t pos = y*dim0 + x;
        return center_qb[pos*2];
    }
    
    void shift_center_1(double x){center[1] += x;}
    void shift_center_2(double x){center[2] += x;}

    void copy3(float *a, float *b){
        a[0] = b[0];a[1] = b[1];a[2] = b[2];
    }
    void reinit_map(vector<int> & index){
        float temp[3];
        float scale = four_pi / wavelength * 1e-10;

        for(int i = 0; i< index.size(); i++){
            int pos = index[i];
            copy3(temp, Detector->center_pos+pos*3);
            add(temp,center);
            r.apply(temp);
            double r2 = sqrt(temp[0] * temp[0] + temp[1] * temp[1]);
            double Theta = atan2(r2, temp[2]);
            center_qb[pos*2] = scale * sin(0.5*Theta);
        }
    }
    
    void init_map()
    {
        float *ptr = Detector->center_pos;
        float scale = four_pi / wavelength * 1e-10;
        float *ptr_qb = center_qb;
        float norm[3];
        norm[0] = 1; norm[1] = 0; norm[2] = 0;
        r.apply(norm);
        float temp[3];
        cout << "center "<<center[0] <<" "<<center[1] <<" "<<center[2] << endl;
        for(int i = 0; i < Detector->image_size; i++, ptr+=3, ptr_qb+=2){
            for(int i = 0; i < 3; i++) temp[i] = ptr[i];
            add(temp, center);
            r.apply(temp);
            double r2 = sqrt(temp[0] * temp[0] + temp[1] * temp[1]);
            double Theta = atan2(r2, temp[2]);
            ptr_qb[0] = scale * sin(0.5*Theta);
            double beta = atan2(temp[0], temp[1]);
            while(beta<0) beta+=two_pi;
            ptr_qb[1] =beta;
            
            double cT = cos(Theta);
            double sT = sin(Theta);
            cT*=cT; sT*=sT;
            float corr = 0.5*(1 + cT + 1*cos(beta * 2) * sT); // phase 90˚ diff from -cos..
            corr = 1./corr;
            double r3 =temp[0] * temp[0] + temp[1] * temp[1] + temp[2]*temp[2];
            corr *= fabs(r3*sqrt(r3)/scalar(temp, norm));
            center_corr[i] = corr;
        }
        init_index();
        sort_index_by_q();
    }
    
    void init_index_nomask()
    {
        for(int i = 0; i!= Detector->image_size; i++) pixel_index[i]=i;
        valid_count = Detector->image_size;
    }
    void init_index_mask(char *mask){
        valid_count = 0;
        for(int i = 0; i!= Detector->image_size; i++) {
            if(mask[i]) continue;
            else{
                pixel_index[valid_count] = i;
                valid_count++;
            }
        }
    }
    void init_index()
    {
        if(Detector->masked){
            init_index_mask(Detector->mask);
        }
        else{
            init_index_nomask();
        }
    }
    
    template<class T>
    int get_data(vector<int> index, T *dat, double tm, double qm, vector<double> &data){
        data.clear();
        for(int i = 0; i < index.size(); i++){
            size_t p = index[i];
            size_t p2 = p+p;
            double q = center_qb[p2] - qm;
            double t = center_qb[p2+1] - tm;
            while (t > M_PI) t-= two_pi;
            while (t < -M_PI) t+= two_pi;
            double Y = dat[p]* correct_response[p];
            if( isnan(Y) || dat[p]< 1 ||mask[p]) continue; // fit_mask[p] removed;
            data.push_back(t);
            data.push_back(q);
            data.push_back((Y));
        }
    }
    
    
    template <class T1, class T2>
    void get_data(vector<T1> &out_data, vector<int> &index,  T2 *in_data)
    {
        out_data.clear();
        for(int i = 0; i < index.size(); i++){
            int pos = index[i];
            int pos2 = pos+pos;
            out_data.push_back( center_qb[pos2+1]);
            out_data.push_back( center_qb[pos2]);
            out_data.push_back( in_data[pos] * correct_response[pos]);
        }
    }

    
    template <class T1, class T2>
    void get_data(T1 *out_data, int *index, int n, T2 *in_data)
    {
        T1  *ptr = out_data;
        for(int i = 0; i < n; i++, ptr +=3){
            int pos = index[i];
            int pos2 = pos+pos;
            ptr[0] = center_qb[pos2];
            ptr[1] = center_qb[pos2+1];
            ptr[2] = in_data[pos] * correct_response[pos];
        }
    }
    
    template <class T1, class T2>
    void get_data(T1 *out_data, int i0, int i1, T2 *in_data)
    {
        size_t *index = pixel_index+i0;
        int n = i1-i0;
        get_data(out_data, index, n, in_data);
    }
    
    template<class T1, class T2>
    double get_data(vector<T1> &out_data, T2 *in_data, int &i, double q0, double q1)//pixel_index has to be sorted
    {
        while(center_qb[pixel_index[i]*2]< q0 && i < valid_count) i++;
        double sum = 0;
        double count = 0;
        while(center_qb[pixel_index[i]*2]< q1 && i < valid_count) {
            int pos = pixel_index[i];
            out_data.push_back(in_data[pos] * correct_response[pos]);
            sum += center_qb[pos*2];
            i++;count++;
        }
        if(count)
        return sum/count;
        else return 0;
    }
    template<class T1, class T2>
    void get_data(vector<T1> &out_data, T2 *in_data, int &i, double q0, double q1, double scale, Tile_Grid &tg)
    {
        while(center_qb[pixel_index[i]*2]< q0 && i < valid_count) i++;
        while(center_qb[pixel_index[i]*2]< q1 && i < valid_count) {
            int pos = pixel_index[i];
  //          cout <<center_qb[pos*2] <<" "<<center_qb[pos*2+1]<<" "<<in_data[pos] * correct_response[pos] * scale<<" "<<tg.evaluate(center_qb[pos*2+1], center_qb[pos*2])<<endl;
            double v = in_data[pos] * correct_response[pos] * scale - tg.evaluate(center_qb[pos*2+1], center_qb[pos*2]);
            out_data.push_back(v);
            i++;
        }

    }

    
    void sort_index_by_q()
    {
        cout << "valid pixels" <<valid_count <<endl;
        Comparator comp(center_qb);
        sort(pixel_index, pixel_index + valid_count, comp);
        cout <<"sorted" <<endl;
        auto [q, t] = std::div((int)pixel_index[0], Detector->ld());
        cout << "q t "<< q<<" "<<t<<" "<<center_qb[pixel_index[0]*2]<<endl;
        auto [q1, t1] = std::div((int)pixel_index[valid_count-1], Detector->ld());
        cout << "q t "<< q1<<" "<<t1<<" "<<center_qb[pixel_index[valid_count-1]*2]<<endl;
    }
    
    void trim_index_azim(double xmin, double xmax, vector<int> &index_in, vector<int> &index_out){
        index_out.clear();
        for(vector<int>::iterator itr = index_in.begin();itr != index_in.end(); itr++){
            double x = center_qb[2*  (*itr) + 1 ];
            if(x>xmin && x < xmax) index_out.push_back(*itr);
        }
    }
    void sort_index_by_sector_q(int n_sector, float offset){
        Comparator_sector comp(center_qb, n_sector, offset);
        cout<<"comparator "<<endl;
        sort(pixel_index, pixel_index+valid_count, comp);
        cout <<"sortd"<<endl;
        if(! v_sector.size()) v_sector.resize(n_sector);
  //      if(! v_sector.size()) v_sector.reserve(n_sector);

//        cout << "init"<<endl;
//        if(comp.sector(pixel_index[0]) == 0) sectors[0].i0;
        
        int curr_sector = comp.sector(pixel_index[0]);
        if (curr_sector) for(int i = 0; i != curr_sector; i++) {
            v_sector[i].set_null();
        }
        cout << curr_sector <<endl;
        for(int i = 1; i < valid_count; i++){
            int s = comp.sector(pixel_index[i]);
            if (s!=curr_sector){
                v_sector[curr_sector].i1 = i;
                v_sector[curr_sector].num_points = i-v_sector[curr_sector].i0;
                cout << s<<endl;
                if(s!=n_sector){
                    v_sector[s].i0 = i;
                    for(int j = curr_sector+1; j!=s; j++) {
                        v_sector[j].set_null();
                    }
                }
                curr_sector = s;
            }
        }
        if(curr_sector < n_sector){
            v_sector[curr_sector].i1 = valid_count;
            v_sector[curr_sector].num_points = valid_count-v_sector[curr_sector].i0;
            for(int j = curr_sector+1; j != n_sector; j++){
                v_sector[j].set_null();
            }
        }
        for(int i = 0; i< n_sector; i++){
            cout << i <<" ";
            v_sector[i].print();
        }
    }
    
    int copy_index(int i0, int i1, int *index){
        int count = 0;
        for(int i = i0; i < i1; i++, count++)index[count] = pixel_index[i];
        return count;
    }
    
    double qmin(){return center_qb[pixel_index[0]*2];}
    double qmax(){//cout << "valid_count "<<valid_count <<endl;
        return center_qb[pixel_index[valid_count-1]*2];}
    
    int find_i_range(int sector, double q0, double q1, int &i0, int &i1){
        return find_i_range(q0, q1, v_sector[sector].i0, v_sector[sector].i1, i0, i1);
    }
    
    int find_i_range(double q0, double q1, int imin, int imax, int &i0, int &i1)
    {
        i0=imin;
        while (center_qb[pixel_index[i0]*2] < q0  && i0<imax) i0++;
        if (i0==imax){ return 0; }//cerr<< "i0==imax"<<endl; return 0;}
        i1 = i0;
        while(center_qb[pixel_index[i1]*2] < q1 &&i1 < imax) i1++;
        return i1-i0;
    }
    
    int find_i_range(double q0, double q1, int &i0, int &i1){
        cout << "find irange"<<endl;
        i0 = 0;
        while (center_qb[pixel_index[i0]*2] < q0 && i0 < valid_count) i0++;
        if (i0==valid_count){cerr<< "i0==valid_count "<< valid_count<<endl; return 0;}
        i1 = i0;
        while (center_qb[pixel_index[i1]*2] < q1 && i1 < valid_count) i1++;
        return i1-i0;
    }
    
    int make_tiles(int bn, int qn, vector<panel> & v_panel){
        double q_min = qmin();
        double qstep = (qmax()-q_min)/qn;
        double bstep = two_pi/bn;
        double q1= qstep + q_min;
        int i0 = 0;
        Comparator_by_b comp_b(center_qb);
        v_panel.resize(qn*bn);
        vector<panel>::iterator panel_ptr = v_panel.begin();
        int i1 = 0;
        for(int i = 0; i < qn; i++, q1 += qstep){
            while (center_qb[pixel_index[i1]*2] < q1  && i1< valid_count) i1++;
            sort(pixel_index+i0, pixel_index+i1, comp_b);
            double b = bstep;
            int j0 = i0;
            i0 = i1;

            for(int j = 0; j < bn; j++, panel_ptr++, b+= bstep){
                panel_ptr->clear();
                while(center_qb[pixel_index[j0]*2+1] < b && j0 < i1 ){
                    panel_ptr->push_back(pixel_index[j0]);
                    j0++;
                }
            }
        }
    }
    
    void add(float *a, float *b)
    {
        a[0] += b[0];
        a[1] += b[1];
        a[2] += b[2];
    }
    float sq(float *x)
    {
        return x[0] * x[0] + x[1]*x[1] + x[2]*x[2];
    }
    float scalar(float *a, float *b)
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    void write(ofstream &fo, double data){
        fo.write(reinterpret_cast<char *>(&data), sizeof(double));
    }
    void write(ofstream &fo, int data){
        fo.write(reinterpret_cast<char *>(&data), sizeof(int));
    }
    
    void select_qrange(double qmin, double qmax, int *dat, float *out_data, int &count)
    {
        check_flat();
        if(Detector->flat) select_qrange(qmin, qmax, dat, out_data, count, &apply_poni::get_flattened_value);
        else select_qrange(qmin, qmax, dat, out_data, count, &apply_poni::get_value);
    }
    
    typedef  float (apply_poni::*FreeFn)(int p);
    float get_value(int pos){return data[pos];}
    float get_flattened_value(int pos){return data[pos] * response[pos];}

    int get_q_range_index(double qmin, double qmax, vector<int> &index){
        int i = 0;
        index.clear();
        while(center_qb[pixel_index[i]*2] < qmin && i < valid_count) i++;
        if(i== valid_count) return 0;
        while(center_qb[pixel_index[i]*2] < qmax && i < valid_count){
            index.push_back(pixel_index[i]);
            i++;
        }
        return index.size();
    }
    
    void select_qrange(double qmin, double qmax, int *dat, float *out_data, int &count, FreeFn func)
    {
        int i = 0;
        float *ptr = out_data;
        int count1 = 0;
        double q;
        while(center_qb[pixel_index[i]*2] < qmin) {
            i++;
            //            if (i> valid_count)  need to implement exception throw
        }
        q = center_qb[pixel_index[i]*2];
        while(q < qmax && i < valid_count){
            int pos =pixel_index[i];
            int pos2 = pos*2;
            ptr[0] = center_qb[pos2];
            ptr[1] = center_qb[pos2+1];
            ptr[2] = (this->*func)(pos);
            count1++;
            //need to check if count1 < count;
        }
        count = count1;
    }
    
    void add_mask(const char *mask1){
        for(int i = 0; i != Detector->image_size; i++) if(mask1[i]) mask[i] = 1;
    }
    
    void mask_pos(int i){mask[i] = 1;}
    
    int mask_line(int n){
        char *ptr = mask + n*dim0;
        for(int i = 0; i != dim0; i++ ) ptr[i] = 1;
        return 0;
    }
    
    template<class T1, class T2>
    int get_data_from_sector(int isec, vector<T1> &dat, T2* data, int &npoint)
    {
        if(isec<0){ cerr<< "sector num < 0"; return 0;}
        cout << isec << " "<< v_sector.size()<<endl;
        if(isec >= v_sector.size()){cerr <<"sector num out of range"<<endl; return 0;}
        npoint = v_sector[isec].num_points;
        cout << "npoint "<<npoint<<endl;
        int size = npoint*3;
        if(size > dat.size() )dat.resize(size);
        int i0 = v_sector[isec].i0;
        for(int i = 0; i < size; i+=3, i0++){
            int p = pixel_index[i0];
            dat[i] = center_qb[2*p];
            dat[i+1] = center_qb[2*p+1];
            dat[i+2] = data[p];
        }
    }
    
    template <class T>
    void make_flat(double (*interpolate)(double q), T *data, float *flat)
    {
        for(int i = 0; i < Detector->image_size; i++){
            if(!data[i] || isnan(data[i])) {
                flat[i] = 0;
            }else{
                flat[i] = interpolate(center_qb[i*2])/(center_corr[i] * data[i]);
            }
        }
    }



protected:
     
    inline void write_float(float *dat){fo.write(reinterpret_cast<char *>(dat), sizeof(float));}
    inline void write_int(int *dat){fo.write(reinterpret_cast<char *>(dat), sizeof(int));}
    inline void increment_q()
    {
        q0 += qstep;
        q1 += qstep;
        current_q += qstep;
    }
    
    
    inline void init_file(const char filename[])
    {
        fo.open(filename);
        do_write=true;
    }
    
    double *work;
    int *data;
    int *bg;
    char *mask;
    


    float *response = NULL;
    float current_q;
    double qstep;
    double q0, q1;
    int nq;
    double span;
    int kmax;
    ofstream fo;
    bool do_write;
    float *correct_response;
    
    float *center_qb = NULL;
    float *center_corr = NULL;
    size_t *pixel_index = NULL;
    size_t valid_count;
    int dim0;
    map<string, double> poni;
    rotator<float> r;
    float rot[3];
    float center[3];
    vector<sector>  v_sector;
    double wavelength;
    detector *Detector;
};


#endif /* apply_poni2003_h */
