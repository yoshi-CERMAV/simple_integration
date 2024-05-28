//
//  tile_grid.hpp
//  
//
//  Created by Yoshiharu Nishiyama on 16/05/2024.
//

#ifndef tile_grid_hpp
#define tile_grid_hpp

#include <stdio.h>
#include <vector>
#include <Accelerate/Accelerate.h>
//#include <apply_poni.h>
#include <iostream>
#include <fstream>
#include <smoothing_spline.h>
class apply_poni;
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
//int calc_cof(double **y, double *out_a,double step[2]);
int calc_cof(double **y, double *out_a,double xstep, double ystep);

using namespace std;
class Tile{
public:
    Tile(double *grid_array, int grid[4], double step[2])
    {
        double *y[4];
        for(int i = 0; i < 4; i++){y[i] = grid_array + 4 * grid[i] ;}
        calc_cof(y, a, step[0], step[1]);
    }
    
    void init(double *grid_array, int grid[4], double xstep, double ystep){
        double *y[4];
        for(int i = 0; i < 4; i++){y[i] = grid_array + 4 * grid[i] ;}
        calc_cof(y, a, xstep, ystep);
    }
    
    inline
    void fill(double x, double *fx){
        fx[0] = 1;
        for(int i = 1; i< 4; i++) fx[i] =fx[i-1] * x;
    }
    double cal(double x, double y){
        double fx[4], fy[4];
        double temp[4];
        fill(x, fx);
        fill(y, fy);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 4, 4, 1., a, 4, fy, 1, 0., temp, 1);
        return cblas_ddot(4, fx, 1, temp, 1);
    }
    void write(ofstream &fo){fo.write(reinterpret_cast<char *>(a), sizeof(double)*16);}
protected:
    double a[16];
};

class Tile_data{
    int grid[4];
};

class Tile_Grid
{
public:
    Tile_Grid(const char filename[]){
        ifstream fi(filename);
        fi.read(reinterpret_cast<char *>(size), sizeof(int)*2);
        fi.read(reinterpret_cast<char *>(step), sizeof(double)*2);
        fi.read(reinterpret_cast<char *>(offset), sizeof(double)*2);
        cout <<"size "<< size[0] << " "<<size[1]<<endl;
        cout <<"step "<< step[0] << " "<<step[1]<<endl;
        cout <<"offset "<< offset[0] << " "<<offset[1]<<endl;

        fi.read(reinterpret_cast<char *>(&num_grids),sizeof(int)); // number of valid grid poits;
        fi.read(reinterpret_cast<char *>(&num_tiles),sizeof(int)); // number of valid tiles;
        cout << num_grids<< " "<<num_tiles <<endl;
        int len = size[0] *size[1];
        tile_num = new int[len]; // number of total tiles;
        fi.read(reinterpret_cast<char *>(tile_num), len*sizeof(int)); // point to the valid tiles;
        tile_corners = new int[num_tiles *4];
        fi.read(reinterpret_cast<char *>(tile_corners), 4*num_tiles*sizeof(int));
        grid_data=new double[num_grids*4];
        fi.read(reinterpret_cast<char *>(&(grid_data[0])), num_grids * 4*sizeof(double));
        int curr = 0;
        tiles.clear();
        for(int i = 0; i < num_tiles; i++){
    //        if(tile_num[i] < 0 ) continue;
            if(curr!= tiles.size()){cerr<< "number of tiles does not match"; cout << curr <<" "<<tiles.size()<<endl; exit(0);}
            tiles.push_back(Tile(grid_data, tile_corners+i*4, step));
//            cout << curr<<" "<<tiles.size()<<endl;;
            curr++;
            
        }
    }
    inline
    int max_tile(vector<panel> &vp){
        int max =0;
        for(vector<panel>::iterator itr = vp.begin(); itr!=vp.end(); itr++)
            if(itr->size() > max) max = itr->size();
        return max;
    }

    template<class T>
    Tile_Grid(apply_poni *poni, T *data, int nx, int ny, int mx, int my, double lambda);
   
    void save(const char filename[]){
        ofstream fo(filename);
        cout << "going to save"<<endl;
        fo.write(reinterpret_cast<char *>(size), sizeof(int)*2);
        fo.write(reinterpret_cast<char *>(step), sizeof(double)*2);
        fo.write(reinterpret_cast<char *>(offset), sizeof(double)*2);
        fo.write(reinterpret_cast<char *>(&num_grids),sizeof(int)); // number of valid tiles;
        fo.write(reinterpret_cast<char *>(&num_tiles),sizeof(int)); // number of valid tiles;
        int len = size[0] *size[1];
        cout << "len "<<len<<endl;
        fo.write(reinterpret_cast<char *>(tile_num), len*sizeof(int)); // point to the valid tiles;
        cout << "tile_num wrote"<<endl;
        fo.write(reinterpret_cast<char *>(tile_corners), num_tiles*4*sizeof(int));
        cout << "grid_num wrote"<<endl;
        fo.write(reinterpret_cast<char *>(grid_data), num_grids*4*sizeof(double));
        cout << "grid_num wrote"<<endl;

    }
    inline
    void divide(const double &x, const double &step, int &i, double &r)
    {
        r = x/step;
        i = int(r);
        r-=i;
    }
    inline
    bool out_of_range(int &i, int &j){
        return(i<0 || i >= size[0] || j < 0 || j >= size[1] );
    }
    double evaluate(double x, double y){
        double xr, yr;
        int xi, yi;
        x-= offset[0];
        y-= offset[1];
        double range = size[0]*step[0];
        while(x < 0) x += range;
        while(x> range) x -=range;
        divide(x, step[0], xi, xr);
        divide(y, step[1], yi, yr);
   //     cout << xi <<" "<<yi<<endl;
   //     cout << xr<<" "<<yr <<endl;
        if(out_of_range(xi, yi)) return 0;
        int n = tile_num[yi*size[0] + xi];
   //     cout <<"n "<<n<<endl;
        if (n < 0) return 0;
        if (n >= num_tiles){ cerr<< "n cannot exceed num_tils"; exit(0);}
        else return tiles[n].cal(xr, yr);
    }
protected:
    int size[2];
    double step[2];
    double offset[2];
    int num_tiles;
    int num_grids;
    int *tile_num;
    double *grid_data;
    vector<Tile> tiles;
    int *tile_corners;

};
#endif /* tile_grid_hpp */
