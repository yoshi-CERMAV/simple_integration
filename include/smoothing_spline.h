//
//  smoothing_spline.hpp
//  
//
//  Created by Yoshiharu Nishiyama on 12/05/2024.
//

#ifndef smoothing_spline_hpp
#define smoothing_spline_hpp

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <Accelerate/Accelerate.h>
using namespace std;

char gpfmt[] = "/opt/homebrew/bin/gnuplot --persist";


class smoothing_spline{
public:
    smoothing_spline(int max_ndat){
        init(max_ndat, 16);
        gp = popen(gpfmt, "w");
    }
    int init(int max_ndat, int m){
        max_nlen = max_ndat + 27;
        jacobian = new double[max_nlen * m];
        A = new double[max_nlen * m];
        z = new double[max_nlen];
        z_temp = new double[max_nlen];
        z_cal = new double[max_nlen];
        z_diff = new double[max_nlen];
        weight = new double[max_nlen];
        lwork = -1;
        cout <<"initializing"<<endl;
        double worksize;
        dgels_("No transpose", &max_nlen, &num_para, &one, A, &max_nlen, z, &max_nlen, &worksize, &lwork, &info);
        lwork = worksize;
        work = new double[lwork];
        cout << "lwork: "<< lwork<<endl;
    }
    
    int mul(int n, double *x1, double *x2, double *x3){
        for(int i = 0; i < n; i++){
            x3[i] = x1[i] * x2[i];
        }
    }
    int mul(int n, double x1, double *x2, double *x3){
        for(int i = 0; i < n; i++){
            x3[i] = x1 * x2[i];
        }
    }
    int mul(int n, double x1, double *x2, double *x3, double *x4){
        for(int i = 0; i < n; i++){
            x4[i] = x1 * x2[i] * x3[i];
        }
    }

    int fill_constraint(double *J, double x[4], double y[4], double lambda){
        double *ptr = J;
        for(int j = 0; j < 4; j++){
            for(int i =0; i < 4; i++, ptr += nld){
                ptr[0] = x[i] * y[j] * lambda;
            }
        }
    }
    int calc_jacobian(double *xyz, int step, int n, double lambda){
        nld = n+27;
        vector<double *> J;
        J.clear();//J1.clear(); J2.clear();J3.clear();
        J.push_back (jacobian);
        double *J1 = jacobian + n;
        double *J2 = J1 + 9;
        double *J3 = J2 + 9;
        cout << "pointer OK"<<endl;
        for(int i = 1; i < 16; i++ ){
            J.push_back(J[i-1] + nld);
        }
        std::fill_n (J[0], n, 1); //1
        cblas_dcopy(n, xyz, step, J[1], 1);//x
        mul(n, J[1], J[1], J[2]); //x**2
        mul(n, J[1], J[2], J[3]); //x**3
        cblas_dcopy(n, xyz+1, step, J[4], 1);//y
        for(int i = 1; i < 4; i++){
            mul(n, J[4], J[i], J[i+4]); //yx, yxx, yxxx
        }
        for(int i = 0; i < 4; i++){
            mul(n, J[4], J[i+4], J[i+8]); //yx, yxx, yxxx
            mul(n, J[4], J[i+8], J[i+12]); //yx, yxx, yxxx
        }
        cout << "B <"<<endl;
        double B[] = {1, -1, 1, -1,  1, 0, 0, 0,  1, 1, 1, 1,
                      0, 1, -2,  3,  0, 1, 0, 0,  0, 1, 2, 3,
                      0, 0,  2, -6,  0, 0, 2, 0,  0, 0, 2, 6};

        int irow = 0;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++, irow++){
                fill_constraint(J1+irow, B+i*4,    B+24+j*4, lambda);
                fill_constraint(J2+irow, B+12+i*4, B+12+j*4, lambda);
                fill_constraint(J3+irow, B+24+i*4, B+j*4,    lambda);
            }
        }
  //      for(int i = 0; i < 27; i++){
  //          for(int j = 0; j < 16; j++) cout << J1[j*nld+i]<<" ";
  //          cout<<endl;
  //      }
    }
    
    int plot(double *x, double *data, int l)
    {
        fprintf(gp,"splot %10.6f + %10.6f * x + %10.6f *x**2 + %10.6f * x**3 +", x[0], x[1], x[2], x[3]);
        fprintf(gp,"%10.6f *y + %10.6f * x * y + %10.6f *x**2 * y  + %10.6f * x**3 * y +", x[4], x[5], x[6], x[7]);
        fprintf(gp,"%10.6f *y**2 + %10.6f * x *y**2 + %10.6f *x**2 *y**2 + %10.6f * x**3 *y**2 +", x[8], x[9], x[10], x[11]);
        fprintf(gp, "%10.6f *y**3+ %10.6f * x*y**3 + %10.6f *x**2*y**3 + %10.6f * x**3*y**3 ", x[12], x[13], x[14], x[15]);
        fputs(", \"-\" binary format = \"%double%double%double\" record = ",gp);
        fprintf(gp, "%d\n", l);
        fwrite(data, 24, l, gp);
 
        fflush(gp);
    }
    
    int fit(double *xyz, int step, int len_, double lambda , double *out){
        data_len = len_;
   //     lambda *= data_len;
        cout << "calculating jacobian"<< endl;
        calc_jacobian(xyz, step, data_len, lambda);
        cblas_dcopy(data_len, xyz+2, step, z, 1);//x
 //       for(int i = 0; i < 10; i++){
 //           for(int j = 0; j < 16; j++){
 //               cout << jacobian[j*nld+i]<< " ";
 //           }
 //           cout <<z[i]<<endl;
 //       }
        std::fill_n(z+data_len, 27, 0);
 //       cout << "filled "<<endl;
        nld = data_len + 27;
        cblas_dcopy(nld* num_para, jacobian, 1, A, 1);
        cblas_dcopy(nld, z, 1, z_temp, 1);
  //      cout << "copied"<<endl;
        cout << num_para << " "<<nld <<endl;
        dgels_("No transpose", &nld, &num_para, &one,A, &nld, z_temp, &nld, work, &lwork, &info);
 //       cout << "info "<<info<<endl;
 //       for(int i = 0; i < 16; i++) cout << z_temp[i]<<" ";
 //       cout << endl;
 //       cout <<"plotting"<<endl;
        cblas_dgemv(CblasColMajor, CblasNoTrans,
                    nld, num_para, 1., jacobian, nld,  z_temp, 1, 0., z_cal, 1 );
        cblas_daxpy(nld, -1, z, 1, z_cal, 1);
        int quarter = data_len/4;
  //      cout << "z_cal "<<endl;
  //    for(int i = 0; i < data_len; i++) cout << z_cal[i]<<endl;
        cblas_dcopy(nld, z_cal, 1, z_diff, 1);
        nth_element(z_cal, z_cal+quarter, z_cal+data_len);
        nth_element(z_cal+quarter, z_cal+data_len-quarter, z_cal+data_len);
        int l = data_len - 2*quarter;
        double sigma = cblas_ddot(l, z_cal+quarter, 1, z_cal + quarter, 1 );
        sigma /= l;
        sigma = sqrt(sigma);
        sigma /= 0.377;
        cout <<"sigma "<<sigma<<endl;
        double sigma3 = sigma * 3;
        fill_n(weight, data_len, 1.);
        int n_outlier = 0;
        for(int i = 0; i < data_len; i++){
            if (fabs(z_diff[i] ) >  sigma3) {
                weight[i] = 0;
                n_outlier++;
            }
        }
        if(n_outlier) cout <<"number of out_liers: "<< n_outlier<<endl;
        cycle();
    //    plot(z_temp,xyz, data_len);

        out[0] = z_temp[0];
        out[1] = z_temp[1];
        out[2] = z_temp[4];
        out[3] = z_temp[5];
        double s_back = cblas_ddot(data_len, z_diff, 1, z_diff, 1);
        return 0;
    }
    
    int cycle()
    {
        cblas_dcopy(nld * num_para, jacobian, 1, A, 1);
        cblas_dcopy(nld, z, 1, z_temp, 1);
        for(int i = 0; i < data_len; i++){
            double *ptr = A+i;
            for(int j = 0; j < num_para; j++, ptr+= nld){
                (*ptr)*=weight[i];
            }
            z_temp[i]*=weight[i];
        }
        // No transpose
        char *trans = "N";

        dgels_(trans, &nld, &num_para, &one,A, &nld, z_temp, &nld, work, &lwork, &info);
        cblas_dgemv(CblasColMajor, CblasNoTrans,
                    nld, num_para, 1., jacobian, nld,  z_temp, 1, 0., z_cal, 1 );
        cblas_daxpy(nld, -1, z, 1, z_cal, 1);
        
        return 0;
    }
protected:
    double *jacobian;
    double *A;
    double *z;
    double *z_temp;
    double *z_cal;
    double *z_diff;
    double *cof;
    double lambda;
    int max_nlen;
    FILE *gp;
    
    int data_len;
    int nld;

    int num_para = 16;
    int one = 1;
    double *work;
    int lwork;
    int info;
    double *weight;
};

#endif /* smoothing_spline_hpp */
