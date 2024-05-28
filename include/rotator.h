/*
 *  rotator.h
 *  ghemical
 *
 *  Created by Yoshiharu Nishiyama on 05/12/13.
 *  Copyright 2005 Cermav. All rights reserved.
 *
 */
 
#ifndef ROTATOR_H
#define ROTATOR_H 
#include <iostream>
#include <math.h>
#include <string.h>

static const double rad2deg = 180/M_PI;
using namespace std;
template <class T>
class rotator
{
protected:

   T *matrix;

public:
   rotator(){
     matrix = new T[9];
       cout << "constructing"<<endl;
	 init();
       cout << "initialized "<<endl;
   }
    const T operator()(int i)const {return matrix[i];}
    void combine(const rotator<T> &r){
        for(int i = 0; i !=3; i++){
            int i3 = i*3;
            T x = matrix[0+i3] * r(0) + matrix[1+i3] *r(3) + matrix[2+i3] * r(6);
            T y = matrix[0+i3] * r(1) + matrix[1+i3] *r(4) + matrix[2+i3] * r(7);
            T z = matrix[0+i3] * r(2) + matrix[1+i3] *r(5) + matrix[2+i3] * r(8);
            matrix[0+i3] = x;
            matrix[1+i3] = y;
            matrix[2+i3] = z;
        }
    }
   void swap(){ //swapping x and z coordinate
      T t= matrix[2]; matrix[2] = matrix[0]; matrix[0] = t;
      t= matrix[5]; matrix[5] = matrix[3]; matrix[3] = t;
      t= matrix[8]; matrix[8] = matrix[6]; matrix[6] = t;
   } 
   void dump(){
   cout << "matrix\n"<<matrix[0]<<" "<<matrix[1]<<" "<<matrix[2]<<endl;
   cout << matrix[3]<<" "<<matrix[4]<<" "<<matrix[5]<<endl;
   cout << matrix[6]<<" "<<matrix[7]<<" "<<matrix[8]<<endl;
   }
   void get_view(){
       cout <<matrix[0] <<", "<<matrix[1] <<", "<<matrix[2]<<",\\"<<endl;
       cout <<matrix[6] <<", "<<matrix[7] <<", "<<matrix[8]<<",\\"<<endl;
       cout <<matrix[3] <<", "<<matrix[4] <<", "<<matrix[5]<<",\\"<<endl;

   }
   
   void init()
   {
  	  for(int i = 0; i < 9; i++){ 
	    if(!(i%4)) matrix[i] = 1; 
	    else matrix[i] = 0;
     }
   }
   
   inline void rot(const T&x, T *a0, T *a1){      
      T c = cos(x);
	  T s = sin(x);
	  for(int i = 0; i < 3; i++){
	    T t2 = c* *(a1+i) - s* *(a0+i);
	    *(a0+i) =  c* *(a0+i) + s* *(a1+i);
	    *(a1+i) = t2;
	  }
  }

   void rotx(T x){
         rot(x, matrix+3, matrix+6);
   }
   void roty(T x){
         rot(x, matrix+6,  matrix);
   }
   void rotz(T x){
         rot(x, matrix, matrix+3);
   }
   void apply(T *x){
     T t0 = x[0]*matrix[0]+x[1]*matrix[1]+x[2]*matrix[2];
     T t1 = x[0]*matrix[3]+x[1]*matrix[4]+x[2]*matrix[5];
	 x[2] = x[0]*matrix[6]+x[1]*matrix[7]+x[2]*matrix[8];
	 x[0] = t0;
	 x[1] = t1;
   }
    void bring_on_x(T *x){
        init();
        double a0 = atan2(x[1], x[0]);
        rotz(a0);
        T temp[3];
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        double a1 = atan2(temp[2], temp[0]);
        roty(-a1);
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        cout << temp[0] <<" "<< temp[1] << " "<< temp[2] <<endl;
    }

    void bring_on_xy (const T *x, const T *y, T *a){
        init();
        a[0] = atan2(x[1], x[0]);
        cout << "x[0] x[1] x[2] a[0]"<< x[0] << " "<< x[1] << " "<< x[2]<<" "<<a[0]* rad2deg<<endl;
        cout << "y[0] y[1] y[2] a[0]"<< y[0] << " "<< y[1] << " "<< y[2]<<" "<<a[0]* rad2deg<<endl;
       rotz(a[0]);
        T temp[3];
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        a[1] = atan2(temp[2], temp[0]);
        roty(-a[1]);
        memcpy(temp, y, sizeof(T)*3);
        apply(temp);
        a[2] = atan2(temp[2], temp[1]);
        rotx(a[2]);
        memcpy(temp, y, sizeof(T)*3);
        apply(temp);
        cout << temp[0] <<" "<<temp[1] <<" "<< temp[2]<<endl;
        memcpy(temp, x, sizeof(T)*3);
        apply(temp);
        cout << temp[0] <<" "<<temp[1] <<" "<< temp[2]<<endl;
        cout << "angles "<<a[0] *rad2deg <<" "<< a[1] * rad2deg <<" "<< a[2] << rad2deg<<endl;
    }
};

template <class T>
class rotator_axe:public rotator<T>
{
protected:
   T *base_matrix;
public:
   rotator_axe():rotator<T>(){
 //    matrix = new T[9];
     base_matrix = new T[9];
 //    init();
   }
    rotator_axe(const T *u):rotator<T>(){
        base_matrix = new T[9];
        set_pole(u);
    }
    void set_pole(const T *u){
        T scale = 1./norm(u);
        base_matrix[0] = u[0]*scale;
        base_matrix[1] = u[1]*scale;
        base_matrix[2] = u[2]*scale;
        base_matrix[3] = base_matrix[0]*base_matrix[0];
        base_matrix[4] = base_matrix[1]*base_matrix[1];
        base_matrix[5] = base_matrix[2]*base_matrix[2];
        base_matrix[6] = base_matrix[0]*base_matrix[1];
        base_matrix[7] = base_matrix[1]*base_matrix[2];
        base_matrix[8] = base_matrix[2]*base_matrix[0];
    }
    T norm (const T *x){
        T sum = 0;
        for(int i = 0; i < 3; i++){
            sum+= x[i]*x[i];
        }
        return sqrt(sum);
    }
    void rotate(T x){
        T s = sin(x);
        T c = cos(x);
        T c1 = 1 - c;
        T u0u1c1 = base_matrix[6] * c1;
        T u1u2c1 = base_matrix[7] * c1;
        T u0u2c1 = base_matrix[8] * c1;
        T u0s = base_matrix[0] * s;
        T u1s = base_matrix[1] * s;
        T u2s = base_matrix[2] * s;
        rotator<T>::matrix[0] = c + base_matrix[3] * c1;
        rotator<T>::matrix[1] = u0u1c1 - u2s;
        rotator<T>::matrix[2] = u0u2c1 + u1s;
        rotator<T>::matrix[3] = u0u1c1 + u2s;
        rotator<T>::matrix[4] = c + base_matrix[4] * c1;
        rotator<T>::matrix[5] = u1u2c1 - u0s ;
        rotator<T>::matrix[6] = u0u2c1 - u1s;
        rotator<T>::matrix[7] = u1u2c1 + u0s;
        rotator<T>::matrix[8] = c + base_matrix[5] * c1;
    }
};





template<class T>
T angle(const T &x, const T &y)  // returns between -pi:pi
{
   return atan2(y, x);
}

template<class T>
void
getangle(const T b[6], T ang[3])
{
  rotator<T> rot;
  if(fabs(b[1])>1e-6||fabs(b[2]) > 1e-6){
  ang[0] = angle(b[2], b[1]);
  rot.rotx(-ang[0]);
  }
  else ang[0] = 0;
  T temp[3];

  memcpy(temp, b,sizeof(T)*3);
  rot.apply(temp);
  if(fabs(temp[2])>1e-6||fabs(temp[0]) > 1e-6){  
  ang[1] = angle(temp[2], temp[0]);
  rot.roty(-ang[1]);
  }
  else ang[1] = 0;

  memcpy(temp,&(b[3]), sizeof(T)*3);
  rot.apply(temp);
  if(fabs(temp[1])>1e-6||fabs(temp[0]) > 1e-6){  
  ang[2] = angle(temp[1], temp[0]);
  }else ang[2] = 0;
}


template<class T>
void
get_angle_zy(const T b0[3],  T ang[2])
{
  rotator<T> rot;
  if(fabs(b0[0])>1e-6||fabs(b0[1]) > 1e-6){
  ang[0] = angle(b0[0], b0[1]);
  rot.rotz(ang[0]);
  }
  else ang[0] = 0;
  T temp[3];

  memcpy(temp, b0,sizeof(T)*3);
//  cout <<"temp "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;  
  rot.apply(temp);
//  cout <<"temp "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;  
  if(fabs(temp[2])>1e-6||fabs(temp[0]) > 1e-6){  
  ang[1] = angle(temp[2], temp[0]);
  rot.roty(ang[1]);
  }
  else ang[1] = 0;
  memcpy(temp, b0,sizeof(T)*3);
  rot.apply(temp);
//  cout <<"temp "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;  
}

template <class T>
void 
get_angle_xy(const T b0[3], T ang[2])
{
   rotator <T> rot;
   ang[0] = atan2(b0[1], b0[2]);
   rot.rotx(ang[0]);
   T temp[3];
   memcpy(temp ,b0 ,sizeof(T)*3);
   rot.apply(temp);
   ang[1] = atan2(b0[0], b0[2]);

}





#endif

