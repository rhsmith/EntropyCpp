#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "math.h"
#include "tnt/tnt.h"
#include "tnt/jama_lu.h"

using namespace std;
using namespace boost;
using namespace TNT;

double LFM3DInterp(double x[], double y[], double z[], double bz[], double x0, double y0, double z0);
double kshell_tri_interp(double x[],double y[], double z[], double bz[], double x0, double y0, double z0, int KK);
int shellIndex(int i, int j, int k);
int kindex(int linearIndex);
int jindex(int linearIndex);
int iindex(int linearIndex);
double dist(double,double,double);
double sign(double x);
const double PI = M_PI;
const double dl = .3;


// Main Program Entry Point
int main()
{
  string positionsPath = "../EntropyCpp/positions.txt";
  const int N = 172250; //number of all points (i.e. including night side points)
  //const int N = 106795; //number of points excluding night side
  //const int nightCutOff = 19; //j value for the night side
  double xpos[N];
  double ypos[N];
  double zpos[N];
  ifstream posFile;

  //read in positions data
  posFile.open(positionsPath.c_str(), ios::in);
  if(!posFile.is_open()){cout << "Error reading file" << endl;exit(1);}
  for(int i = 0; i < N; i++)
  {
    posFile >> xpos[i];
    posFile >> ypos[i];
    posFile >> zpos[i];
  }

  for(int i=0; i < 10; i++)
  {
    cout<<xpos[i]<<", "<<ypos[i]<<", "<<zpos[i]<<endl;
  }
  cout << "Printing values" << endl;
  cout << xpos[shellIndex(4,4,0)] << endl;
  cout << ypos[shellIndex(4,4,0)] << endl;
  cout << zpos[shellIndex(4,4,0)] << endl;


  //read in hdf data (in text file format)
  string hdfPath = "../EntropyCpp/300";
  double pressure[N];
  double bx[N];
  double by[N];
  double bz[N];
  double unused;

  ifstream hdfFile;
  hdfFile.open(hdfPath.c_str(), ios::in);
  if(!hdfFile.is_open()){cout << "Error reading HDF file: " << hdfPath << endl; exit(1);}
  for(int i = 0; i < N; i++)
  {
    hdfFile >> pressure[i];
    hdfFile >> bx[i];
    hdfFile >> by[i];
    hdfFile >> bz[i];
    hdfFile >> unused;
    hdfFile >> unused;
    hdfFile >> unused;
  }

  cout << "done reading data" << endl;

  //create meshgrid
  int Nx = 26;
  int Ny = 28;
  double X[Nx][Ny];
  double Y[Nx][Ny];
//vector< vector <double> > X(26, vector<double> (26, 0));
//vector< vector <double> > Y(28, vector<double> (28, 0));
  for(int i=0; i < Nx; i++)
  {
    for(int j=0; j < Ny; j++)
    {
      X[i][j] = -2*j-5;
      Y[i][j] = 2*i-25;
    }
  }

  cout << X[4][2] << endl;

  //Main loop over Nx and Ny
  double volumetotal[Nx*Ny];
  double pressuretotal[Nx*Ny];
  double entrop[Nx*Ny];
  for(int j=0; j<Nx; j++)
  {
    for(int k=0; k<Ny; k++)
    {
      cout << j << ", " << k << endl;
      double xstart = X[j][k];
      double ystart = Y[j][k];
      double zstart = 0.1;

      double xtrack[10] = {xstart,0,0,0,0,0,0,0,0,0};
      double ytrack[10] = {ystart,0,0,0,0,0,0,0,0,0};
      double ztrack[10] = {zstart,0,0,0,0,0,0,0,0,0};

      double xtrack2[10] = {xstart,0,0,0,0,0,0,0,0,0};
      double ytrack2[10] = {ystart,0,0,0,0,0,0,0,0,0};
      double ztrack2[10] = {zstart,0,0,0,0,0,0,0,0,0};

      double bztest = LFM3DInterp(xpos,ypos,zpos,bz,xstart,ystart,zstart);

      if(bztest > 0)
      {
        double volume = 0;
        double press = LFM3DInterp(xpos,ypos,zpos,pressure,xstart,ystart,zstart);
        double bxtrack,bytrack,bztrack,btrack;
        double radi;
        for(int i =0; i< 1000; i++)
        {
          bxtrack=LFM3DInterp(xpos,ypos,zpos,bx,xtrack[i],ytrack[i],ztrack[i]);
          bytrack=LFM3DInterp(xpos,ypos,zpos,by,xtrack[i],ytrack[i],ztrack[i]);
          bztrack=LFM3DInterp(xpos,ypos,zpos,bz,xtrack[i],ytrack[i],ztrack[i]);
          btrack = sqrt(pow(bxtrack,2)+pow(bytrack,2)+pow(bztrack,2));

          xtrack[i+1]=bxtrack/btrack*dl+xtrack[i];
          ytrack[i+1]=bytrack/btrack*dl+ytrack[i];
          ztrack[i+1]=bztrack/btrack*dl+ztrack[i];

          volume = dl/btrack+volume;

          radi = sqrt(pow(xtrack[i+1],2)+pow(ytrack[i+1],2)+pow(ztrack[i+1],2));

          if(radi < 3.5) break;

          if(radi > 2*dist(xstart,ystart,zstart) || (sqrt(pow(ytrack[i],2)) > 45))
          {
            volume = 0;
            break;
          }

        }
        
        for(int i =0; i< 1000; i++)
        {
          bxtrack=LFM3DInterp(xpos,ypos,zpos,bx,xtrack2[i],ytrack2[i],ztrack2[i]);
          bytrack=LFM3DInterp(xpos,ypos,zpos,by,xtrack2[i],ytrack2[i],ztrack2[i]);
          bztrack=LFM3DInterp(xpos,ypos,zpos,bz,xtrack2[i],ytrack2[i],ztrack2[i]);
          btrack = sqrt(pow(bxtrack,2)+pow(bytrack,2)+pow(bztrack,2));

          xtrack2[i+1]=bxtrack/btrack*(-dl)+xtrack2[i];
          ytrack2[i+1]=bytrack/btrack*(-dl)+ytrack2[i];
          ztrack2[i+1]=bztrack/btrack*(-dl)+ztrack2[i];

          volume = dl/btrack+volume;

          radi = sqrt(pow(xtrack2[i+1],2)+pow(ytrack2[i+1],2)+pow(ztrack[i+1],2));

          if(radi < 3.5) break;

          if(radi > 2*dist(xstart,ystart,zstart) || (sqrt(pow(ytrack2[i],2)) > 45))
          {
            volume = 0;
            break;
          }

        }
        
        //volumetotal[Ny*j+k]=volume;
        //pressuretotal[Ny*j+k]=press;
        //entrop[Ny*j+k]=pow(press*volume,1.6666);
      }else
      {
        //volumetotal[Ny*j+k]=0;
        //pressuretotal[Ny*j+k]=0;
        //entrop[Ny*j+k]=0;
      }
    }
  }

  for(int j=1;j<10;j++)
  {
    //cout << volumetotal[j][1] << endl;
    //cout << pressuretotal[j][1] << endl;
    //cout << entrop[j][1] << endl;
  }

}


///// HELPER FUNCTIONS /////////////////////////////////////////////////////////////////////////////////////////////

double LFM3DInterp(double x[], double y[], double z[], double bz[], double x0, double y0, double z0)
{
  //cout << "Entering LFM3DInterp" << endl;
  int NI = 53;
  int NJ = 60;
  int NK = 65;
  double theta, theta1, theta2;
  int KK;

  for(int k=0; k<NK-1; k++)
  {
    theta1 = PI - sign(z[shellIndex(4,4,k)])*acos(-y[shellIndex(4,4,k)]/sqrt(pow(y[shellIndex(4,4,k)],2.) + pow(z[shellIndex(4,4,k)],2.)));
    theta2 = PI - sign(z[shellIndex(4,4,k+1)])*acos(-y[shellIndex(4,4,k+1)]/sqrt(pow(y[shellIndex(4,4,k+1)],2.) + pow(z[shellIndex(4,4,k+1)],2.)));
    
    if (z0==0){
      theta = PI/2 - PI/2*sign(y0);
    }else{
      theta = PI - sign(z0)*acos(-y0/sqrt(pow(y0,2)+pow(z0,2)));
    }

    KK=64;
    if((theta-theta1)*(theta-theta2)<0){
      KK = k;
      break;
    }
  }

  double rho0=sqrt(pow(y0,2)+pow(z0,2));
  double y1=rho0*cos(theta1);
  double z1=rho0*sin(theta1);
  double y2=rho0*cos(theta2);
  double z2=rho0*sin(theta2);
  double d1,d2,d;

  if(KK==65){
    d1=kshell_tri_interp(x,y,z,bz,x0,y1,z1,KK);
    d2=kshell_tri_interp(x,y,z,bz,x0,y2,z2,1);
  } else {
    d1=kshell_tri_interp(x,y,z,bz,x0,y1,z1,KK);
    d2=kshell_tri_interp(x,y,z,bz,x0,y2,z2,KK);
  }

  if (KK<64){
    d=(d2-d1)*(theta-theta1)/(theta2-theta1)+d1;
  }else{
    theta1=theta1-2*PI;
    d=(d2-d1)*(theta-theta1)/(theta2-theta1)+d1;
  }

  return d;

}

double kshell_tri_interp(double x[],double y[], double z[], double data[], double x0, double y1, double z1, int KK)
{
  //cout << "Entering kshel" << endl;
  int NI = 53;
  int NJ = 60;
  int NK = 65;
  double p[NI][NJ];
  double q[NI][NJ];
  double data1[NI][NJ];
  double d = 0.0;
  double p1,q1;

  for(int i = 0; i < NI; i++)
  {
    for(int j = 0; j < NJ; j++)
    {
      p[i][j] = x[shellIndex(i,j,KK)];
      q[i][j] = sqrt(pow(y[shellIndex(i,j,KK)],2.0)+pow(z[shellIndex(i,j,KK)],2.0));
      p1=x0;
      q1=sqrt(pow(y1,2.0)+pow(z1,2.0)); //Why is this line here? Seems like waste of computational time
      data1[i][j] = data[shellIndex(i,j,KK)];
    }
  }

  //Find out which triangle is (p1.q1) in
  //search through (i,j) pairs, each cell is divided into two triangles
  // 1 (i,j) (i+1,j),(i,j+1)
  // 2 (i+1,j+1) (i+1,j), (i,j+1)

  double s1[2];
  double s2[2];
  double s3[2];
  double s4[2];
  double xx1,yy1,ff1,xx2,yy2,ff2,xx3,yy3,ff3;

  for(int i = 0; i < NI-1; i++)
  {
    for(int j = 0; j < NJ-1; j++)
    {
      s1[0] = p[i][j]-p1;
      s1[1] = q[i][j]-q1;
      s2[0] = p[i+1][j]-p1;
      s2[1] = q[i+1][j]-q1;
      s3[0] = p[i+1][j+1]-p1;
      s3[1] = q[i+1][j+1]-q1;
      s4[0] = p[i][j+1]-p1;
      s4[1] = q[i][j+1]-q1;

      //Triangle 1, ANG(12)+ANG(24)+ANG(41)=2*pi
      double theta12, theta24, theta41;
      theta12=acos((s1[0]*s2[0]+s1[1]*s2[1])/sqrt((pow(s1[0],2)+pow(s1[1],2))*(pow(s2[0],2)+pow(s2[1],2))));
      theta24=acos((s2[0]*s4[0]+s2[1]*s4[1])/sqrt((pow(s2[0],2)+pow(s2[1],2))*(pow(s4[0],2)+pow(s4[1],2))));
      theta41=acos((s4[0]*s1[0]+s4[1]*s1[1])/sqrt((pow(s4[0],2)+pow(s4[1],2))*(pow(s1[0],2)+pow(s1[1],2))));

      if(abs(theta12+theta24+theta41-2.0*PI) < 0.001)
      {
        xx1=p[i][j];
        yy1=q[i][j];
        ff1=data1[i][j];
        xx2=p[i+1][j];
        yy2=q[i+1][j];
        ff2=data1[i+1][j];
        xx3=p[i][j+1];
        yy3=q[i][j+1];
        ff3=data1[i][j+1];
        break;
      }

      //Triangle 2, ANG(23)+ANG(34)+ANG(42)=2*pi
      double theta23, theta34, theta42;
      theta23=acos((s2[0]*s3[0]+s2[1]*s3[1])/sqrt((pow(s2[0],2)+pow(s2[1],2))*(pow(s3[0],2)+pow(s3[1],2))));
      theta34=acos((s3[0]*s4[0]+s3[1]*s4[1])/sqrt((pow(s3[0],2)+pow(s3[1],2))*(pow(s4[0],2)+pow(s4[1],2))));
      theta42=acos((s4[0]*s2[0]+s4[1]*s2[1])/sqrt((pow(s4[0],2)+pow(s4[1],2))*(pow(s2[0],2)+pow(s2[1],2))));


      if(abs(theta23+theta34+theta42-2.0*PI) < 0.001)
      {
        xx1=p[i+1][j+1];
        yy1=q[i+1][j+1];
        ff1=data1[i+1][j+1];
        xx2=p[i+1][j];
        yy2=q[i+1][j];
        ff2=data1[i+1][j];
        xx3=p[i][j+1];
        yy3=q[i][j+1];
        ff3=data1[i][j+1];
        break;
      }

    }
  }

  Array2D< double > arr1(3,3);
  Array2D< double > arr2(3,3);
  Array2D< double > arr3(3,3);
  Array2D< double > arr(3,3);

  double temp[3][3] = { {xx1, yy1, 1},
                       {xx2, yy2, 1},
                       {xx3, yy3, 1}};
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      arr[i][j] = temp[i][j];
    }
  }

  double temp1[3][3] = { {p1, q1, 1},
                        {xx2, yy2, 1},
                        {xx3, yy3, 1}};
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      arr1[i][j] = temp1[i][j];
    }
  }

  double temp2[3][3] = { {p1, q1, 1},
                        {xx1, yy1, 1},
                        {xx3, yy3, 1}};
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      arr2[i][j] = temp2[i][j];
    }
  }

  double temp3[3][3] = { {p1, q1, 1},
                        {xx1, yy1, 1},
                        {xx2, yy2, 1}};
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      arr3[i][j] = temp3[i][j];
    }
  }

  JAMA::LU< double > compute(arr);
  JAMA::LU< double > compute1(arr1);
  JAMA::LU< double > compute2(arr2);
  JAMA::LU< double > compute3(arr3);

  d = (ff1*compute1.det() - ff2*compute2.det() + ff3*compute3.det())/compute.det();


  return d;


}


double dist(double x, double y, double z)
{
  return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}

double dist(double x, double y)
{
  return sqrt(pow(x,2)+pow(y,2));
}









//// FUNCTIONS FOR SIMULATING RESHAPING OF DATA ///////////////////////////////////////////////////////////////////

int shellIndex(int i, int j, int k)
{
  //gets single index value for the point specified by i,j,k shell numbering
  //use this instead of "reshaping" the array for effeciency.
  int index = k*50*53+j*53+i;
  return index;
}

int kindex(int linearIndex)
{
  //gets the k shell index given an index in the linear array
  return (linearIndex/(50*53));
}


int jindex(int linearIndex)
{
  //gets the j shell index given an index in the linear array
  return (linearIndex%(50*53))/50;
}

int iindex(int linearIndex)
{
  //gets the i shell index given an index in the linear array
  return ((linearIndex%(50*53))%53);
}

double sign(double x)
{
  int sign = signbit(x);
  if(sign) return -1.0;
  else return 1.0;
}


