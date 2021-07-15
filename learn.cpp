#include <iostream>
#include <cstdlib>
#include <random>
#include <ctime>
#include <cmath>

using namespace std;


//Declerations                                                                                                              
int generate(int size);

//Main class                                                                                                                
int main()
{
  //Define starting parameters                                                                                              
  int size = 16; //size of position matrix                                                                                  

  //Define channel matrix                                                                                                   
  double x_weight[size][size];
  double y_weight[size][size];

  for (int i=1; i<size; i++) {
      for (int j=1; j<size; j++) {
        x_weight[i][j]=j;
        y_weight[i][j]=j;
      }
  }

  //Generate position channels                                                                                              
  int chan_x = generate(size);
  int chan_y = generate(size);

  int iteration = 1000;
  double precision =0.01;

  double x;
  double y;

  //begin                                                                                                                   
  for (int j=1; j<iteration; j++) {
    cout<<"Iteration "<<j<<endl;
    x=0; y=0;

    //train                                                                                                                 
    for (int i=1; i<size; i++) {
          x_weight[chan_y][i] = x_weight[chan_y][i]/2;
          y_weight[chan_x][i] = y_weight[chan_x][i]/2;
    }

    x_weight[chan_y][chan_x] = x_weight[chan_y][chan_x]*4;
    y_weight[chan_x][chan_y] = y_weight[chan_x][chan_y]*4;

    //guess                                                                                                                 
    for (int i=1; i<size; i++) {
          x += (1-exp(-x_weight[chan_y][i]))*chan_x;
          y += (1-exp(-y_weight[chan_x][i]))*chan_y;
    }



    if (abs(chan_x-x)<precision&&abs(chan_y-y)<precision) {

      cout<<chan_x<<","<<chan_y<<endl;
      cout<<x<<","<<y<<endl;
      break;
    }

    cout<<chan_x<<","<<chan_y<<endl;
    cout<<x<<","<<y<<endl;

  }



  return 0;
}

int generate(int size)
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis(1,size-1);
  return dis(gen);
}

