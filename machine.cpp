#include <iostream>
#include <cstdlib>
#include <random>
#include <ctime>
#include <cmath>

using namespace std;


//Declerations                                                                                                              
int gen(int size);


//Main class                                                                                                                
int main()
{
  //Define starting parameters                                                                                              
  int size = 5; //size of position matrix                                                                                   

  //Define channel matrix                                                                                                   
  double old_x_weight[size][size];
  double old_y_weight[size][size];
  double new_x_weight[size][size];
  double new_y_weight[size][size];


  for (int i=1; i<size; i++) {
      for (int j=1; j<size; j++) {
        old_x_weight[i][j]=j;
        old_y_weight[i][j]=j;
        new_x_weight[i][j]=j;
        new_y_weight[i][j]=j;
      }
  }

  //Generate position channels                                                                                              


  int iteration = size*5000;
  double precision =0.1;

  double x;
  double y;
  double guess_x;
  double guess_y;


  for (int chan_x=1; chan_x<size; chan_x++) {
    for (int chan_y=1; chan_y<size; chan_y++) {

      //begin                                                                                                               
      for (int j=1; j<iteration; j++) {
        x=0; y=0;
        guess_x=0; guess_y=0;
//Generate new weights                                                                                              
        int xychoice = gen(1);
        int wchoice = gen(size);
        int flip = gen(1);

        if (xychoice == 0) {
          if (flip==0) {
            new_x_weight[chan_y][wchoice] = old_x_weight[chan_y][wchoice]*2;
          }
          if (flip==1) {
            new_x_weight[chan_y][wchoice] = old_x_weight[chan_y][wchoice]*0.5;
          }
        }

        if (xychoice == 1) {
          if (flip==0) {
            new_y_weight[chan_x][wchoice] = old_y_weight[chan_x][wchoice]*2;
          }
          if (flip==1) {
            new_y_weight[chan_x][wchoice] = old_y_weight[chan_x][wchoice]*0.5;
          }
        }

//guess                                                                                                             
        for (int i=1; i<size; i++) {
          x += (1-exp(-old_x_weight[chan_y][i]))*chan_x;
          y += (1-exp(-old_y_weight[chan_x][i]))*chan_y;
          guess_x += (1-exp(-new_x_weight[chan_y][i]))*chan_x;
          guess_y += (1-exp(-new_y_weight[chan_x][i]))*chan_y;
        }

        if ((abs(chan_x-x)+abs(chan_y-y))>(abs(chan_x-guess_x)+abs(chan_y-guess_y))) {
          for (int i=1; i<size; i++) {
            for (int j=1; j<size; j++) {
              old_x_weight[i][j]=new_x_weight[i][j];
              old_y_weight[i][j]=new_y_weight[i][j];
            }
          }

          cout<<"New weights accepted, error = "<<abs(chan_x-guess_x)+abs(chan_y-guess_y)<<endl;

        }

        if ((abs(chan_x-guess_x)+abs(chan_y-guess_y))<precision) {
          cout<<"error<"<<precision<<" so finished"<<endl;
          break;}

      }
      cout<<chan_x<<","<<chan_y<<" error = "<<abs(chan_x-guess_x)+abs(chan_y-guess_y)<<endl;
    }
  }

cout<<"all combinations optimised";
 return 0;
}

int gen(int size)
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis(0,size);
  return dis(gen);
}
