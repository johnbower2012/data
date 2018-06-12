#include<iostream>
#include<random>
#include<cmath>

int main(int argc, char* argv[]){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0.75,1.5);

  double mean=0.0,std=0.0,x=0.0;
  int steps=10000;
  for(int i=0;i<steps;i++){
    x = distribution(generator);
    mean += x;
    std += (x-mean/(double) (i+1))*(x-mean/(double) (i+1));
  }
  mean /= (double) steps;
  std /= (double) steps;
  std -= mean*mean;
  printf("%f +/- %f\n",mean,std);

  return 0;
}
