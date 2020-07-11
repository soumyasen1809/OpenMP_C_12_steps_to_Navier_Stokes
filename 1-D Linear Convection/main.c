#include<stdio.h>
#include<omp.h>
int main(){

  // Define the domain
  float x_len = 2.0;    // Length of the domain
  int x_points = 101;    // Number of points to consider
  float del_x = x_len/x_points;     // Length of an element
  float x[x_points];

  #pragma omp parallel for
  for (int i = 0; i < x_points; i++){
    x[i] = i * del_x;       // x co-ordinates
  }

  printf("\n The value of x \n");
  for (int i = 0; i < x_points; i++){
      printf("%f \t", x[i]);
  }

  // Define the parameters
  int t_itrs = 2500;  // number of time iterations
  float del_t = 0.001;
  float c = 1.0;  // speed of wave

  float u[x_points];                // Velocity at current time
  float u_new[x_points];       // Velocity at next time interval
  for (int i = 0; i < x_points; i++){
    if (x[i] > 0.5 && x[i] < 1.0){
      u[i] = 2.0;
      u_new[i] = 2.0;
    }
    else{
      u[i] = 1.0;
      u_new[i] = 1.0;
    }
  }

  printf("\n The initial value of u is \n");
  for (int i = 0; i < x_points; i++){
      printf("%f \t", u[i]);
  }

  // Loop iterations
  #pragma omp parallel
  for (int it = 0; it < t_itrs; it++){
    #pragma omp for nowait
    for (int i = 1; i < x_points; i++){
      u_new[i] = u[i] - (c*del_t/del_x)*(u[i] - u[i-1]);
    }

    #pragma omp for
    for (int i = 0; i < x_points; i++){
      u[i] = u_new[i];
    }
  }

  printf("\n The value of u at the end of the iterations \n");
  for (int i = 0; i < x_points; i++){
      printf("%f \t", u[i]);
  }

}
