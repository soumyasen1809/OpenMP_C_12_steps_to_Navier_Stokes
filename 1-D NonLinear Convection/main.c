#include<stdio.h>
#include<omp.h>
int main(){

  // Define the domain
  double x_len = 2.0;    // Length of the domain
  int x_points = 101;    // Number of points to consider
  double del_x = x_len/x_points;     // Length of an element
  double x[x_points];

  #pragma omp parallel for firstprivate(del_x)
  for (int i = 0; i < x_points; i++){
    x[i] = i * del_x;       // x co-ordinates
  }

  printf("\n The value of x \n");
  for (int i = 0; i < x_points; i++){
      printf("%f \t", x[i]);
  }

  // Define the parameters
  int t_itrs = 500;  // number of time iterations
  double del_t = 0.001;
  // double c = 1.0;  // speed of wave

  double u[x_points];                // Velocity at current time
  double u_new[x_points];       // Velocity at next time interval
  #pragma omp parallel for shared(x)
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
  double par_start_time = omp_get_wtime();
  #pragma omp parallel firstprivate(del_x, del_t)
  for (int it = 0; it < t_itrs; it++){

      #pragma omp for nowait
      for (int i = 1; i < x_points; i++){
          u_new[i] = u[i] - (u[i]*del_t/del_x)*(u[i] - u[i-1]);
      }

      #pragma omp for
      for (int i = 0; i < x_points; i++){
            u[i] = u_new[i];
      }
  }
  double par_end_time = omp_get_wtime();

  printf("\n The value of u at the end of the iterations \n");
  for (int i = 0; i < x_points; i++){
      printf("%f \t", u[i]);
  }

  printf("\n Time taken for parallel execution is: %f \t", par_end_time-par_start_time);


  // Serial execution of loop for comparing time

  // Initialize back to initial u value
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


  double ser_start_time = omp_get_wtime();

  for (int it = 0; it < t_itrs; it++){

      for (int i = 1; i < x_points; i++){
          u_new[i] = u[i] - (u[i]*del_t/del_x)*(u[i] - u[i-1]);
      }

      for (int i = 0; i < x_points; i++){
            u[i] = u_new[i];
      }
  }
  double ser_end_time = omp_get_wtime();

  printf("\n Time taken for serial execution is: %f \t", ser_end_time-ser_start_time);

  return 0;

}
