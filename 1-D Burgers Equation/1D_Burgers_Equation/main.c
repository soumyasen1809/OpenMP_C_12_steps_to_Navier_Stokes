#include <stdio.h>
#include<omp.h>
#include <stdlib.h>
#include<math.h>

int main()
{
    double pi = 3.1416;

    // Define the domain
    double x_len = 2.0 * pi;
    int x_points = 1001;
    double del_x = x_len/(x_points-1);
    double x[x_points];

    #pragma omp parallel for
    for (int i = 0; i < x_points; i++){
            x[i] = i * del_x;
    }

    int num_itrs = 500;
    double nu = 0.07;
    double del_t = nu * del_x;

    double u[x_points], u_new[x_points];

    // Initial value of u
    #pragma omp parallel for
    for (int i = 0; i < x_points; i++){
        u[i] =  - 2.0 * nu * ( - (2.0 * x[i]) * exp( - (x[i] * x[i]) / (4.0 * nu)) / (4.0 * nu)  -  (2.0 * x[i]  -  4.0 * pi) * exp( - (x[i]  -  2.0 * pi) * (x[i]  -  2.0 * pi) / (4.0 * nu)) / (4.0 * nu)) / (exp( - (x[i]  -  2.0 * pi) * (x[i]  -  2.0 * pi) / (4.0 * nu)) + exp( - (x[i] * x[i]) / (4.0 * nu))) + 4.0;
        u_new[i] = u[i];
    }

//    printf("\n The initial u velocity is : \n \t");
//    for (int i = 0; i < x_points; i++){
//        printf("%f \t", u[i]);
//    }

    double par_start_time = omp_get_wtime();

    #pragma omp parallel
    for (int it = 0; it < num_itrs; it++){

        #pragma omp for nowait
        for (int i = 1; i < x_points-1; i++){
                u_new[i] = u[i] - u[i] * (del_t/del_x)* (u[i] - u[i-1]) + nu * (del_t/(del_x*del_x)) * (u[i+1] + u[i-1] - 2*u[i]);
        }

        #pragma omp for
        for (int i = 0; i < x_points; i++){
            u[i] = u_new[i];
        }

        u_new[0] = u[0] - u[0] * (del_t/del_x)* (u[0] - u[x_points-2]) + nu * (del_t/(del_x*del_x)) * (u[1] + u[x_points-2] - 2*u[0]);
        u_new[x_points-1] = u_new[0];
    }

    double par_end_time = omp_get_wtime();

//    printf("\n The final u velocity is : \n \t");
//    for (int i = 0; i < x_points; i++){
//        printf("%f \t", u[i]);
//    }

    printf("\n The time taken for parallel execution is : %f \t", par_end_time - par_start_time);


    // Serial Execution - for time comparison

    // Initial value of u
    for (int i = 0; i < x_points; i++){
        u[i] =  - 2.0 * nu * ( - (2.0 * x[i]) * exp( - (x[i] * x[i]) / (4.0 * nu)) / (4.0 * nu)  -  (2.0 * x[i]  -  4.0 * pi) * exp( - (x[i]  -  2.0 * pi) * (x[i]  -  2.0 * pi) / (4.0 * nu)) / (4.0 * nu)) / (exp( - (x[i]  -  2.0 * pi) * (x[i]  -  2.0 * pi) / (4.0 * nu)) + exp( - (x[i] * x[i]) / (4.0 * nu))) + 4.0;
        u_new[i] = u[i];
    }

    double ser_start_time = omp_get_wtime();

    for (int it = 0; it < num_itrs; it++){

        for (int i = 1; i < x_points-1; i++){
                u_new[i] = u[i] - u[i] * (del_t/del_x)* (u[i] - u[i-1]) + nu * (del_t/(del_x*del_x)) * (u[i+1] + u[i-1] - 2*u[i]);
        }

        for (int i = 0; i < x_points; i++){
            u[i] = u_new[i];
        }

        u_new[0] = u[0] - u[0] * (del_t/del_x)* (u[0] - u[x_points-2]) + nu * (del_t/(del_x*del_x)) * (u[1] + u[x_points-2] - 2*u[0]);
        u_new[x_points-1] = u_new[0];
    }

    double ser_end_time = omp_get_wtime();

    printf("\n The time taken for serial execution is : %f \t", ser_end_time - ser_start_time);

    printf("\n The speedup is : %f \t", (ser_end_time - ser_start_time)/ (par_end_time - par_start_time));

    return 0;
}
