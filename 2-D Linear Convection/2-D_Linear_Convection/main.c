#include <stdio.h>
#include<omp.h>
#include <stdlib.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    double y_len = 2.0;
    int x_points = 351;
    int y_points = 351;
    double del_x = x_len/(x_points-1);
    double del_y = y_len/(y_points-1);

    double x[x_points], y[y_points];
    #pragma omp parallel
    {
        #pragma omp for nowait
        for(int i = 0; i < x_points; i++){
            x[i] = i * del_x;
        }
        #pragma omp for
        for(int j = 0; j < y_points; j++){
            y[j] = j * del_y;
        }
    }

//    printf("\n The x and y points are : \n \t");
//    for(int i = 0; i < x_points; i++){
//            for (int j = 0; j < y_points; j++){
//                printf("%f; %f \n \t", x[i], y[j]);
//            }
//    }

    double u[y_points][x_points];
    double u_new[y_points][x_points];

    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            u[i][j] = 1.0;
            u_new[i][j] = 1.0;
            if(x[i] > 0.5 && x[i] < 1.0 && y[i] > 0.5 && y[i] < 1.0){
                    u[i][j] = 2.0;
                    u_new[i][j] = 2.0;
            }
        }
    }
//
//    printf("\n The initial velocity is : \n \t");
//    for(int i = 0; i < y_points; i++){
//            for(int j = 0; j < x_points; j++){
//                printf("%f \t", u[i][j]);
//        }
//    }

    // Define the parameters
    int num_itrs = 700;
    double c = 1.0;
    double sigma = 0.2;
    double del_t = sigma * del_x;       // CFL criteria

    // Iterations
    double par_start_time = omp_get_wtime();

    #pragma omp parallel
    for (int itr = 0; itr < num_itrs; itr ++){

        #pragma omp for nowait
        for (int i = 1; i < y_points; i++){
            for (int j = 1; j < x_points; j++){

                    u_new[i][j] = u[i][j] - (c * del_t/del_x * (u[i][j] - u[i][j-1])) - (c * del_t/del_y * (u[i][j] - u[i-1][j]));
            }
        }

        #pragma omp for nowait
        for (int i = 0; i < y_points; i++){
            for (int j = 0; j < x_points; j++){

                    u[i][j] = u_new[i][j];
            }
        }

        // Setting the boundary value
        #pragma omp for nowait
        for (int i = 0; i < y_points; i++){
            u[i][0] = 1.0;
            u[i][x_points-1] = 1.0;
        }
        #pragma omp for
        for (int j = 0; j < x_points; j++){
            u[0][j] = 1.0;
            u[y_points-1][j] = 1.0;
        }

    }
    double par_end_time = omp_get_wtime();

//    printf("\n The final velocity is : \n \t");
//    for(int i = 0; i < y_points; i++){
//            for(int j = 0; j < x_points; j++){
//                printf("%f \t", u[i][j]);
//        }
//    }

    printf("\n Time taken for parallel computing is: %f \t", par_end_time - par_start_time);


    // ------------------------------------------------------------------------- //
    // Serial computing - for time comparison

    // Reinitializing u velocity
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            u[i][j] = 1.0;
            u_new[i][j] = 1.0;
            if(x[i] > 0.5 && x[i] < 1.0 && y[i] > 0.5 && y[i] < 1.0){
                    u[i][j] = 2.0;
                    u_new[i][j] = 2.0;
            }
        }
    }


    // Iteration
    double ser_start_time = omp_get_wtime();

    for (int itr = 0; itr < num_itrs; itr ++){

        for (int i = 1; i < y_points; i++){
            for (int j = 1; j < x_points; j++){

                    u_new[i][j] = u[i][j] - (c * del_t/del_x * (u[i][j] - u[i][j-1])) - (c * del_t/del_y * (u[i][j] - u[i-1][j]));
            }
        }

        for (int i = 0; i < y_points; i++){
            for (int j = 0; j < x_points; j++){

                    u[i][j] = u_new[i][j];
            }
        }

        // Setting the boundary value
        for (int i = 0; i < y_points; i++){
            u[i][0] = 1.0;
            u[i][x_points-1] = 1.0;
        }
        for (int j = 0; j < x_points; j++){
            u[0][j] = 1.0;
            u[y_points-1][j] = 1.0;
        }

    }
    double ser_end_time = omp_get_wtime();

//    printf("\n The final velocity is : \n \t");
//    for(int i = 0; i < y_points; i++){
//            for(int j = 0; j < x_points; j++){
//                printf("%f \t", u[i][j]);
//        }
//    }

    printf("\n Time taken for serial computing is: %f \t", ser_end_time - ser_start_time);

    printf("\n Speedup is : %f \t", (ser_end_time - ser_start_time)/(par_end_time - par_start_time));

    return 0;
}
