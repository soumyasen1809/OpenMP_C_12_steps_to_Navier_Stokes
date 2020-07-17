#include <stdio.h>
#include<omp.h>
#include <stdlib.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    double y_len = 2.0;
    int x_points = 251;
    int y_points = 251;
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

//    printf("\n The domain coordinate is <x,y> \n \t");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f \n \t", x[j], y[i]);
//        }
//    }

    // Define the parameters
    int num_itrs = 40;
    double nu = 0.05;
    double sigma = 0.25;
    double del_t = sigma * del_x * del_y / nu;

    double u[y_points][x_points], u_new[y_points][x_points];
    double v[y_points][x_points], v_new[y_points][x_points];
    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            u[i][j] = 1.0;
            v[i][j] = 1.0;
            u_new[i][j] = 1.0;
            v_new[i][j] = 1.0;

            if(x[j] > 0.5 && x[j] < 1.0 && y[i] > 0.5 && y[i] < 1.0){
                u[i][j] = 2.0;
                v[i][j] = 2.0;
                u_new[i][j] = 2.0;
                v_new[i][j] = 2.0;
            }
        }
    }

//    printf("\n The initial velocity is <u,v> \n \t");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f \n \t", u[i][j], v[i][j]);
//        }
//    }

    // Iteration (parallel)
    double par_start_time = omp_get_wtime();

    #pragma omp parallel
    for(int itr = 0; itr < num_itrs; itr++){

        #pragma omp for nowait
        for(int i = 1; i < y_points-1; i++){
            for(int j = 1; j < x_points-1; j++){
                u_new[i][j] = u[i][j] + (nu*del_t/(del_x*del_x))*(u[i][j+1] + u[i][j-1] -2*u[i][j]) + (nu*del_t/(del_y*del_y))*(u[i+1][j] + u[i-1][j] -2*u[i][j]);
                v_new[i][j] = v[i][j] + (nu*del_t/(del_x*del_x))*(v[i][j+1] + v[i][j-1] -2*v[i][j]) + (nu*del_t/(del_y*del_y))*(v[i+1][j] + v[i-1][j] -2*v[i][j]);
            }
        }

        // Boundary conditions assign
        #pragma omp for nowait
        for(int i = 0; i < x_points; i++){
            u_new[0][i] = 1.0;
            v_new[0][i] = 1.0;
            u_new[x_points-1][i] = 1.0;
            v_new[x_points-1][i] = 1.0;
        }
        #pragma omp for nowait
        for(int j = 0; j < y_points; j++){
            u_new[j][0] = 1.0;
            v_new[j][0] = 1.0;
            u_new[j][y_points-1] = 1.0;
            v_new[j][y_points-1] = 1.0;
        }

        // Updating older values to newer ones
        #pragma omp for
        for(int i = 0; i < y_points; i++){
            for(int j = 0; j < x_points; j++){
                u[i][j] = u_new[i][j];
                v[i][j] = v_new[i][j];
            }
        }

    }

    double par_end_time = omp_get_wtime();

//    printf("\n The final velocity is <u,v> \n \t");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f \n \t", u[i][j], v[i][j]);
//        }
//    }

    printf("\n Time taken for parallel computing is: %f", par_end_time - par_start_time);


    // Serial computing - to compare time

    // Redefining velocities
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            u[i][j] = 1.0;
            v[i][j] = 1.0;
            u_new[i][j] = 1.0;
            v_new[i][j] = 1.0;

            if(x[j] > 0.5 && x[j] < 1.0 && y[i] > 0.5 && y[i] < 1.0){
                u[i][j] = 2.0;
                v[i][j] = 2.0;
                u_new[i][j] = 2.0;
                v_new[i][j] = 2.0;
            }
        }
    }


    // Iteration (parallel)
    double ser_start_time = omp_get_wtime();

    for(int itr = 0; itr < num_itrs; itr++){

        for(int i = 1; i < y_points-1; i++){
            for(int j = 1; j < x_points-1; j++){
                u_new[i][j] = u[i][j] + (nu*del_t/(del_x*del_x))*(u[i][j+1] + u[i][j-1] -2*u[i][j]) + (nu*del_t/(del_y*del_y))*(u[i+1][j] + u[i-1][j] -2*u[i][j]);
                v_new[i][j] = v[i][j] + (nu*del_t/(del_x*del_x))*(v[i][j+1] + v[i][j-1] -2*v[i][j]) + (nu*del_t/(del_y*del_y))*(v[i+1][j] + v[i-1][j] -2*v[i][j]);
            }
        }

        // Boundary conditions assign
        for(int i = 0; i < x_points; i++){
            u_new[0][i] = 1.0;
            v_new[0][i] = 1.0;
            u_new[x_points-1][i] = 1.0;
            v_new[x_points-1][i] = 1.0;
        }
        for(int j = 0; j < y_points; j++){
            u_new[j][0] = 1.0;
            v_new[j][0] = 1.0;
            u_new[j][y_points-1] = 1.0;
            v_new[j][y_points-1] = 1.0;
        }

        // Updating older values to newer ones
        for(int i = 0; i < y_points; i++){
            for(int j = 0; j < x_points; j++){
                u[i][j] = u_new[i][j];
                v[i][j] = v_new[i][j];
            }
        }

    }

    double ser_end_time = omp_get_wtime();

    printf("\n Time taken for serial computing is: %f", ser_end_time - ser_start_time);
    printf("\n Speedup is \t : %f", (ser_end_time - ser_start_time)/(par_end_time - par_start_time));

    return 0;
}
