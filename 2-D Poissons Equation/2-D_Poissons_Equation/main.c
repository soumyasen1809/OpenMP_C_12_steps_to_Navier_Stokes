#include <stdio.h>
#include<omp.h>
#include<math.h>
#include <stdlib.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    double y_len = 1.0;
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

//    printf("\n <x,y> co-ordinates are: \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f; %f \n", x[j], y[i]);
//        }
//    }

    // Define the parameters
    int num_itr = 100;                              // Pseudo time iteration
    double p[y_points][x_points], p_new[y_points][x_points];
    double b[y_points][x_points];           // source term

    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            p[i][j] = 0.0;
            p_new[i][j] = 0.0;
        }
    }

    // Initialize source term - add spikes at 1/4th and 3/4th length
    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            b[i][j] = 0.0;
            if(i == abs(0.25*x_points) && j == abs(0.25*y_points)){
                b[i][j] = 100;
            }
            if(i == abs(0.75*x_points) && j == abs(0.75*y_points)){
                b[i][j] = -100;
            }
        }
    }

    // Iteration
    double par_start_time = omp_get_wtime();
    #pragma omp parallel
    {
        for(int it = 0; it < num_itr; it++){
            #pragma omp for nowait
            for(int i = 1; i < y_points-1; i++){
                for(int j = 1; j < x_points-1; j++){
                    p_new[i][j] = (del_y*del_y*(p[i][j+1] + p[i][j-1]) + del_x*del_x*(p[i+1][j] + p[i-1][j]) - del_x*del_x*del_y*del_y*b[i][j]) / (2*(del_x*del_x + del_y*del_y));
                }
            }

            // Boundary conditions
            #pragma omp for nowait
            for(int i = 0; i < y_points; i++){
                p_new[i][0] = 0;
                p_new[i][x_points-1] = 0;
            }
            #pragma omp for nowait
            for(int j = 0; j < x_points; j++){
                p_new[0][j] = 0;
                p_new[y_points-1][j] = 0;
            }


            // Assigning the new values to the previous values
            #pragma omp for
            for(int i = 0; i < y_points; i++){
                for(int j = 0; j < x_points; j++){
                    p[i][j] = p_new[i][j];
                }
            }

        }

    }
    double par_end_time = omp_get_wtime();

//    printf("\n The final pressure conditions are: \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f \n", p[i][j]);
//        }
//    }

    printf("\n Parallel execution time taken: %f \t", par_end_time - par_start_time);



    // Serial computation - for comparison

    // Defining the initial pressure
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            p[i][j] = 0.0;
            p_new[i][j] = 0.0;
        }
    }

    // Iterations
    double ser_start_time = omp_get_wtime();
    for(int it = 0; it < num_itr; it++){
        for(int i = 1; i < y_points-1; i++){
            for(int j = 1; j < x_points-1; j++){
                p_new[i][j] = (del_y*del_y*(p[i][j+1] + p[i][j-1]) + del_x*del_x*(p[i+1][j] + p[i-1][j]) - del_x*del_x*del_y*del_y*b[i][j]) / (2*(del_x*del_x + del_y*del_y));
            }
        }

        // Boundary conditions
        for(int i = 0; i < y_points; i++){
            p_new[i][0] = 0;
            p_new[i][x_points-1] = 0;
        }
         for(int j = 0; j < x_points; j++){
            p_new[0][j] = 0;
            p_new[y_points-1][j] = 0;
        }

        // Assigning the new values to the previous values
        for(int i = 0; i < y_points; i++){
            for(int j = 0; j < x_points; j++){
                p[i][j] = p_new[i][j];
            }
        }

    }
    double ser_end_time = omp_get_wtime();

    printf("\n Serial execution time taken: %f \t", ser_end_time - ser_start_time);

    printf("\n Speedup is : %f \t", (ser_end_time - ser_start_time)/(par_end_time - par_start_time));


    return 0;
}
