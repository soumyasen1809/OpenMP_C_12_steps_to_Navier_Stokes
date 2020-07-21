#include <stdio.h>
#include<omp.h>
#include<math.h>
#include <stdlib.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    double y_len = 1.0;
    int x_points = 301;
    int y_points = 301;
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

        printf("\n <x,y> co-ordinates are: \n");
        for(int i = 0; i < y_points; i++){
            for(int j = 0; j < x_points; j++){
                printf("%f; %f \n", x[j], y[i]);
            }
        }

    // Define the parameters
    double p[y_points][x_points], p_new[y_points][x_points];
    double l1norm = 1.0;
    double l1norm_limit = 0.0001;
    double sum_num, sum_den;

    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            p[i][j] = 0.0;
            p_new[i][j] = 0.0;
        }
    }

    // Add boundary conditions
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < y_points; i++){
            p[i][0] = 0;                                        // p = 0 at x = 0
            p[i][x_points-1] = y[i];                      // p = y at x = 2
            p_new[i][0] = 0;
            p_new[i][x_points-1] = y[i];
        }
        #pragma omp for
        for(int j = 0; j < x_points; j++){
            p[0][j] = p[1][j];                                      // dp/dy = 0 at y = 0
            p[y_points-1][j] = p[y_points-2][j];       // dp/dy = 0 at y = 1
            p_new[0][j] = p_new[1][j];
            p_new[y_points-1][j] = p_new[y_points-2][j];
        }
    }

//    printf("\n The initial pressure condition is: \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f \n", p[i][j]);
//        }
//    }

    // Iteration
    double par_start_time = omp_get_wtime();
    #pragma omp parallel shared(sum_num, sum_den)
    {
        while(l1norm > l1norm_limit){
            #pragma omp for nowait
            for(int i = 1; i < y_points-1; i++){
                for(int j = 1; j < x_points-1; j++){
                    p_new[i][j] = (del_y*del_y*(p[i][j+1] + p[i][j-1]) + del_x*del_x*(p[i+1][j] + p[i-1][j])) / (2*(del_x*del_x + del_y*del_y));
                }
            }

            // Boundary conditions
            #pragma omp for
            for(int i = 0; i < y_points; i++){
                // p[i][0] = 0;
                // p[i][x_points-1] = y[i];
                p_new[i][0] = 0;
                p_new[i][x_points-1] = y[i];
            }
            #pragma omp for
            for(int j = 0; j < x_points; j++){
                // p[0][j] = p[1][j];
                // p[y_points-1][j] = p[y_points-2][j];
                p_new[0][j] = p_new[1][j];
                p_new[y_points-1][j] = p_new[y_points-2][j];
            }

            // Find l1 norm
            sum_num = 0;
            sum_den = 0;
            #pragma omp for reduction(+: sum_num, sum_den)
            for(int i = 0; i < y_points; i++){
                for(int j = 0; j < x_points; j++){
                    sum_num += fabs(p_new[i][j]) - fabs(p[i][j]) ;
                    sum_den += fabs(p[i][j]);
                }
            }
            l1norm = sum_num/sum_den;


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

    printf("\n The final pressure conditions are: \n");
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            printf("%f \n", p[i][j]);
        }
    }

    printf("\n Parallel execution time taken: %f \t", par_end_time - par_start_time);



    // Serial computation - for comparison
    l1norm = 1.0;
    // Defining the initial pressure
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            p[i][j] = 0.0;
            p_new[i][j] = 0.0;
        }
    }
    // Defining B.C.s for initial pressure
    for(int i = 0; i < y_points; i++){
        p[i][0] = 0;
        p[i][x_points-1] = y[i];
        p_new[i][0] = 0;
        p_new[i][x_points-1] = y[i];
    }
    for(int j = 0; j < x_points; j++){
        p[0][j] = p[1][j];
        p[y_points-1][j] = p[y_points-2][j];
        p_new[0][j] = p_new[1][j];
        p_new[y_points-1][j] = p_new[y_points-2][j];
    }

    // Iterations
    double ser_start_time = omp_get_wtime();

        while(l1norm > l1norm_limit){
            for(int i = 1; i < y_points-1; i++){
                for(int j = 1; j < x_points-1; j++){
                    p_new[i][j] = (del_y*del_y*(p[i][j+1] + p[i][j-1]) + del_x*del_x*(p[i+1][j] + p[i-1][j])) / (2*(del_x*del_x + del_y*del_y));
                }
            }

            // Boundary conditions
            for(int i = 0; i < y_points; i++){
                // p[i][0] = 0;
                // p[i][x_points-1] = y[i];
                p_new[i][0] = 0;
                p_new[i][x_points-1] = y[i];
            }
            for(int j = 0; j < x_points; j++){
                // p[0][j] = p[1][j];
                // p[y_points-1][j] = p[y_points-2][j];
                p_new[0][j] = p_new[1][j];
                p_new[y_points-1][j] = p_new[y_points-2][j];
            }

            // Find l1 norm
            sum_num = 0;
            sum_den = 0;
            for(int i = 0; i < y_points; i++){
                for(int j = 0; j < x_points; j++){
                    sum_num += fabs(p_new[i][j]) - fabs(p[i][j]) ;
                    sum_den += fabs(p[i][j]);
                }
            }
            l1norm = sum_num/sum_den;


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
