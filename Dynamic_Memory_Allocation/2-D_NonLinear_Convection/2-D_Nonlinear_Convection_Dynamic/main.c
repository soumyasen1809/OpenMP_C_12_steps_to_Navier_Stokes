#include <stdio.h>
#include<omp.h>
#include <stdlib.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    double y_len = 2.0;
    int x_points = 201;
    int y_points = 201;
    printf("\n The grid is of %d x %d points \n", x_points, y_points);
    double del_x = x_len/(x_points-1);
    double del_y = y_len/(y_points-1);

    double *x, *y;
    x = (double*)malloc(x_points*sizeof(double));
    y = (double*)malloc(x_points*sizeof(double));
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i = 0; i < x_points; i++){
            *(x+i) = i * del_x;
        }
        #pragma omp for
        for (int j = 0; j < y_points; j++){
            *(y+j) = j * del_y;
        }
    }

//        printf("\n The x and y points are : \n \t");
//        for(int i = 0; i < y_points; i++){
//                for (int j = 0; j < x_points; j++){
//                    printf("%f; %f \n \t", *(x+j), *(y+i));
//                }
//        }


    // Define the parameters
    int num_itrs = 50;
    double sigma = 0.2;
    double del_t = sigma * del_x;

    double *u, *u_new, *v, *v_new;
    u = (double *)malloc(y_points*x_points*sizeof(double));
    u_new = (double *)malloc(y_points*x_points*sizeof(double));
    v = (double *)malloc(y_points*x_points*sizeof(double));
    v_new = (double *)malloc(y_points*x_points*sizeof(double));
    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            *(u+i*x_points+j) = 1.0;
            *(v+i*x_points+j) = 1.0;
            *(u_new+i*x_points+j) = 1.0;
            *(v_new+i*x_points+j) = 1.0;

            if(*(x+j) > 0.5 && *(x+j) < 1.0 && *(y+i) > 0.5 && *(y+i) < 1.0){
                *(u+i*x_points+j) = 2.0;
                *(v+i*x_points+j) = 2.0;
                *(u_new+i*x_points+j) = 2.0;
                *(v_new+i*x_points+j) = 2.0;
            }
        }
    }

//    printf("\n Initial value of (u,v) at the grid points: \n \t");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f \n", *(u+i*x_points+j), *(v+i*x_points+j));
//        }
//    }

    // Iteration
    double par_start_time = omp_get_wtime();

    #pragma omp parallel
    for(int itr = 0; itr < num_itrs; itr++){
        #pragma omp for nowait
        for(int i = 1; i < y_points; i++){
            for(int j = 1; j < x_points; j++){
                *(u_new+i*x_points+j) = *(u+i*x_points+j) - *(u+i*x_points+j) * (del_t/del_x) * (*(u+i*x_points+j) - *(u+i*x_points+(j-1))) - *(v+i*x_points+j) * (del_t/del_x) * (*(u+i*x_points+j) - *(u+(i-1)*x_points+j));
                *(v_new+i*x_points+j) = *(v+i*x_points+j) - *(u+i*x_points+j) * (del_t/del_x) * (*(v+i*x_points+j) - *(v+i*x_points+(j-1))) - *(v+i*x_points+j) * (del_t/del_x) * (*(v+i*x_points+j) - *(v+(i-1)*x_points+j));
            }
        }

        // Boundary conditions apply
        #pragma omp for nowait
        for(int i = 0; i < x_points; i++){
            *(u_new+0*x_points+i) = 1.0;
            *(v_new+0*x_points+i) = 1.0;
            *(u_new+(y_points-1)*x_points+i) = 1.0;
            *(v_new+(y_points-1)*x_points+i) = 1.0;
        }
        #pragma omp for nowait
        for(int j = 0; j < y_points; j++){
            *(u_new+j*x_points+0) = 1.0;
            *(v_new+j*x_points+0) = 1.0;
            *(u_new+j*x_points+(x_points-1)) = 1.0;
            *(v_new+j*x_points+(x_points-1)) = 1.0;
        }

        // Assigning newer values to the older ones
        #pragma omp for
        for(int i = 0; i < y_points; i++){
            for(int j = 0; j < x_points; j++){
                *(u+i*x_points+j) = *(u_new+i*x_points+j);
                *(v+i*x_points+j) = *(v_new+i*x_points+j);
            }
        }
    }

    double par_end_time = omp_get_wtime();

//    printf("\n Final value of (u,v) at the grid points: \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f \n", *(u+i*x_points+j), *(v+i*x_points+j));
//        }
//    }

    printf("\n Time taken for parallel execution is: %f \t", par_end_time - par_start_time);



    // ----------------------------------------------------------------- //
    // Serial execution - for comparison
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            *(u+i*x_points+j) = 1.0;
            *(v+i*x_points+j) = 1.0;
            *(u_new+i*x_points+j) = 1.0;
            *(v_new+i*x_points+j) = 1.0;

            if(*(x+j) > 0.5 && *(x+j) < 1.0 && *(y+i) > 0.5 && *(y+i) < 1.0){
                *(u+i*x_points+j) = 2.0;
                *(v+i*x_points+j) = 2.0;
                *(u_new+i*x_points+j) = 2.0;
                *(v_new+i*x_points+j) = 2.0;
            }
        }
    }

    // Iteration
    double ser_start_time = omp_get_wtime();

    for(int itr = 0; itr < num_itrs; itr++){
        for(int i = 1; i < y_points; i++){
            for(int j = 1; j < x_points; j++){
                *(u_new+i*x_points+j) = *(u+i*x_points+j) - *(u+i*x_points+j) * (del_t/del_x) * (*(u+i*x_points+j) - *(u+i*x_points+(j-1))) - *(v+i*x_points+j) * (del_t/del_x) * (*(u+i*x_points+j) - *(u+(i-1)*x_points+j));
                *(v_new+i*x_points+j) = *(v+i*x_points+j) - *(u+i*x_points+j) * (del_t/del_x) * (*(v+i*x_points+j) - *(v+i*x_points+(j-1))) - *(v+i*x_points+j) * (del_t/del_x) * (*(v+i*x_points+j) - *(v+(i-1)*x_points+j));
            }
        }

        // Boundary conditions apply
        for(int i = 0; i < x_points; i++){
            *(u_new+0*x_points+i) = 1.0;
            *(v_new+0*x_points+i) = 1.0;
            *(u_new+(y_points-1)*x_points+i) = 1.0;
            *(v_new+(y_points-1)*x_points+i) = 1.0;
        }
        for(int j = 0; j < y_points; j++){
            *(u_new+j*x_points+0) = 1.0;
            *(v_new+j*x_points+0) = 1.0;
            *(u_new+j*x_points+(x_points-1)) = 1.0;
            *(v_new+j*x_points+(x_points-1)) = 1.0;
        }

        // Assigning newer values to the older ones
        for(int i = 0; i < y_points; i++){
            for(int j = 0; j < x_points; j++){
                *(u+i*x_points+j) = *(u_new+i*x_points+j);
                *(v+i*x_points+j) = *(v_new+i*x_points+j);
            }
        }
    }

    double ser_end_time = omp_get_wtime();

    printf("\n Time taken for serial execution is: %f \t", ser_end_time - ser_start_time);
    printf("\n Speedup is: %f \t", (ser_end_time - ser_start_time)/(par_end_time - par_start_time));

    free(u);
    free(u_new);
    free(v);
    free(v_new);

    return 0;
}
