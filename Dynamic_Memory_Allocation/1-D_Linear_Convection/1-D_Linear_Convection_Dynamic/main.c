#include <stdio.h>
#include<stdlib.h>
#include<omp.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    int x_points = 100001;
    double del_x = x_len/(x_points-1);
    double *x;
    x = (double*)malloc(x_points*sizeof(double));
    for(int i = 0; i < x_points; i++){
        *(x+i) = i * del_x;
    }
//    printf("\n Co-ordinates of x are: \n");
//    for(int i = 0; i < x_points; i++){
//        printf("%f \n", *(x+i));
//    }

    // Define the parameters
    double c = 1.0;
    double del_t = 0.025;
    int num_itr = 25;

    double *u, *u_new;
    u = (double*)malloc(x_points*sizeof(double));
    u_new = (double*)malloc(x_points*sizeof(double));
    for(int i = 0; i < x_points; i++){
        *(u+i) = 1.0;
        *(u_new+i) = 1.0;
        if(*(x+i) > 0.5 && *(x+i) < 1.0){
            *(u+i) = 2.0;
            *(u_new+i) = 2.0;
        }
    }

//    printf("\n Initial value of u is: \n");
//    for(int i = 0; i < x_points; i++){
//        printf("%f \n", *(u+i));
//    }

    // Iterations
    double ser_start_time = omp_get_wtime();
    for(int it = 0; it < num_itr; it++){
        for(int i = 1; i < x_points; i++){
            *(u_new+i) = *(u+i) - (c*del_t/del_x)*(*(u+i) - *(u+i-1));
        }
        for(int i = 0; i < x_points; i++){
            *(u+i) = *(u_new+i);
        }
    }
    double ser_end_time = omp_get_wtime();
//    printf("\n Final value of u is: \n");
//    for(int i = 0; i < x_points; i++){
//        printf("%f \n", *(u+i));
//    }
    printf("\n Time taken by serial execution is: %f \n", ser_end_time - ser_start_time);


    // Parallel execution

    // Defining the initial conditions
    #pragma omp parallel for
    for(int i = 0; i < x_points; i++){
        *(u+i) = 1.0;
        *(u_new+i) = 1.0;
        if(*(x+i) > 0.5 && *(x+i) < 1.0){
            *(u+i) = 2.0;
            *(u_new+i) = 2.0;
        }
    }

    // Iterations
    double par_start_time = omp_get_wtime();
    #pragma omp parallel
    {
        for(int it = 0; it < num_itr; it++){
            #pragma omp for
            for(int i = 1; i < x_points; i++){
                *(u_new+i) = *(u+i) - (c*del_t/del_x)*(*(u+i) - *(u+i-1));
            }
            #pragma omp for
            for(int i = 0; i < x_points; i++){
                *(u+i) = *(u_new+i);
            }
        }
    }
    double par_end_time = omp_get_wtime();

    printf("\n Time taken by parallel execution is: %f \n", par_end_time - par_start_time);

    printf("\n Speedup is: %f \n", (ser_end_time - ser_start_time)/ (par_end_time - par_start_time));

    free(x);
    free(u);
    free(u_new);


    return 0;
}
