#include <stdio.h>
#include<omp.h>
#include<math.h>
#include <stdlib.h>

int main()
{
    // Define the domain
    double x_len = 2.0;
    double y_len = 2.0;
    int x_points = 151;
    int y_points = 151;
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

//    printf("\n The <x,y> co-ordinates are \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f \n", x[j], y[i]);
//        }
//    }

    // Define parameters
    int num_time_itrs = 100;        // Number of time iterations
    int num_pres_itrs = 50;         // Number of pseudo-time iterations for pressure calculation (Poisson's equations)
    double rho = 1.0;
    double nu = 0.1;
    double del_t = 0.001;

    double u[y_points][x_points], v[y_points][x_points], p[y_points][x_points];
    double u_new[y_points][x_points], v_new[y_points][x_points], p_new[y_points][x_points];
    #pragma omp parallel for
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            u_new[i][j] = 0.0;
            v_new[i][j] = 0.0;
            p_new[i][j] = 0.0;
        }
    }

    // Iterations
    double par_start_time = omp_get_wtime();
    #pragma omp parallel
    {
        for(int it1 = 0; it1 < num_time_itrs; it1++){
            // Velocity calculations
            #pragma omp for nowait
            for(int i = 1; i < y_points-1; i++){
                for(int j = 1; j < x_points-1; j++){
                    u_new[i][j] = u[i][j] - u[i][j]*(del_t/del_x)*(u[i][j] - u[i][j-1]) - v[i][j]*(del_t/del_y)*(u[i][j] - u[i-1][j]) - (del_t/(2*rho*del_x))*(p[i][j+1] - p[i][j-1]) + nu*( ((del_t/(del_x*del_x)) * (u[i][j+1] + u[i][j-1] - 2*u[i][j])) + ((del_t/(del_y*del_y))*(u[i+1][j] + u[i-1][j] - 2*u[i][j])) );
                    v_new[i][j] = v[i][j] - u[i][j]*(del_t/del_x)*(v[i][j] - v[i][j-1]) - v[i][j]*(del_t/del_y)*(v[i][j] - v[i-1][j]) - (del_t/(2*rho*del_y))*(p[i+1][j] - p[i-1][j]) + nu*( ((del_t/(del_x*del_x)) * (v[i][j+1] + v[i][j-1] - 2*v[i][j])) + ((del_t/(del_y*del_y))*(v[i+1][j] + v[i-1][j] - 2*v[i][j])) );
                }
            }

            // Assign boundary conditions in u and v
            #pragma omp for nowait
            for(int i = 0; i < y_points; i++){
                u_new[i][0] = 0;                                          // u = 0 at x = 0
                u_new[i][x_points-1] = 0;                          // u = 0 at x = 2
                v_new[i][0] = 0;                                          // v = 0 at x = 0
                v_new[i][x_points-1] = 0;                          // v = 0 at x = 2
            }
            #pragma omp for
            for(int j = 0; j < x_points; j++){
                u_new[0][j] = 0;                                          // u = 0 at y = 0
                u_new[y_points-1][j] = 1.0;                       // u = 1 at y = 2
                v_new[0][j] = 0;                                          // v = 0 at y = 0
                v_new[y_points-1][j] = 0;                          // v = 0 at y = 2
            }

            // Assign new velocity values to old ones
            #pragma omp for
            for(int i = 0; i < y_points; i++){
                for(int j = 0; j < x_points; j++){
                    u[i][j] = u_new[i][j];
                    v[i][j] = v_new[i][j];
                }
            }
            // End velocity calculations


            // Pressure calculations
            // #pragma omp for
            for(int it2 = 0; it2 < num_pres_itrs; it2++){
                #pragma omp for nowait
                for(int i = 1; i < y_points-1; i++){
                    for(int j = 1; j < x_points-1; j++){
                        p_new[i][j] = ( ( (del_y*del_y*(p[i][j+1] + p[i][j-1])) + (del_x*del_x*(p[i+1][j] + p[i-1][j])) ) / (2 * ((del_x*del_x)+(del_y*del_y))) ) - ( ((rho*del_x*del_x*del_y*del_y) / (2 * ((del_x*del_x)+(del_y*del_y)))) * ( ( (1/del_t)*( ( (u[i][j+1] - u[i][j-1])/(2*del_x) ) + ( (v[i+1][j] - v[i-1][j])/(2*del_y) ) ) ) - ( ( (u[i][j+1] - u[i][j-1])/(2*del_x) ) * ( (u[i][j+1] - u[i][j-1])/(2*del_x) ) ) - 2.0*( ((u[i+1][j] - u[i-1][j])/(2*del_y)) * ((v[i][j+1] - v[i][j-1])/(2*del_x)) ) - (( (v[i+1][j] - v[i-1][j])/(2*del_y) ) * ( (v[i+1][j] - v[i-1][j])/(2*del_y) )) ));
                    }
                }

                // Boundary conditions in p
                #pragma omp for
                for(int i = 0; i < y_points; i++){
                    p_new[i][0] = p_new[i][1];                                          // dp/dx = 0 at x = 0
                    p_new[i][x_points-1] = p_new[i][x_points-2];           // dp/dx = 0 at x = 2
                }
                #pragma omp for
                for(int j = 0; j < x_points; j++){
                    p_new[0][j] = p_new[1][j];                                          // dp/dy = 0 at y = 0
                    p_new[y_points-1][j] = 0.0;                                       // p = 0 at y = 2
                }

                // Assign new value of p to old p
                #pragma omp for
                for(int i = 0; i < y_points; i++){
                    for(int j = 0; j < x_points; j++){
                        p[i][j] = p_new[i][j];
                    }
                }

            }
            // End pressure calculation

        }
    }
    double par_end_time = omp_get_wtime();

//    printf("\n New velocity and pressure conditions <u;v;p> \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f ; %f \n", u[i][j], v[i][j], p[i][j]);
//        }
//    }

    printf("\n Parallel execution time taken is : %f \n", par_end_time - par_start_time);



    // Serial execution - for comparison
    // Initializing all velocities and pressure conditions to 0
    for(int i = 0; i < y_points; i++){
        for(int j = 0; j < x_points; j++){
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            u_new[i][j] = 0.0;
            v_new[i][j] = 0.0;
            p_new[i][j] = 0.0;
        }
    }

    // Iterations
    double ser_start_time = omp_get_wtime();

        for(int it1 = 0; it1 < num_time_itrs; it1++){
            // Velocity calculations
            for(int i = 1; i < y_points-1; i++){
                for(int j = 1; j < x_points-1; j++){
                    u_new[i][j] = u[i][j] - u[i][j]*(del_t/del_x)*(u[i][j] - u[i][j-1]) - v[i][j]*(del_t/del_y)*(u[i][j] - u[i-1][j]) - (del_t/(2*rho*del_x))*(p[i][j+1] - p[i][j-1]) + nu*( ((del_t/(del_x*del_x)) * (u[i][j+1] + u[i][j-1] - 2*u[i][j])) + ((del_t/(del_y*del_y))*(u[i+1][j] + u[i-1][j] - 2*u[i][j])) );
                    v_new[i][j] = v[i][j] - u[i][j]*(del_t/del_x)*(v[i][j] - v[i][j-1]) - v[i][j]*(del_t/del_y)*(v[i][j] - v[i-1][j]) - (del_t/(2*rho*del_y))*(p[i+1][j] - p[i-1][j]) + nu*( ((del_t/(del_x*del_x)) * (v[i][j+1] + v[i][j-1] - 2*v[i][j])) + ((del_t/(del_y*del_y))*(v[i+1][j] + v[i-1][j] - 2*v[i][j])) );
                }
            }

            // Assign boundary conditions in u and v
            for(int i = 0; i < y_points; i++){
                u_new[i][0] = 0;                                          // u = 0 at x = 0
                u_new[i][x_points-1] = 0;                          // u = 0 at x = 2
                v_new[i][0] = 0;                                          // v = 0 at x = 0
                v_new[i][x_points-1] = 0;                          // v = 0 at x = 2
            }
            for(int j = 0; j < x_points; j++){
                u_new[0][j] = 0;                                          // u = 0 at y = 0
                u_new[y_points-1][j] = 1.0;                       // u = 1 at y = 2
                v_new[0][j] = 0;                                          // v = 0 at y = 0
                v_new[y_points-1][j] = 0;                          // v = 0 at y = 2
            }

            // Assign new velocity values to old ones
            for(int i = 0; i < y_points; i++){
                for(int j = 0; j < x_points; j++){
                    u[i][j] = u_new[i][j];
                    v[i][j] = v_new[i][j];
                }
            }
            // End velocity calculations


            // Pressure calculations
            for(int it2 = 0; it2 < num_pres_itrs; it2++){
                for(int i = 1; i < y_points-1; i++){
                    for(int j = 1; j < x_points-1; j++){
                        p_new[i][j] = ( ( (del_y*del_y*(p[i][j+1] + p[i][j-1])) + (del_x*del_x*(p[i+1][j] + p[i-1][j])) ) / (2 * ((del_x*del_x)+(del_y*del_y))) ) - ( ((rho*del_x*del_x*del_y*del_y) / (2 * ((del_x*del_x)+(del_y*del_y)))) * ( ( (1/del_t)*( ( (u[i][j+1] - u[i][j-1])/(2*del_x) ) + ( (v[i+1][j] - v[i-1][j])/(2*del_y) ) ) ) - ( ( (u[i][j+1] - u[i][j-1])/(2*del_x) ) * ( (u[i][j+1] - u[i][j-1])/(2*del_x) ) ) - 2.0*( ((u[i+1][j] - u[i-1][j])/(2*del_y)) * ((v[i][j+1] - v[i][j-1])/(2*del_x)) ) - (( (v[i+1][j] - v[i-1][j])/(2*del_y) ) * ( (v[i+1][j] - v[i-1][j])/(2*del_y) )) ));
                    }
                }

                // Boundary conditions in p
                for(int i = 0; i < y_points; i++){
                    p_new[i][0] = p_new[i][1];                                          // dp/dx = 0 at x = 0
                    p_new[i][x_points-1] = p_new[i][x_points-2];           // dp/dx = 0 at x = 2
                }
                for(int j = 0; j < x_points; j++){
                    p_new[0][j] = p_new[1][j];                                          // dp/dy = 0 at y = 0
                    p_new[y_points-1][j] = 0.0;                                       // p = 0 at y = 2
                }

                // Assign new value of p to old p
                for(int i = 0; i < y_points; i++){
                    for(int j = 0; j < x_points; j++){
                        p[i][j] = p_new[i][j];
                    }
                }

            }
            // End pressure calculation

        }
    double ser_end_time = omp_get_wtime();

//    printf("\n New velocity and pressure conditions <u;v;p> \n");
//    for(int i = 0; i < y_points; i++){
//        for(int j = 0; j < x_points; j++){
//            printf("%f ; %f ; %f \n", u[i][j], v[i][j], p[i][j]);
//        }
//    }

    printf("\n Serial execution time taken is : %f \n", ser_end_time - ser_start_time);
    printf("\n Speedup is : %f \n", (ser_end_time - ser_start_time)/(par_end_time - par_start_time));

    return 0;
}
