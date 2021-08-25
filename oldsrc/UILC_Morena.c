#include "../include/UILC.h"

double UILC_f_Morena_Linear(
    const double x, 
    void* p
)
{
    UILC_fparams_linear * params = (UILC_fparams_linear *)p;
    const unsigned int m =(params->m);
    const unsigned int n =(params->n);
    double y =0.00;
    for(int i=1; i<n+1;i++)
    {
        y += (1-(m+3)*gsl_pow_2( (n+1-2*(i))*(x/2) ) ) / pow(( 1+ gsl_pow_2((x/2)*(n+1-2*(i)))),(m/2+3));
    }
    return y;
}
double UILC_f_Morena_SquareGrid(
    const double x, 
    void * p
)
{
    UILC_fparams_Rectangle * params = (UILC_fparams_Rectangle *)p;
    const unsigned int m = (params->m);
    const unsigned int N = (params->N);
    const unsigned int M = (params->M);

    double y =0.00;
    for(int i =1; i < N+1; i++ )
    {
        for(int j=1; j< M+1; j++)
        {
            y+= (1-((m+3)*gsl_pow_2(N+1-2*i)-gsl_pow_2(M+1-2*j))*gsl_pow_2(x/2))/pow((gsl_pow_2(N+1-2*i)+gsl_pow_2(M+1-2*j))*gsl_pow_2(x/2)+1,(m/2+3));
        }    
    }
    return y;
}

double UILC_f_Morena_getdm_Linear(
    const UILC_Lamber_LED l, 
    const int led_n, 
    const unsigned int itetnum, 
    const unsigned int min_selector,
    const unsigned int roo_selector,
    const double precison
)
{
    double dm =0.5;
    double l_x_lower = 0.0;
    double l_x_upper = 1.0;
    int status =0;
    int iter =0, max_iter = itetnum;
    UILC_fparams_linear l_params = {l.m,led_n};
    gsl_function F;
    F.function = &UILC_f_Morena_Linear;
    F.params = &l_params;

    if(fmod(l.m , 2.0)<DBL_EPSILON){ dm = 0.8;}
    else{dm = 0.6;}
/*
    if(l.m <3.0)
    {
        dm =0.8;
    }
    else if (l.m <4.0)
    {
        dm = 0.6;
    }
    else if(l.m < 6.0)
    {
        dm = 0.5;
    }
    else if(l.m < 12)
    {
        dm = 0.4;
    }
    else
    {
        dm = 0.2;
    }
*/


    if(GSL_IS_ODD(led_n))
    { // led_n = odd case we need to find the local minimum case
  
        /*
        There are 3 algorithms are provided in GSL library for find minimization of the function.
        - gsl_min_fminimizer_goldensection : The simplest method of bracketing the minimum of a function. It is the slowest algorithm provided by the library, with linear convergence.
        - gsl_min_fminimizer_brent : Using parabolic interpolation with the golden section algorithm. This produces a fast algorithm which is still robust.
        - gsl_min_fminimizer_quad_golden : This is a variant of Brent’s algorithm which uses the safeguarded step-length algorithm of Gill and Murray
        */
        const gsl_min_fminimizer_type * T = gsl_min_fminimizer_goldensection;
        switch(min_selector)
        {
            case 1: break;
            case 2: T = gsl_min_fminimizer_brent; break;
            case 3: T = gsl_min_fminimizer_quad_golden; break;
        }
        gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
        double f_min = UILC_f_Morena_Linear(dm,&l_params);
        double f_lower = UILC_f_Morena_Linear(l_x_lower,&l_params);
        double f_upper = UILC_f_Morena_Linear(l_x_upper,&l_params);
        gsl_min_fminimizer_set_with_values(s,&F,dm,f_min,l_x_lower,f_lower,l_x_upper,f_upper);

        /*
         printf ("using %s method\n",
          gsl_min_fminimizer_name (s));

         printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "min",
          "err", "err(est)");

         printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
          iter, l_x_lower, l_x_upper,
          dm, dm - 0.7, l_x_upper - l_x_lower);
        */

        do
        {
            iter++;
            status = gsl_min_fminimizer_iterate (s);
            dm = gsl_min_fminimizer_x_minimum (s);
            l_x_lower = gsl_root_fsolver_x_lower (s);
            l_x_upper = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval(l_x_lower, l_x_upper, 0.0001, precison );
        
            if (status == GSL_SUCCESS)
            {
                printf ("Converged:\n");
            }
                
            /*
            printf ("%5d [%.7f, %.7f] "
              "%.7f %+.7f %.7f\n",
              iter, l_x_lower, l_x_upper,
              dm, dm - 0.7, l_x_upper - l_x_lower);
            */
        }
        while(status == GSL_CONTINUE && iter < max_iter);
        
        gsl_min_fminimizer_free(s);
    }
    else
    { // n = even then we need to find the root of the function.
        /*
        There are 3 algorithms are provided in GSL library for find root without derivative.
        - gsl_root_fsolver_bisection : simplest method, slowest algorithms with linear convergence.
        - gsl_root_fsolver_falsepos : using linear interpolation, linear convergence, faster than bisection.
        - gsl_root_fsolver_brent : combines the interpolation strategy with the bisection algorithm. This produces a fast alorithm which is still robust.
        */
        const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;

        switch(roo_selector)
        {
            case 1: break;
            case 2: T = gsl_root_fsolver_falsepos; break;
            case 3: T = gsl_root_fsolver_brent; break;
        }
        gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
        gsl_root_fsolver_set(s,&F,l_x_lower,l_x_upper);

        do
        {
            iter++;
            status = gsl_root_fsolver_iterate (s);
            dm = gsl_root_fsolver_root (s);
            l_x_lower = gsl_root_fsolver_x_lower (s);
            l_x_upper = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval(l_x_lower, l_x_upper, 0, precison );
            if (status == GSL_SUCCESS)
            {
                printf ("Converged:\n");
            }
        }
        while(status == GSL_CONTINUE && iter < max_iter);

        gsl_root_fsolver_free(s);
        
    }

    return(dm);
}

double UILC_f_Morena_getdm_SquareGrid( // return the dm for Square Grid
    const UILC_Lamber_LED l, 
    const unsigned int led_n, 
    const unsigned int N, 
    const unsigned int M, 
    const unsigned int itetnum, 
    const unsigned int min_selector,
    const unsigned int roo_selector,
    const double precison
)
{
     // this function won't return the negative double vlaue unless there is an error
    if(led_n != N*M)
    {
        return -1.0;
    }

    double dm =0.0;
    double l_x_lower = 0.0;
    double l_x_upper = 1.0;
    int status =0;
    int iter =0, max_iter = itetnum;
    UILC_fparams_Rectangle R_param = {l.m, N, M};
    gsl_function F;
    F.function = &UILC_f_Morena_SquareGrid;
    F.params = &R_param;

    if( GSL_IS_ODD(N) && GSL_IS_ODD(M))
    { // both odd Minimization
        const gsl_min_fminimizer_type * T = gsl_min_fminimizer_goldensection;

        switch(min_selector)
        {
            case 0: break;
            case 1: T = gsl_min_fminimizer_brent; break;
            case 2: T = gsl_min_fminimizer_quad_golden; break;
        }

        gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
        gsl_min_fminimizer_set(s,&F,dm,l_x_lower,l_x_upper);

        do
        {
            iter++;
            status = gsl_min_fminimizer_iterate (s);
            dm = gsl_min_fminimizer_x_minimum (s);
            l_x_lower = gsl_root_fsolver_x_lower (s);
            l_x_upper = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval(l_x_lower, l_x_upper, 0, precison );
        }
        while(status == GSL_CONTINUE && iter < max_iter);
        gsl_min_fminimizer_free(s);
    } 
    else if(GSL_IS_EVEN(N) && GSL_IS_EVEN(M))
    { // both even: Root
        const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;

        switch(roo_selector)
        {
            case 0: break;
            case 1: T = gsl_root_fsolver_falsepos; break;
            case 2: T = gsl_root_fsolver_brent; break;
        }
        gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
        gsl_root_fsolver_set(s,&F,l_x_lower,l_x_upper);

        do
        {
            iter++;
            status = gsl_root_fsolver_iterate (s);
            dm = gsl_root_fsolver_root (s);
            l_x_lower = gsl_root_fsolver_x_lower (s);
            l_x_upper = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval(l_x_lower, l_x_upper, 0, precison );
        }
        while(status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fsolver_free(s);

    }
    else
    { // Minimum point or the root
        const gsl_min_fminimizer_type * T = gsl_min_fminimizer_goldensection;
        
        switch(min_selector)
        {
            case 0: break;
            case 1: T = gsl_min_fminimizer_brent; break;
            case 2: T = gsl_min_fminimizer_quad_golden; break;
        }

        gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
        gsl_min_fminimizer_set(s,&F,dm,l_x_lower,l_x_upper);

        do
        {
            iter++;
            status = gsl_min_fminimizer_iterate (s);
            dm = gsl_min_fminimizer_x_minimum (s);
            l_x_lower = gsl_root_fsolver_x_lower (s);
            l_x_upper = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval(l_x_lower, l_x_upper, 0, precison );
        }
        while(status == GSL_CONTINUE && iter < max_iter);
        gsl_min_fminimizer_free(s);


        /*
        There are three case, if local minimum or root exist in region.
        i) Minimum exists: Fine the dm has a minimum value.
        ii) Minimum = root : Fine also
        iii) Minimum <0: We nned to find a first root.
        */
        if(UILC_f_Morena_SquareGrid(dm,F.params) <0.0) 
        {
            const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;
            
            switch(roo_selector)
            {
            case 0: break;
            case 1: T = gsl_root_fsolver_falsepos; break;
            case 2: T = gsl_root_fsolver_brent; break;
            }
            
            gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
            gsl_root_fsolver_set(s,&F,l_x_lower,l_x_upper);
               
            l_x_upper = dm; // if dm <0.0 then there is a first root on interval and it is enought to serach 0 to dm;minium point.
            
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                dm = gsl_root_fsolver_root (s);
                l_x_lower = gsl_root_fsolver_x_lower (s);
                l_x_upper = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval(l_x_lower, l_x_upper, 0, precison );
            }
            while(status == GSL_CONTINUE && iter < max_iter);
            gsl_root_fsolver_free(s);

        }
    }

    return(dm);
}

UILC_LED_Arr UILC_f_Morena_get_Arr(
    const double dm, 
    double height,
    const unsigned int N, 
    const unsigned int M
)
{   
    double x=0.0;
    double y=0.0;
    double d = dm*height;
    //printf(" ARR DM: %le \n", d);
    gsl_vector * arr = gsl_vector_calloc( N * M *3);

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<M ; j++)
        {
            x=((double)i-((double)N-1.0)/2);
            y=((double)j-((double)M-1.0)/2);
            gsl_vector_set(arr,i*3 + 3*j+0,x*d) ;
            gsl_vector_set(arr,i*3 + 3*j+1,y*d) ;
            gsl_vector_set(arr,i*3 + 3*j+2, 0.0) ;
            //printf("%le x dm = %le\n",y,y*dm);
        }
    }
    UILC_LED_Arr Arr = {dm,height, d,0.0, N, M, arr};
    return(Arr);
}


inline double UILC_f_Morena_get_Morena_Boundary(
    UILC_LED_Arr arr,
    const int selector,
    const double md1,
    const double md2
)
{
    return(UILC_f_get_arr_target_Area(arr,selector)-md1*md2);
}

