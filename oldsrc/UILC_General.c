#include "../include/UILC.h"



/*-------------------------------------------------*/



/*-------------------------------------------------*/
extern inline double UILC_f_get_intensity_Lamber(
    UILC_Lamber_LED led, 
    const double theta
)
{
    return led.intensity * pow(cos(theta),led.m);
}

extern inline double UILC_f_get_intensity_Poly(
    UILC_Poly_LED led, 
    const double theta
)
{
    double y=0.0;
    for(int i =0; i < led.n+1;i++)
    {
        y+= *(led.param+i)*pow(theta,i);
    }
}

extern inline double UILC_f_get_intensity_Lamber_target(
    UILC_Lamber_LED ledmodel,
    gsl_vector * led,
    gsl_vector * target
)
{
    double vec[3] = {
        gsl_vector_get(led,0)-gsl_vector_get(target,0),
        gsl_vector_get(led,1)-gsl_vector_get(target,1),
        gsl_vector_get(led,2)-gsl_vector_get(target,2)};
    double r= gsl_hypot3(vec[0], vec[1], vec[2]);
    double pr = gsl_hypot(vec[0],vec[1]);
    double theta = atan2(pr,fabs(vec[2]));
    return (UILC_f_get_intensity_Lamber(ledmodel,theta) )/gsl_pow_2(r);
}

extern inline double UILC_f_get_intensity_Poly_target(
    UILC_Poly_LED ledmodel,
    gsl_vector * led,
    gsl_vector * target
)
{
    double vec[3] = {
        fabs(gsl_vector_get(led,0)-gsl_vector_get(target,0)),
        fabs(gsl_vector_get(led,1)-gsl_vector_get(target,1)),
        fabs(gsl_vector_get(led,2)-gsl_vector_get(target,2))};
    double r= gsl_hypot3(vec[0], vec[1], vec[2]);
    double pr = gsl_hypot(vec[0],vec[1]);
    double theta = atan2(pr,vec[2]);
    return (UILC_f_get_intensity_Poly(ledmodel,theta) )/gsl_pow_2(r);
}

//=>LED array function-------------------------------------------
UILC_LED_Arr UILC_f_Arr_calloc(
    const double dm,
    const double height,
    const unsigned int N, 
    const unsigned int M
)
{
    UILC_LED_Arr arr = {dm, height, height * dm, 0.0, N, M, gsl_vector_calloc(3*N*M)};
    return(arr);
}

int UILC_f_Arr_Lamber_set_max(
    UILC_LED_Arr * arr,
    UILC_Lamber_LED led
)
{
    gsl_vector * vec = gsl_vector_calloc(3);
    arr->Max_I = UILC_f_get_intensity_arr_Lamber_target(*arr, led, vec);
    return(0);
}

int UILC_f_Arr_free(
    UILC_LED_Arr arr
)
{
    gsl_vector_free(arr.coor);
    return(0);
}

gsl_vector * UILC_f_get_ArrCoordinate(
    UILC_LED_Arr arr,
    const unsigned int i,
    const unsigned int j
)
{
    int index =3*((i-1)*arr.M+j-1);
    gsl_vector_view vec = gsl_vector_subvector(arr.coor,index,3);

    return(&vec.vector);
}

int UILC_f_set_ArrCoordinate(
    UILC_LED_Arr arr,
    const unsigned int i,
    const unsigned int j,
    const gsl_vector * value
)
{

    int index =3*((i-1)*arr.M+j-1);
    gsl_vector_set(arr.coor,index, gsl_vector_get(value,0));
    gsl_vector_set(arr.coor,index+1, gsl_vector_get(value,1));
    gsl_vector_set(arr.coor,index+2, gsl_vector_get(value,2));
    
    return(0);
}

int UILC_f_set_AllArrCoordinate(
    UILC_LED_Arr arr,
    gsl_vector_view (*fill)(unsigned int, unsigned int)
) // Beaware of the 'fill' function form.
{
    int index=0;
    gsl_vector_view vec;
    for(int i=0; i< arr.N; i++ )
    {
        for(int j=0; j < arr.M; j++)
        {
            index =3*((i-1)*arr.M+j-1);
            vec = fill(i,j);
            if(vec.vector.data == NULL){ return(1);}

            gsl_vector_set(arr.coor,index, *(vec.vector.data));
            gsl_vector_set(arr.coor,index+1, *(vec.vector.data+1));
            gsl_vector_set(arr.coor,index+2, *(vec.vector.data+2));
        }
    }
    return(0);
}
extern inline double UILC_f_get_intensity_arr_Lamber_target(
    UILC_LED_Arr arr,
    UILC_Lamber_LED led,
    gsl_vector * target
)
{
    gsl_vector * vec = gsl_vector_calloc(3);
    double y =0.0;
    int index =0;
    for(int i=0;i< arr.N; i++)
    {
        for(int j=0; j< arr.M; j++)
        {
            index =3*((i)*arr.M+j);
            gsl_vector_set(vec,0,gsl_vector_get(arr.coor,index+0));
            gsl_vector_set(vec,1,gsl_vector_get(arr.coor,index+1));
            gsl_vector_set(vec,2,gsl_vector_get(arr.coor,index+2));

            y += UILC_f_get_intensity_Lamber_target(led,vec, target);
            /*
            printf("LED: (%le,%le,%le) Target: (%le,%le,%le) y= %le \n",
                    gsl_vector_get(vec,0),
                    gsl_vector_get(vec,1),
                    gsl_vector_get(vec,2),
                    gsl_vector_get(target,0),
                    gsl_vector_get(target,1), 
                    gsl_vector_get(target,2),
                    y
                    );
            */
        }
        
    }
    gsl_vector_free(vec);

    
    return(y);
}

extern inline double UILC_f_get_Normal_intensity_arr_Lamber_target(
    UILC_LED_Arr arr,
    UILC_Lamber_LED led,
    gsl_vector * target
)
{
    gsl_vector * vec = gsl_vector_calloc(3);
    gsl_vector_set(vec,2,arr.height);
    double max = UILC_f_get_intensity_arr_Lamber_target(arr,led, vec);
    double y =0.0;
    int index =0;
    for(int i=0;i< arr.N; i++)
    {
        for(int j=0; j< arr.M; j++)
        {
            index =3*((i)*arr.M+j);
            gsl_vector_set(vec,0,gsl_vector_get(arr.coor,index+0));
            gsl_vector_set(vec,1,gsl_vector_get(arr.coor,index+1));
            gsl_vector_set(vec,2,gsl_vector_get(arr.coor,index+2));

            y += UILC_f_get_intensity_Lamber_target(led,vec, target) / max ;
            /*
            printf("LED: (%le,%le,%le) Target: (%le,%le,%le) y= %le \n",
                    gsl_vector_get(vec,0),
                    gsl_vector_get(vec,1),
                    gsl_vector_get(vec,2),
                    gsl_vector_get(target,0),
                    gsl_vector_get(target,1), 
                    gsl_vector_get(target,2),
                    y
                    );
            */
        }
        
    }
    gsl_vector_free(vec);

    
    return(y);
}

extern inline double UILC_f_get_intensity_arr_Poly_target(
    UILC_LED_Arr arr,
    UILC_Poly_LED led,
    gsl_vector * target
)
{
    double y =0;
    gsl_vector_view vec;
    int index =0;
    for(int i=0;i< arr.N; i++)
    {
        for(int j=0; j< arr.M; j++)
        {
            index =3*((i)*arr.M+j);
            vec = gsl_vector_subvector(arr.coor,index,3);
            y += UILC_f_get_intensity_Poly_target(led,&vec.vector, target);
        }
    }
    gsl_vector_free(&vec.vector);
    return(y);
}

double  UILC_f_get_arr_target_Area(
    UILC_LED_Arr arr,
    const int selector
)
{ 
    double h = arr.d * (arr.M-1);
    //double h=gsl_vector_get(arr.coor,1); 
    //double w=fabs(gsl_vector_get(arr.coor,0)-gsl_vector_get(arr.coor,3*(arr.M-1)));
    double w= arr.d * (arr.N-1);

    if(arr.M == 1){
        h =1.0;
    }
    if(arr.N ==1){
        w = 1.0;
    }

    //printf("%d \n h0: %le h1: %le h: %le w: %le dm: %le\n",3*arr.M*(arr.N-1)+1,h0,h1, h , w ,dm);
    switch(selector)
    {
        case BC: // Set the edge LED location as the boundary edge
        return(h*w);
        break;
        case HDM: // add edge square which side is 1/2 dm
        return((h+arr.d/2)*(w+arr.d/2));
        break;
        case FDM: // add edge square which side is dm
        return((h+2*arr.d)*(w+2*arr.d));
    }

    
}

int UILC_f_print_arr(
    UILC_LED_Arr arr
)
{
    for(int i=0; i< arr.N; i++){
        for(int j=0; j< arr.M; j++){
            printf(" (%le,%le,%le)\n ", gsl_vector_get(arr.coor,3*((i)*arr.M+j)+0),gsl_vector_get(arr.coor,3*((i)*arr.M+j)+1),gsl_vector_get(arr.coor,3*((i)*arr.M+j)+2) );
        }
        printf("\n");
    }
    return(0);
}

//-----------------------------------------------------------------------

extern inline double UILC_df_get_intensity_arr_Lamber_target(
    const double x, 
    const double h, 
    UILC_LED_Arr arr,
    UILC_Lamber_LED led,
    const int axis
)
{
    int index =0;
    double result =0.0;
    if(axis ==0){    }
    else if(axis ==1){ index =1;}
    gsl_vector * vec = gsl_vector_calloc(3);
    gsl_vector_set(vec,2, arr.height);


    gsl_vector_set(vec,index,x-h);
    double fm1 = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
    gsl_vector_set(vec,index,x+h);
    double fp1 = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
    gsl_vector_set(vec,index,x-h/2);
    double fmh = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
    gsl_vector_set(vec,index,x+h/2);
    double fph =UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
    
    result = (4.0 / 3.0) * (fph - fmh) - (1.0 / 6.0) * (fp1 - fm1);
    gsl_vector_free(vec);
    return(result/h);
}

extern inline double UILC_df_get_Normal_intensity_arr_Lamber_target(
    const double x, 
    const double h, 
    UILC_LED_Arr arr,
    UILC_Lamber_LED led,
    const int axis
)
{
    int index =0;
    double result =0.0;
    if(axis ==0){    }
    else if(axis ==1){ index =1;}
    gsl_vector * vec = gsl_vector_calloc(3);
    gsl_vector_set(vec,2, arr.height);


    gsl_vector_set(vec,index,x-h);
    double fm1 = UILC_f_get_Normal_intensity_arr_Lamber_target(arr,led,vec);
    gsl_vector_set(vec,index,x+h);
    double fp1 = UILC_f_get_Normal_intensity_arr_Lamber_target(arr,led,vec);
    gsl_vector_set(vec,index,x-h/2);
    double fmh = UILC_f_get_Normal_intensity_arr_Lamber_target(arr,led,vec);
    gsl_vector_set(vec,index,x+h/2);
    double fph =UILC_f_get_Normal_intensity_arr_Lamber_target(arr,led,vec);
    
    result = (4.0 / 3.0) * (fph - fmh) - (1.0 / 6.0) * (fp1 - fm1);
    gsl_vector_free(vec);
    return(result/h);
}


//-----------------------------------------------------------------------

extern inline double UILC_f_find_derivative_Lamber(
    const int axis,
    UILC_LED_Arr arr,
    UILC_Lamber_LED led
)
{
    double x_i = 0.0;
    double x_i1 =0.0;
    double h = 0.0;
    double df_result =0.0; 
    double df_result_tem = 0.0;
    double df_abserr =0.0;
    int index = 0;

    const double height = arr.height;

    double return_value =0.0;

    if(axis == 0){    }
    else if (axis == 1){
        index = 1;
    }
    else{
        return(0.0);
    }
    x_i = gsl_vector_get(arr.coor,index);
    h = (gsl_vector_get(arr.coor,index+3)-gsl_vector_get(arr.coor,index))/1000;

    UILC_df_Lamber_param params = {height, h, axis, arr, led};

// Part 1. Find the x range (x_i, x_i1)such that f'(x) =1.0

    //printf("\n dm: %le, x_0 = %le ,\n (%le, %le, %le)\n", h *10.0, x_i, gsl_vector_get(arr.coor,0),gsl_vector_get(arr.coor,1), gsl_vector_get(arr.coor,2) );
    for(int i=0; x_i <0.0;i++)
    {
        df_result = UILC_df_get_intensity_arr_Lamber_target(x_i,h, arr,led,axis);
        if( i==0 && df_result < 1.0){
            printf("\nf'(x_init) < 1.0, df_result at %le = %le \n",x_i,df_result );
            return(0.0);
        }

        if(df_result >1.0){
            x_i += h;
        }
        else if(df_result <1.0 ){
            if(df_result <0.0){
                x_i -= h;
                h= h/10;
                x_i += h;
            }
            else{
                x_i1 = x_i;
                x_i -= h;
                break;
            }
            
        }
        else if(fabs(df_result-1.0)<DBL_EPSILON){
            printf("%le =1.0 at x= ", df_result,x_i);
            return(x_i);
        }
        printf("x: %le, df_dx: %le\n", x_i,df_result);

    }
// Part 2. Find the x such that f'(x) =1.0 using Bracket method.
    double size = x_i1 -x_i;
    double x_p = (x_i1 + x_i)/2;
    /*
    printf("x_i = %le, df=",x_i,  UILC_df_get_intensity_arr_Lamber_target(x_i,h, arr, led,axis));
    printf("x_p = %le, df=",x_p,  UILC_df_get_intensity_arr_Lamber_target(x_p,h, arr, led,axis));
    printf("x_i+1 = %le, df=",x_i1,  UILC_df_get_intensity_arr_Lamber_target(x_i1,h, arr, led,axis));

    printf("\n ------------------------------------------------------------------------------------------\n");
    printf("\n x_i: %le x_i1: %le x_p: %le df/dx(x_p): %le size: %le\n", x_i, x_i1,x_p,df_result,size);
    printf("\n ------------------------------------------------------------------------------------------\n");
    */
    do
    {
        //printf("==========================================================================================\n");
        //printf("x_i = %le, df=%le\n",x_i,  UILC_df_get_intensity_arr_Lamber_target(x_i,h, arr, led,axis));
        //printf("x_p = %le, df=%le\n",x_p,  UILC_df_get_intensity_arr_Lamber_target(x_p,h, arr, led,axis));
        //printf("x_i+1 = %le, df=%le\n",x_i1,  UILC_df_get_intensity_arr_Lamber_target(x_i1,h, arr, led,axis));
        //printf("==========================================================================================\n");
        
        df_result = UILC_df_get_intensity_arr_Lamber_target(x_p,h, arr, led,axis);
        if(fabs(df_result-1.0)<FLT_EPSILON){
            return(x_p);
        }
        else if(df_result >1.0){
            x_i = x_p;
            
        }
        else if(df_result <1.0){
            x_i1 = x_p;
            
        }
        

        size = x_i1 -x_i;
        x_p = (x_i1 +x_i)/2;

        params.h = params.h<size/5 ? params.h: size/5;
    }
    while(size > DBL_EPSILON);

    return_value = (x_i1 +x_i)/2;
    return(return_value);

}

extern inline double UILC_f_find_Normal_derivative_Lamber(
    const int axis,
    UILC_LED_Arr arr,
    UILC_Lamber_LED led
)
{
    double x_i = 0.0;
    double x_i1 =0.0;
    double h = 0.0;
    double df_result =0.0; 
    double df_result_tem = 0.0;
    double df_abserr =0.0;
    int index = 0;
    const double height = arr.height;

    double return_value =0.0;

    if(axis == 0){    }
    else if (axis == 1){
        index = 1;
    }
    else{
        return(0.0);
    }
    x_i = gsl_vector_get(arr.coor,index);
    h = (gsl_vector_get(arr.coor,index+3)-gsl_vector_get(arr.coor,index))/1000;

    UILC_df_Lamber_param params = {height, h, axis, arr, led};

// Part 1. Find the x range (x_i, x_i1)such that f'(x) =1.0

    //printf("\n dm: %le, x_0 = %le ,\n (%le, %le, %le)\n", h *10.0, x_i, gsl_vector_get(arr.coor,0),gsl_vector_get(arr.coor,1), gsl_vector_get(arr.coor,2) );
    for(int i=0; x_i <0.0;i++)
    {
        df_result = UILC_df_get_Normal_intensity_arr_Lamber_target(x_i,h, arr,led,axis);
        if( i==0 && df_result < 1.0){
            printf("\nf'(x_init) < 1.0, df_result at %le = %le \n",x_i,df_result );
            return(0.0);
        }

        if(df_result >1.0){
            x_i += h;
        }
        else if(df_result <1.0 ){
            if(df_result <0.0){
                x_i -= h;
                h= h/10;
                x_i += h;
            }
            else{
                x_i1 = x_i;
                x_i -= h;
                break;
            }
            
        }
        else if(fabs(df_result-1.0)<DBL_EPSILON){
            printf("%le =1.0 at x= ", df_result,x_i);
            return(x_i);
        }
        //printf("x: %le, df_dx: %le\n", x_i,df_result);

    }
// Part 2. Find the x such that f'(x) =1.0 using Bracket method.
    double size = x_i1 -x_i;
    double x_p = (x_i1 + x_i)/2;
    
    //printf("x_i = %le, df=",x_i,  UILC_df_get_intensity_arr_Lamber_target(x_i,h, arr, led,axis));
    //printf("x_p = %le, df=",x_p,  UILC_df_get_intensity_arr_Lamber_target(x_p,h, arr, led,axis));
    //printf("x_i+1 = %le, df=",x_i1,  UILC_df_get_intensity_arr_Lamber_target(x_i1,h, arr, led,axis));
    //printf("\n ------------------------------------------------------------------------------------------\n");
    //printf("\n x_i: %le x_i1: %le x_p: %le df/dx(x_p): %le size: %le\n", x_i, x_i1,x_p,df_result,size);
    //printf("\n ------------------------------------------------------------------------------------------\n");
    
    do
    {
        printf("==========================================================================================\n");
        printf("x_i = %le, df=%le\n",x_i,  UILC_df_get_Normal_intensity_arr_Lamber_target(x_i,h, arr, led,axis));
        printf("x_p = %le, df=%le\n",x_p,  UILC_df_get_Normal_intensity_arr_Lamber_target(x_p,h, arr, led,axis));
        printf("x_i+1 = %le, df=%le\n",x_i1,  UILC_df_get_Normal_intensity_arr_Lamber_target(x_i1,h, arr, led,axis));
        printf("==========================================================================================\n");
        
        df_result = UILC_df_get_Normal_intensity_arr_Lamber_target(x_p,h, arr, led,axis);
        if(fabs(df_result-1.0)<FLT_EPSILON){
            return(x_p);
        }
        else if(df_result >1.0){
            x_i = x_p;
            
        }
        else if(df_result <1.0){
            x_i1 = x_p;
            
        }
        

        size = x_i1 -x_i;
        x_p = (x_i1 +x_i)/2;

        params.h = params.h<size/5 ? params.h: size/5;
    }
    while(size > DBL_EPSILON);

    return_value = (x_i1 +x_i)/2;
    return(return_value);

}

extern inline double UILC_f_find_Normal_value_Lamber(
    const int axis,
    UILC_LED_Arr arr,
    UILC_Lamber_LED led,
    const double td
)
{
    printf("FIND %le point \n", td);
    
    double x_i = 0.0;
    double x_i1 =0.0;
    double h = 0.0;
    double df_result =0.0; 
    double df_result_tem = 0.0;
    double df_abserr =0.0;
    int index = 0;
    const double height = arr.height;

    double return_value =0.0;

    if(axis == 0){    }
    else if (axis == 1){
        index = 1;
    }
    else{
        return(0.0);
    }
    x_i = gsl_vector_get(arr.coor,index)-arr.d /2;
    h = (gsl_vector_get(arr.coor,index+3)-gsl_vector_get(arr.coor,index))/1000;

    gsl_vector * vec = gsl_vector_calloc(3);
    gsl_vector_set(vec,2,height);

    UILC_df_Lamber_param params = {height, h, axis, arr, led};

// Part 1. Find the x range (x_i, x_i1)such that f'(x) =1.0

    printf("\n dm: %le, x_0 = %le ,\n (%le, %le, %le)\n", h *10.0, x_i, gsl_vector_get(arr.coor,0),gsl_vector_get(arr.coor,1), gsl_vector_get(arr.coor,2) );
    for(int i=0; x_i <0.0;i++)
    {
        gsl_vector_set(vec,index,x_i);
        df_result = UILC_f_get_Normal_intensity_arr_Lamber_target(arr,led,vec);
        if( i==0 && df_result > td){
            printf("%le, %le\n", gsl_vector_get(vec,1),gsl_vector_get(vec,2));
            printf("\nf'(x_init) < 1.0, df_result at %le = %le \n",x_i,df_result );
            return(0.0);
        }
        //printf("x: %le f(x): %le \n", x_i, df_result);
        if(fabs(df_result-td)<FLT_EPSILON)
        {
            printf("%le =%le at x= %le \n", df_result,td,x_i);
            return(x_i);
        }
        else if(df_result >td){
            //printf("big\n");
            x_i -= h;
            h = h/10;
            x_i += h;
        }
        else if(df_result <td ){
               // printf("small\n");
                x_i += h;
        }
        //printf("x: %le, df_dx: %le\n", x_i,df_result);

    }
// Part 2. Find the x such that f'(x) =1.0 using Bracket method.
    double size = x_i1 -x_i;
    double x_p = (x_i1 + x_i)/2;
    
    //printf("x_i = %le, df=",x_i,  UILC_df_get_intensity_arr_Lamber_target(x_i,h, arr, led,axis));
    //printf("x_p = %le, df=",x_p,  UILC_df_get_intensity_arr_Lamber_target(x_p,h, arr, led,axis));
    //printf("x_i+1 = %le, df=",x_i1,  UILC_df_get_intensity_arr_Lamber_target(x_i1,h, arr, led,axis));
    //printf("\n ------------------------------------------------------------------------------------------\n");
    //printf("\n x_i: %le x_i1: %le x_p: %le df/dx(x_p): %le size: %le\n", x_i, x_i1,x_p,df_result,size);
    //printf("\n ------------------------------------------------------------------------------------------\n");
    
    do
    {
        gsl_vector_set(vec,index,x_p);
        df_result = UILC_f_get_Normal_intensity_arr_Lamber_target(arr,led,vec);
        if(fabs(df_result-td)<FLT_EPSILON){
            return(x_p);
        }
        else if(df_result >td){
            x_i1 = x_p;
            
        }
        else if(df_result <td){
            x_i = x_p;
        }
        

        size = x_i1 -x_i;
        x_p = (x_i1 +x_i)/2;

        params.h = params.h<size/5 ? params.h: size/5;
    }
    while(size > DBL_EPSILON);

    return_value = (x_i1 +x_i)/2;
    return(return_value);

}


/*-------------------------------------------------*/
/*
extern inline double UILC_f_find_derivative_Lamber(
    const double df_dx,
    const int axis,
    const double initialpoint, 
    const double endpoint,
    double step,
    UILC_LED_Arr arr,
    UILC_Lamber_LED led,
    const double height
)
{
    int N = (int)fabs(endpoint - initialpoint) / step;
    double h = step/2;
    step = (endpoint - initialpoint)< 0.0 ? -step: step;
    double fm1 = 0.0;
    double fp1 = 0.0;
    double fmh = 0.0;
    double fph = 0.0;
    double x =0.0;
    double d=0.0;
    gsl_vector * vec = gsl_vector_calloc(3);
    gsl_vector_set(vec,2, height);
    if(axis == 0)
    {// get x axis case
        for(int i=0; i<N;i++) 
        {   
            x = initialpoint+i*step;
            gsl_vector_set(vec,0, x-h);
            fm1 = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            gsl_vector_set(vec,0, x+h);
            fp1 = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            gsl_vector_set(vec,0, x-h/2);
            fmh = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            gsl_vector_set(vec,0, x+h/2);
            fph = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
        
            double result = (4.0/3.0) *(fph - fmh) - (1.0 / 3.0) *(0.5 * (fp1 - fm1));

            if(result >df_dx){
                d = x;
                break;
            }
        }

    }
    else if(axis == 1)
    {// get y axis case
        for(int i=0; i<N;i++) 
        {   
            
            x = initialpoint+i*step;
            printf("%d, %le \n",i+1,x);
            gsl_vector_set(vec,1, x-h);
            fm1 = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            gsl_vector_set(vec,1, x+h);
            fp1 = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            gsl_vector_set(vec,1, x-h/2);
            fmh = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            gsl_vector_set(vec,1, x+h/2);
            fph = UILC_f_get_intensity_arr_Lamber_target(arr,led,vec);
            printf("fm1: %le\nfp1: %le\nfmh: %le\nfph: %le\n", fm1,fp1,fmh,fph);
            double result = (4.0/3.0) *(fph - fmh) - (1.0 / 3.0) *(0.5 * (fp1 - fm1));
            printf("result: %le\n",result);
            if(result >df_dx){
                d = x;
                break;
            }
        }
    }

    return(d);
}
*/
inline double UILC_f_find_derivative_Poly(
    const double df_dx,
    const int axis,
    const double initialpoint, 
    const double endpoint,
    const double step,
    UILC_LED_Arr arr,
    UILC_Poly_LED led
)
{
    int N = (int)(endpoint - initialpoint) / step;
    double h = step/2;
    double fm1 = 0.0;
    double fp1 = 0.0;
    double fmh = 0.0;
    double fph = 0.0;
    double x =0.0;
    double d=0.0;
    gsl_vector * vec = gsl_vector_calloc(3);

    if(axis == 0)
    {// get x axis case
        for(int i=0; i<N;i++) 
        {   
            x = initialpoint+i*step;
            gsl_vector_set(vec,0, x-h);
            fm1 = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
            gsl_vector_set(vec,0, x+h);
            fp1 = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
            gsl_vector_set(vec,0, x-h/2);
            fmh = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
            gsl_vector_set(vec,0, x+h/2);
            fph = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
        
            double result = (4.0/3.0) *(fph - fmh) - (1.0 / 3.0) *(0.5 * (fp1 - fm1));

            if(fabs(result - df_dx) < DBL_EPSILON){
                d = result;
                break;
            }
        }

    }
    else if(axis == 1)
    {// get y axis case
        for(int i=0; i<N;i++) 
        {   
            x = initialpoint+i*step;
            gsl_vector_set(vec,1, x-h);
            fm1 = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
            gsl_vector_set(vec,1, x+h);
            fp1 = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
            gsl_vector_set(vec,1, x-h/2);
            fmh = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
            gsl_vector_set(vec,1, x+h/2);
            fph = UILC_f_get_intensity_arr_Poly_target(arr,led,vec);
        
            double result = (4.0/3.0) *(fph - fmh) - (1.0 / 3.0) *(0.5 * (fp1 - fm1));

            if(fabs(result - df_dx) < DBL_EPSILON){
                d = result;
                break;
            }
        }
    }

    return(d);
}


/*----------------------------------------------------------------------------------------*/

