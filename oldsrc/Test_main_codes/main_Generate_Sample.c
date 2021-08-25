#include "../include/UILC.h"
#include <string.h>


#define CLEARINBUFFER while(getchar() != '\n');

int main(void)
{
    
    //

/*
    printf("Enter the Intensity and the Lamber number:");
    scanf("%le %le",&(led.intensity), &(led.m));
    // printf("Enter the line length and target distance:");
    // scanf("%le %le", &Width, &D);
    printf("Enter the initial Morena led number(as many as possible):");
    scanf("%d",&init_led);
    printf("Select the min/root method sep with space, Default is 1,1\n   Root      Minimum\n1: BISECTION GOLDENSECTION\n2: FALSEPOS  BRENT_Min\n3: BRENT_Roo QUADGOLDEN\nExample, Root: FASLEPOS, Min: FOLDENSECTION: \"2 1\"\nenter: ");
    scanf("%d %d", &roo_selector, &min_selector);
    roo_selector = ( (roo_selector==ROO_BISECTION) || (roo_selector==ROO_FALSEPOS) || (roo_selector==ROO_BRENT) ) ? roo_selector: ROO_BISECTION;
    min_selector = ( (min_selector==MIN_GOLDENSECTION)||(min_selector==MIN_BRENT)||(min_selector==MIN_QUADGOLDEN)  ) ? min_selector : MIN_GOLDENSECTION;
    printf("Algorithm iteration (default 1000):");
    scanf("%d",&iter);
    iter = iter <0 ? 1000 : iter;
*/
    printf("Start\n\n");
    UILC_Lamber_LED led = {0.0,1};
    double height = 0.5;
    unsigned int init_led = 1;
    int min_selector=MIN_GOLDENSECTION;
    int roo_selector=ROO_BISECTION;
    int iter = 1000;

    led.intensity =1000000;

    roo_selector =1;
    min_selector =1;
    iter =1000;
    gsl_vector * vec = gsl_vector_calloc(3);

    UILC_LED_Arr arr;

    double start = 0.0;
    double end = 0.0;
    double step = 0.0;

    FILE *fp;
    FILE *fp2;
    char filename[100]="./data/data.csv";
    char filename2[100]="";
    
    fp = fopen(filename, "w+");
    
    fprintf(fp,"m, LED_intensity, Middle_Intensity, #_led, dm, height, D, Area_L, Area_HD, Area_FD, 90%\n");
    //fprintf(fp,"x,f(x)\n");
    /*
    height = 0.05;
    init_led  =5;
    double dm = UILC_f_Morena_getdm_Linear(led,init_led,iter,min_selector,roo_selector,DBL_EPSILON);
    arr = UILC_f_Morena_get_Arr(dm,height,1,init_led);
    gsl_vector_set(vec,2,height);
    arr.Max_I = UILC_f_get_intensity_arr_Lamber_target(arr, led, vec);

    printf("%le %le %le ->  %le\n",
            gsl_vector_get(vec,0),
            gsl_vector_get(vec,1),
            gsl_vector_get(vec,2),
            arr.Max_I
            );


    double start = gsl_vector_get(arr.coor,1);
    double end = -1.0*start;
    double step = arr.d/10000;

    printf("start: %le\n End: %le \n", start, end);
    while(start < end){
        gsl_vector_set(vec,1,start);
        fprintf(fp,"%le, %le\n",start, UILC_f_get_intensity_arr_Lamber_target(arr, led, vec));
        start += step;
    }
    */
    led.m =2.0;
    init_led = 2;
    while(led.m<101.0)
    {
        height = 0.05;
        while(height < 0.5)
        {
            init_led =2;
            while(init_led<30)
            {
                printf("      m,  #_led, height\n");
                printf("case: %le, %d, %le\n", led.m, init_led, height);

                double dm = UILC_f_Morena_getdm_Linear(led,init_led,iter,min_selector,roo_selector,DBL_EPSILON);
                arr = UILC_f_Morena_get_Arr(dm,height,1,init_led);
                gsl_vector_set(vec,2,height);
                arr.Max_I = UILC_f_get_intensity_arr_Lamber_target(arr, led, vec);

                sprintf(filename2,"./data/plot/p_%d_%d_%.5f.csv",(int)led.m,init_led,height);

                fp2 = fopen(filename2, "w+");
                fprintf(fp2, "x, f(x)\n");
                start = gsl_vector_get(arr.coor,1)-arr.d;
                end = -1.0*start;
                step = arr.d/10000;

                printf("start: %le\n End: %le \n", start, end);
                while(start < end){
                    gsl_vector_set(vec,1,start);
                    fprintf(fp2,"%le, %le\n",start, UILC_f_get_intensity_arr_Lamber_target(arr, led, vec));
                    start += step;
                }
                fclose(fp2);
                
                //int i = 1;
                //int j = 1;
                //int index =3*((i-1)*arr.M+j-1);
                //gsl_vector_view v = gsl_vector_subvector(arr.coor,index,3); //이건 새로운 벡터를 만드는게 아니라 그냥 원래 백터에서 index 부터 + strid 까지를 가르킨 주소만 가져오는듯
                //gsl_vector * vec = gsl_vector_calloc(3);
                //gsl_vector_memcpy(ve, &v.vector);
                //gsl_vector_set(ve,2,0.3);
                //gsl_vector_set(ve,1,0.0);
                //UILC_f_print_arr(arr);
                
                fprintf(fp,"%le, %le, %le, %d, %le, %le, %le, %le, %le, %le, %le\n",
                            led.m,
                            led.intensity,
                            arr.Max_I,
                            init_led,
                            arr.dm,
                            arr.height,
                            arr.d,
                            UILC_f_get_arr_target_Area(arr,BC),
                            UILC_f_get_arr_target_Area(arr,HDM),
                            UILC_f_get_arr_target_Area(arr,FDM),
                            2*fabs(UILC_f_find_Normal_value_Lamber(1, arr, led, 0.9))
                            );
            
                

                UILC_f_Arr_free(arr);
                printf("Done\n");

                init_led++;
            }
            height += 0.05;
        }
        led.m++;
    }   

    printf("End\n");
    gsl_vector_free(vec);

    fclose(fp);

    /*

    double start = gsl_vector_get(arr.coor,1)-dm;
                double end = fabs(start);
                double dmh = dm/10000;
                gsl_vector_set(ve,2,height);
    printf("start:%le, dmh = %le \n",start, dmh);

    gsl_vector* vec = gsl_vector_calloc(3); // location of each led
    double y =0.0;
    index =0;
    for(int i=0;i< arr.N; i++)
    {
        for(int j=0; j< arr.M; j++)
        {
            index =3*((i)*arr.M+j);
            gsl_vector_set(vec,0,gsl_vector_get(arr.coor,index+0));
            gsl_vector_set(vec,1,gsl_vector_get(arr.coor,index+1));
            gsl_vector_set(vec,2,gsl_vector_get(arr.coor,index+2));

            double vv[3] = {
                gsl_vector_get(vec,0),
                gsl_vector_get(vec,1),
                gsl_vector_get(vec,2)};
            double vvv[3] = {
                gsl_vector_get(ve,0),
                gsl_vector_get(ve,1),
                gsl_vector_get(ve,2)
            };   

            printf("(x,y,z) \n1: (%le, %le, %le)\n2: (%le, %le, %le)\n",
                vv[0],vv[1],vv[2],
                vvv[0], vvv[1],vvv[2]);

            double vect[3] = {
                (gsl_vector_get(vec,0)-gsl_vector_get(ve,0)),
                (gsl_vector_get(vec,1)-gsl_vector_get(ve,1)),
                (gsl_vector_get(vec,2)-gsl_vector_get(ve,2))};

            printf("Vect[i]: %le %le %le\n", vect[0], vect[1], vect[2]);
            double r= gsl_hypot3(vect[0], vect[1], vect[2]);
            double pr = gsl_hypot(vect[0],vect[1]);
            double theta = atan2(pr,fabs(vect[2]));
            
            double p= pow(cos(theta),led.m);
            double r2 = gsl_pow_2(r);
            y += led.intensity * p/r2;
            printf("inf test \n p = %le\n r^2 = %le\n", p, r2);
            printf("LED: (%le,%le,%le) Target: (%le,%le,%le) y= %le \n",
                    gsl_vector_get(vec,0),
                    gsl_vector_get(vec,1),
                    gsl_vector_get(vec,2),
                    gsl_vector_get(ve,0),
                    gsl_vector_get(ve,1),
                    gsl_vector_get(ve,2),
                    y
                    );
        }
        
    }
    gsl_vector_free(vec);

    //printf("%le \n",y);
    // 미분계수 1인 함수 찾는 과정이 제대로 안됨 확인 필요.
    //double df_dx_1 = UILC_f_find_derivative_Lamber(1.0,0,gsl_vector_get(vec,0)- dm,0.0,0.00001,arr,led);
    //double df_dy_1 = UILC_f_find_derivative_Lamber(1.0,1,gsl_vector_get(vec,1)+ dm,0.0,0.00001,arr,led);
    
    //printf("df/dx=1: %le, df/dy =1: %le \n",df_dx_1,df_dy_1);
    
    //double Morena_Bound = UILC_f_Morena_get_Morena_Boundary(arr,FDM,df_dx_1,df_dy_1);

*/
 
 /*  
    int i=0;
    while(i>=0)
    {
        printf("Discrete step:");
        scanf("%d",&i);
        
        
        
    }
    */
   return(0);

}