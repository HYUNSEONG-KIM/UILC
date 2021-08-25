#include "../include/UILC.h"
#include <string.h>


#define CLEARINBUFFER while(getchar() != '\n');


// Generate the Sample data of give m and H

/*
given -n- and Morena solution

1. Determine the method wheter the simple inverse matrix or Active set method
2. Show the graph of the function including 
   - the vertical lines that indicate the location of LED
   - The boundary of the LED set: last lled + 1/2 * dm
   - The boundary of the Given area
   - 
*/
int main(void)
{
   // LED setting
   UILC_Lamber_LED led = {1.0E10, 4}
   int led_n = 1;
   // Area setting
   double height =30 *1.0E-3;
   double width = 200 * 1.0E-3;
   //Parameters for Morena dm
   int min_selector = MIN_GOLDENSECTION;
   int roo_selector = ROO_BISECTION;
   int iter = 1000;

   int selector =0;

   double dm = UILC_f_Morena_getdm_Linear(led,led_n,iter,min_selector,roo_selector,DBL_EPSILON)
   UILC_LED_Arr arr = UILC_f_Morena_get_Arr(dm,height,1,init_led);


   printf("Setting: Lamber LED \n");
   printf("LED Intensity: %10E, LED cos number: %f \n", led.intensity, led.m);
   printf("Width: %10e, Height: %10e\n", height, width);
   printf("led number : %d\n dm : %e, d : %e", led_n, dm, dm * height);

   //Setting the parameter or Calculate

   printf("1. Reset the parameters\n");
   printf("2. Calculate the graph\n");
   CLEARINBUFFER
   scanf("%d",&selctor);
   if(selector == 1)
   {
      printf("Set Lamber LED:");
      scanf("%e %e", &led.intensity, &led.m);
      printf("Width Height:");
      scanf("%e %e",&height, &width);
      printf("LED number:");
      scanf("%d", &led_n);
   }
   else if
   {
      
   }
   
   return(0);

}