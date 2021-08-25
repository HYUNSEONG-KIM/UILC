#include "uic.h"

#include <math.h>
#include <stdblib.h>
#include <float.h>


double __intensity_2d(double i, double s, double x, double t, double phi, double h)
{
    return i*pow(pow(h,2)+pow((t-x),2), -(s+2)/2) * pow(h*cos(phi)+(t-x)*sin(phi),s);
}

double __intensity_center(double i, double s, doubel x, double phi, double h)
{
    return 2*__intensity_2d(i,s,x,0,phi,h);
}

double __intensity_boundary(double i, double s, doubel x, double phi, double h, double w)
{
    return __intensity_2d(i,s,x,w/2,phi,h)+__intensity_2d(i,s,-x,w/2,phi,h);
}

double __D(double a, double s, double d )
{
    return pow((1+((1/2)*a + d)^2),-((s+2)/2)) + pow((1+((1/2)*a - d)^2),-((s+2)/2))-2*pow((1+(d)^2),-((s+2)/2));
}

double __Di(double I0,  double s, double x, double h, double w)
{
    return (I0/pow(h,2)) *__D(w/h, s, x/h);
}

double _uic_esc_linear_cof(const double x, void * p)
{
    _uic_esc_fparams_linear_cof * params = (_uic_esc_fparams_linear_cof *)p;
    const unsigned int N = p->N;
    const double s = p->s;
    double result = 0.0;
    for(int i=1; i<N+1, i++)
    {
        result += (1-(s+3)*pow2( (N+1-2*(i))*(x/2) ) ) / pow(( 1+ pow2((x/2)*(N+1-2*(i)))),(s/2+3));
    }

}

double _uic_esc_rectangular_cof(const double x, void * p)
{
    _uic_esc_fparams_rectangular_cof * params = (_uic_esc_fparams_rectangular_cof  *)p;
    const unsinged int N = p->N;
    const unsinged int M = p->M;
    double result =0.0;

    for(int i= 1;i < N+1 ; i++)
    {
        for(int j=1;j< M+1 ;j++)
        {
            result += ;
        }
    }

    
}
/*=======================================================================================*/

/*esc: Expanded Sparrow's Criterion------------------------------*/


double uic_esc_coefficient_linear(const double s, const int n)
{

}

double uic_esc_coefficient_Rectangular(const double s, const int n, const int m)
{
    
}

LEDmatrix * uic_esc_array(const double s, const int n, const int m)
{
    
}


/*bcmt: Boundary Center Matching method------------------------------*/


/*Normalized functions------------------------------*/
double uic_bcm_de(double a, double s)
{

}

double uic_bcm_dm(double de, double a, double s, double h)
{
    
}

double uic_bcm_dm_approx(double s)
{
    return sqrt(pow(2, 2/(s+2))-1);
}

double uic_bcm_de_approx(double a, double s)
{
    return 1/6 * uic_bcm_dm_approx(s) + 0.25*a;
}



/*General functions------------------------------*/
double uic_bcm_xe(double s, double h, double w)
{
    return h*uic_bcm_de(w/h,s);
}

double uic_bcm_xm(double s, double h, double w, double xe)
{
    double a = w/h;
    return h*uic_bcm_de(a,s);
}

double uic_bcm_xe_approx()
{
    double a = w/h;
    return h * uic_bcm_de_approx(a, s);
}

double uic_bcm_xm_approx()
{
    double a = w/h;
    return h*uic_bcm_dm_approx(s);
}

LEDmatrix * uic_bcmt_array(const double s, ,const double I, const double w1, const double w2, const double h)
{


}

