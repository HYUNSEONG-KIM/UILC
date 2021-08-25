/*  uic.h 
 * 
 * Copyright (C) 2007 Kim, Hyeansung
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

//


#ifndef __UIC_H__
#define __UIC_H__


#include <float.h>
#include <math.h>

#include<gsl/gsl_min.h>
#include<gsl/errno.h>
#include<gsl/gsl_roots.h>
#include<gsl/gsl_min.h>



#ifdef __cplusplus
#define __BEGIN_DECLS extern "C"{
#define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


#ifdef HAVE_INLINE
#endif /* HAVE_INLINE */

#define pow2(x) (x)*(x)

typedef struct 
{
    /* data */
    double intensity;
    double s;
    double phi;
    double theta;
    gsl_vector * location;
    gsl_vector * orietation;
}
LED;


typedef struct 
{
    int n;
    int N;
    int M;
    double s;
    LED * ledlist;
}
LEDmatrix;

typedef struct 
{
    /* data */
}_uic_esc_fparams_linear_cof;

typedef struct 
{
    /* data */
}_uic_esc_fparams_rectangular_cof;





