#include "../include/UILC.h"

// Verification is needed

inline static int s_f_sort_compare(void *a, void * b)
{
  return (((int *)a > (int *)b)? 1:0); 
}

inline static int s_ASM_R_initiaize(int * R, n)
{
    for(init i=0; i < n ; i++)
    {
        R[i] = i+1;
    }
    return(0);
}

inline static s_ASM_w_setting(
    gsl_vector * w, 
    const gsl_matrix * A, 
    const gsl_vector * b, 
    const gsl_vector *x, 
    int n
)
{
    gsl_vector * y = gsl_vector_calloc(n);
    gsl_vector_memcpy(y ,b);

    gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, y);
    gsl_blas_dgemv(CblasTrans, 1.0, A, y, 0.0, w);

    return(0);
}

inline static int s_ASM_get_max_w_R(
    gsl_vector * w, 
    int * R, 
    &max_w_R, 
    int n//length of the R
)
{
    double num = gsl_vector_get(w,0);
    double inum =0.0;
    int j =0;
    for(init i=0; i < n ; i++)
    {
        if(i ==j )
        {
            inum =gsl_vector_get(w,R[i]);
            if(inum > num )
            {
                num = inum;
            }
            j = i;
        }

    }
    *max_w_R = num;
    return(j);
}

//-------------------------------------------------
inline static int s_ASM_R2P(int * R, int * P, int j,int nr, int np)
{
  P[np] = j;
  np ++;
  qsort(P, np, sizeof(int), s_f_sort_compare);
  
  int k=0;
  for(int i =0; i< nr, i++)
  {
    if((R[i] == j || k ==1))
    {
      k=1;
      
      if(i != nr-1)
      {
        R[i] = R[i+1];
      }
      else
      {
        R[i] = 0;
      }
      
    }
  }
  
  return(0);
  
}
//-------------------------------------------------

//-------------------------------------------------
inline static gsl_matrix * s_ASM_get_A_P(gsl_matrix * A, int * P, int np, int n)
{
    gsl_matrix * A_P = gsl_matrix_calloc(n,np);
    for(int i=0; i < np; i ++)
    {
        gsl_matrix_set_col(A_P, i , gsl_matrix_get_col(A,P[i]))
    }
    
    return(A_P);
}
//-------------------------------------------------

//-------------------------------------------------
inline static s_ASM_s_P_setting(s,P,A_P,b)
{
    gsl_vector * y = gsl_vecotr_calloc(n);
    gsl_matrix * As  = gsl_matrix_calloc(n,n);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0,A_P,A_P,0.0,As);
    gsl_blas_dgemv(CblasTrans,1.0,A_P,y,0.0,y);
    gsl_blas_dgemv(CblasNoTrans,1.0,As,y,0.0,y);
    
    
    for(int i =0; i < np, i++)
    {
        gsl_vector_set(s, P[i], gsl_vector_get(y,i));
    }
    return(0);
}
//-------------------------------------------------

//-------------------------------------------------
inline static s_ASM_s_R_setting(
  gsl_vector *s,
  int *R
  int nr)
{
  for(int i=0; i<nr; i++)
  {
    gsl_vector_set(s,R[i],0.0); 
  }
}

//-------------------------------------------------

//-------------------------------------------------
inline static double s_ASM_a_min(
  gsl_vector *x,
  gsl_vector *s,
  int *P
  int np
)
{
  double min = x[0]/(x[0]-s[0]);
  double g = 0.0;
  for(int i=0; i<np;i++)
  {
    if(s[P[i]] < 0 || fabs(s[P[i]] - 0.0) < DBL_EPSILON)
    {
      g = x[i]/(x[i]-s[i]);
      
      if(min > g)
      {
        min = g;
      }
    }
  }
  return(min);
}
//-------------------------------------------------



//-------------------------------------------------
inline static s_ASM_R2P_all(int *R,int *P, int nr, int np )
{
  int j=0;
    for(int i=0; i<np ; i++)
    {
      if(x[P[i]] <0 || fabs(x[P[i]] - 0.0 ) <DBL_EPSILON)
      {
        R[nr+j]=P[i];
        j++;
      }
    }
    qsort(R, nr, sizeof(int), s_f_sort_compare);
    return(0);
}
//-------------------------------------------------

inline static int s_ASM_get_size_RorP(int *R, n)
{
    int i =0; 
    int num =0;
    for(i; i<n ; i++)
    {
        if(R[0] != 0)
        {
            num++;
        }
    }
    return(num);
}


inline int UILC_f_Hyeon_s_ASM(
    gsl_matrix * A, 
    gsl_vector * b, 
    gsl_vector * result, 
    double epsilon)//Active Set method
{


    //Initiaize:
    int m = A->size1;
    int n = A->size2;
    if(b->size != m){
        printf("Not a right system, dimension is differernt\n");
        return(0);
    }

    int * P = malloc(sizeof(int) * n);
    int * R = malloc(sizeof(int) * n);

    s_ASM_R_initiaize(R,n);

    gsl_vector * x = gsl_vector_calloc(n);
    gsl_vector * s = gsl_vector_calloc(n);
    gsl_vector * w = gsl_vector_calloc(n);
    s_ASM_w_setting(w, A, b, x, n);

    double max_w_R=0.0;
    int j =0.0;

    do
    {
        j = s_ASM_get_max_w_R(w, R, &max_w_R );
        s_ASM_R2P(R,P,j);

        gsl_matrix * A_P = s_ASM_get_A_P(A, P);

        s_ASM_s_P_setting(s,A_P,b);
        s_ASM_s_R_setting(s,A_P,b);
        
        do
        {
            double a = s_ASM_min_a(x,s,P);
            gsl_vector_axpby(a,s, 1-a,x);
            s_ASM_R2P_all(R,P);
            s_ASM_s_P_setting(s,A_P,b);
            s_ASM_s_R_setting(s,A_P,b);

        } while (/* condition */);
    
        gsl_vector_memcpy(x,s);
        s_ASM_w_setting(A, b, x, n);
        
    } while ( s_ASM_size(R) != 0 && s_ASM_max_w_R(w,R) > epsilon);

    gsl_vector_memcpy(result,x);
}