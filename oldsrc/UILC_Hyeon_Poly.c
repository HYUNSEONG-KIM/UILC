inline double UILC_f_s_matrix_value_Poly(UILC_Poly_LED led, int i, int j, double d, double H)
{
    double theta = arctan(abs(i-j)*d /H);
    return(UILC_f_get_intensity_Poly(led, theta));    
}

extern inline int UILC_f_s_matrix_setting_Poly(gsl_matrix * A,  double (*f)(UILC_Lamber_LED led, int i, int j))
{
    int m = A->size1;
    int n = A->size2;

    for(int i=0; i<m, ;i++)
    {
        for(int j=0; j <n; j++)
        {
            gsl_matrix_set(A,i, j, f(led,i,j));
        }
    }
    return(GSL_SUCESS);
}