#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>

int main ()
{
  gsl_ieee_env_setup ();
  float a = 1e-33;
  while (a > 0){
    a = a/(2.0);
    gsl_ieee_printf_float(&a);
    printf("\n");
  }
  return 0;
}