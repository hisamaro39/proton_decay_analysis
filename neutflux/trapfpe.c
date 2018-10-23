#define _GNU_SOURCE 1
#include <fenv.h>
#include <fpu_control.h>
static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */

   fpu_control_t cw = (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE;
   _FPU_SETCW(cw);

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);


}

trapfpe_()
{
  trapfpe();
}
