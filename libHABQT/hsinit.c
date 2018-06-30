#include <stdlib.h>
#include "HsFFI.h"

static void HsStart(void) __attribute__((constructor));
static void HsStart(void)
{
  static char *argv[] = { "+RTS", "-A32m", NULL }, **argv_ = argv;
  static int argc = 2;
  hs_init(NULL,NULL);
}


static void HsEnd(void) __attribute__((destructor));
extern void HsEnd(void)
{
  hs_exit();
}
