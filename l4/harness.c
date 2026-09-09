/* Experimental AENET API smoke check; MPL-2.0, see src/license-header.txt. */
#include "aenet.h"
#include <stdio.h>
int main(void) {
  char *types[]={"Cu","Au"}; char *incoming[]={"Au","Cu"};
  int status=-1, input[]={1,2}, output[]={0,0};
  aenet_init(2,types,&status); if(status!=AENET_OK)return 1;
  if(aenet_all_loaded())return 2;
  aenet_convert_atom_types(2,incoming,2,input,output,&status);
  if(status!=AENET_OK || output[0]!=2 || output[1]!=1)return 3;
  aenet_final(&status);if(status!=AENET_OK)return 4;
  puts("C API initialization/type mapping/finalization passed");return 0;
}
