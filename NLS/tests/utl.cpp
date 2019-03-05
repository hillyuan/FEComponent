#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utl.h"

/* =========================== read_string =============================== */

int read_string ( FILE *in, char string[80] )
{
   int i, c;

   while ( ( c = fgetc ( in ) ) != '\'' ) {   /* finds first '        */
      if ( c == EOF )
      return ( 0 );
   }
   for(i=0;(c=fgetc(in)) != '\'';i++)  {      /* fill string until next '  */
      if ( c == EOF )
         return ( 0 );
      string[i] = c;
   }
   string[i] = '\0';
   return ( 1 );

}  /* End of Read_string  */
