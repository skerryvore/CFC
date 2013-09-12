
#include <stdlib.h>
#include <string.h>
#include "shapefil.h"


void CSHPpoints( int nSHPType, const char *pszFilename, double sampx[], double sampy[], int nsamp)

{
    SHPHandle	hSHPHandle;
    SHPObject	*psShape;
    double	x, y, z, m;
    int i;

    hSHPHandle = SHPCreate( pszFilename, nSHPType );

    for (i = 0; i < nsamp; i++)
    {
      x = sampx[i];
      y = sampy[i];
      z = 0.0;
      m = 0.0;
      psShape = SHPCreateObject( nSHPType, -1, 0, NULL, NULL,
                               1, &x, &y, &z, &m );
      SHPWriteObject( hSHPHandle, -1, psShape );
      SHPDestroyObject( psShape );
    }

    SHPClose( hSHPHandle );
}
