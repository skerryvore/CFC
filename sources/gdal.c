#include "gdal.h"
#include "cpl_conv.h" /* for CPLMalloc() */
#include "cpl_string.h"
#include "ogr_api.h"


void create_tiff(char *tiffName, double adfGeoTransform[6], int vals[],\
                int nx, int ny, int utm, char *geoid)
{

char **papszMetadata;
char **papszOptions = NULL;
GDALDatasetH hDstDS;
const char *pszDstFilename = tiffName;
OGRSpatialReferenceH hSRS;
char *pszSRS_WKT = NULL;
GDALRasterBandH hBand;
GInt16 abyRaster[nx*ny];
int i;

for(i=0; i<nx*ny; i++)
{
abyRaster[i] = vals[i]; 
}


GDALAllRegister();

const char *pszFormat = "GTiff";
GDALDriverH hDriver = GDALGetDriverByName( pszFormat );

papszMetadata = GDALGetMetadata( hDriver, NULL );

hDstDS = GDALCreate (hDriver, pszDstFilename, nx, ny, 1, GDT_Int16, papszOptions );

GDALSetGeoTransform( hDstDS, adfGeoTransform );

hSRS = OSRNewSpatialReference( NULL );

OSRSetUTM( hSRS, utm, FALSE );

OSRSetWellKnownGeogCS( hSRS, geoid );

OSRExportToWkt( hSRS, &pszSRS_WKT );
OSRDestroySpatialReference( hSRS );

GDALSetProjection( hDstDS, pszSRS_WKT);

CPLFree( pszSRS_WKT );

hBand = GDALGetRasterBand( hDstDS, 1 );
GDALRasterIO(hBand, GF_Write, 0, 0, nx, ny, abyRaster, nx, ny, GDT_Int16, 0, 0);

GDALClose( hDstDS );
}


