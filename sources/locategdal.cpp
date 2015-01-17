#include "gdal.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "cpl_minixml.h"

extern "C"
{
void getElevation(double *dfGeoX, double *dfGeoY, char *Filename, double *value);
char *SanitizeSRS( const char *pszUserInput );
}



void getElevation(double *dfGeoX, double *dfGeoY, char *Filename, double *value)

{

  char          *pszSourceSRS = NULL;
  const char    *pszSrcFilename = NULL;


  GDALAllRegister();
  /* -------------------------------------------------------------------- */
  /*      Open source file.                                               */
  /* -------------------------------------------------------------------- */
  pszSrcFilename = Filename;
  GDALDatasetH hSrcDS = NULL;
  hSrcDS = GDALOpen( pszSrcFilename, GA_ReadOnly );

  CPLFree(pszSourceSRS);
  pszSourceSRS = SanitizeSRS("WGS84");
  OGRSpatialReferenceH hSrcSRS = NULL, hTrgSRS = NULL;
  OGRCoordinateTransformationH hCT = NULL;
  hSrcSRS = OSRNewSpatialReference( pszSourceSRS );
  hTrgSRS = OSRNewSpatialReference( GDALGetProjectionRef( hSrcDS ) );
  hCT = OCTNewCoordinateTransformation( hSrcSRS, hTrgSRS );

  if (hCT)
  {
  OCTTransform( hCT, 1, dfGeoX, dfGeoY, NULL );
  }

  double adfGeoTransform[6], adfInvGeoTransform[6];

  GDALGetGeoTransform( hSrcDS, adfGeoTransform );
  GDALInvGeoTransform( adfGeoTransform, adfInvGeoTransform );

  int iPixel, iLine;

  iPixel = (int) floor(
      adfInvGeoTransform[0] 
      + adfInvGeoTransform[1] * *dfGeoX
      + adfInvGeoTransform[2] * *dfGeoY );
  iLine = (int) floor(
      adfInvGeoTransform[3] 
      + adfInvGeoTransform[4] * *dfGeoX
      + adfInvGeoTransform[5] * *dfGeoY );

  GDALRasterBandH hBand = GDALGetRasterBand( hSrcDS, 1 );

  int iPixelToQuery = iPixel;
  int iLineToQuery = iLine;


  double adfPixel[2];
  GDALRasterIO( hBand, GF_Read, iPixelToQuery, iLineToQuery, 1, 1,
      adfPixel, 1, 1, GDT_CFloat64, 0, 0);

  *value = adfPixel[0];

    if (hCT) {
        OSRDestroySpatialReference( hSrcSRS );
        OSRDestroySpatialReference( hTrgSRS );
        OCTDestroyCoordinateTransformation( hCT );
    }
  GDALClose(hSrcDS);

}


char *SanitizeSRS( const char *pszUserInput )

{
    OGRSpatialReferenceH hSRS;
    char *pszResult = NULL;

    CPLErrorReset();
    
    hSRS = OSRNewSpatialReference( NULL );
    if( OSRSetFromUserInput( hSRS, pszUserInput ) == OGRERR_NONE )
        OSRExportToWkt( hSRS, &pszResult );
    else
    {
        CPLError( CE_Failure, CPLE_AppDefined,
                  "Translating source or target SRS failed:\n%s",
                  pszUserInput );
        exit( 1 );
    }
    
    OSRDestroySpatialReference( hSRS );

    return pszResult;
}
