#include "ogr_api.h"
#include "ogr_srs_api.h"

void GetDistanceBetweenPoints(double *X1, double *Y1, double *X2, double *Y2, double *distance)
      {
      int zone;
      //printf("%f %f %f %f\n", *X1, *Y1, *X2, *Y2);

      OGRGeometryH pt1, pt2;
      OGRSpatialReferenceH SREF;
      OGRSpatialReferenceH UTM;

      SREF = OSRNewSpatialReference(NULL);
      OSRImportFromEPSG(SREF, 4326);

      zone = (int)(((*X1 + 180.0) / 6.0) + 1.0);
      //printf("Zone: %i\n", zone);
      UTM = OSRNewSpatialReference(NULL);
      if(*Y1 < 0) OSRSetUTM(UTM,zone,FALSE);
      if(*Y1 > 0) OSRSetUTM(UTM,zone,TRUE);

      pt1 = OGR_G_CreateGeometry(wkbPoint);
      pt2 = OGR_G_CreateGeometry(wkbPoint);

      OGR_G_SetPoint_2D(pt1, 0, *X1, *Y1);
      OGR_G_SetPoint_2D(pt2, 0, *X2, *Y2);

      OGR_G_AssignSpatialReference(pt1, SREF);
      OGR_G_AssignSpatialReference(pt2, SREF);

      // Convert to UTM
      OGR_G_TransformTo(pt1,UTM);
      OGR_G_TransformTo(pt2,UTM);


      *distance = OGR_G_Distance(pt1, pt2);

      }

void LatLon2UTM(double *X1, double *Y1, double *X2, double *Y2)
      {
      int zone;

      OGRGeometryH pt1;
      OGRSpatialReferenceH SREF;
      OGRSpatialReferenceH UTM;

      SREF = OSRNewSpatialReference(NULL);
      OSRImportFromEPSG(SREF, 4326);

      zone = (int)(((*X1 + 180.0) / 6.0) + 1.0);
      UTM = OSRNewSpatialReference(NULL);
      if(*Y1 < 0) OSRSetUTM(UTM,zone,FALSE);
      if(*Y1 > 0) OSRSetUTM(UTM,zone,TRUE);

      pt1 = OGR_G_CreateGeometry(wkbPoint);

      OGR_G_SetPoint_2D(pt1, 0, *X1, *Y1);

      OGR_G_AssignSpatialReference(pt1, SREF);

      // Convert to UTM
      OGR_G_TransformTo(pt1,UTM);

      *X2 = OGR_G_GetX(pt1, 0);
      *Y2 = OGR_G_GetY(pt1, 0);
      }

