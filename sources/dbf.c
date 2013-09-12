
#include <stdlib.h>
#include <string.h>
#include "shapefil.h"


void createDBF(const char *DBFfilename, int id[], int nsamp, double Lon[], double Lat[],\
               double Elevation[], char **SampName, double Temperature[], double ErosionRate[], \
	       double Geotherm[], double ExhumationRate[])
{

    DBFHandle	hDBF;
    int		i;
    int         iRecord;
    char* test;


    hDBF = DBFCreate( DBFfilename );

    // Warning: 10 character limit for field names imposed by dbf standard
    DBFAddField( hDBF , "SName" , FTString , 40 , 0 );
    DBFAddField( hDBF , "SID."  , FTInteger, 10 , 0 );
    DBFAddField( hDBF , "Lon."  , FTDouble , 10 , 2 );
    DBFAddField( hDBF , "Lat."  , FTDouble , 10 , 2 );
    DBFAddField( hDBF , "Elev." , FTDouble , 10 , 1 );
    DBFAddField( hDBF , "Temp." , FTDouble , 10 , 1 );
    DBFAddField( hDBF , "E.Rate" , FTDouble , 10 , 1 );
    DBFAddField( hDBF , "Geoth." , FTDouble , 10 , 1 );
    DBFAddField( hDBF , "Ex.Rate" , FTDouble , 10 , 1 );

    for(i=0; i < nsamp; i++)
    {
    iRecord = DBFGetRecordCount( hDBF );

    DBFWriteStringAttribute(hDBF, iRecord, 0, SampName[i] );
    DBFWriteIntegerAttribute(hDBF, iRecord, 1, id[i] );
    DBFWriteDoubleAttribute(hDBF, iRecord, 2, Lon[i]);
    DBFWriteDoubleAttribute(hDBF, iRecord, 3, Lat[i]);
    DBFWriteDoubleAttribute(hDBF, iRecord, 4, Elevation[i]);
    DBFWriteDoubleAttribute(hDBF, iRecord, 5, Temperature[i]);
    DBFWriteDoubleAttribute(hDBF, iRecord, 6, ErosionRate[i]);
    DBFWriteDoubleAttribute(hDBF, iRecord, 7, Geotherm[i]);
    DBFWriteDoubleAttribute(hDBF, iRecord, 8, ExhumationRate[i]);
    }

    DBFClose( hDBF );

}

