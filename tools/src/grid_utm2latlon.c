#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void utm2latlon(double, double, int, double *, double *);

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:
  This program reprojects a UTM grid to geographic coordinates.  User must specify
  output domain via min/max lat/lon, and user must specify output cell size in
  degrees.  All UTM pixels falling within an output pixel are averaged together.

  Modifications:

******************************************************************************/
{
  FILE *fin,*fout;
  int i,j,k;
  int ii,jj;
  char argstr[10];
  int int_data;
  int coord_prec;
  int data_prec;
  int nrows,ncols;
  int nrows_out,ncols_out;
  double xllcorner,yllcorner;
  double cellsize;
  float nodata;
  double cellsize_out;
  double minlon,maxlon,minlat,maxlat;
  float **data_out;
  int **count;
  float data;
  double x,y;
  double lat,lon;
  int tmp_int, tmp_count, tmp_sum;
  int zone;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];

  /* Usage */
  if(argc!=12) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <zone> <resolution> <minlon> <maxlon> <minlat> <maxlat> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <zone>       UTM longitudinal zone\n");
    fprintf(stdout,"  <resolution> Output resolution (cellsize, in degrees)\n");
    fprintf(stdout,"  <minlon>     Western boundary of output grid (degrees)\n");
    fprintf(stdout,"  <maxlon>     Eastern boundary of output grid (degrees)\n");
    fprintf(stdout,"  <minlat>     Southern boundary of output grid (degrees)\n");
    fprintf(stdout,"  <maxlat>     Northern boundary of output grid (degrees)\n");
    fprintf(stdout,"  <coord_prec> Precision of coordinates, i.e. number of decimal places\n");
    fprintf(stdout,"  <data_prec>  Precision of data, i.e. number of decimal places (ignored for \"int\" data)\n");
    fprintf(stdout,"  <out_grid>   Output grid file name\n");
    exit(0);
  }

  /* Open the input file */
  if((fin=fopen(argv[1],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[1]);
    exit(1);
  }

  /* Read data type */
  strcpy(argstr,argv[2]);
  int_data = 0;
  if(!strcasecmp(argstr,"int"))   {
    int_data = 1;
  }

  /* Read coord and data precision */
  zone = atoi(argv[3]);
  if (zone < 0) {
    fprintf(stderr,"%s: ERROR: Invalid UTM zone (%d)\n",argv[0],zone);
    exit(1);
  }

  /* Read output cellsize */
  cellsize_out = atof(argv[4]);
  if (cellsize_out <= 0) {
    fprintf(stderr,"%s: ERROR: Invalid output cellsize (%s)\n",argv[0],argv[4]);
    exit(1);
  }

  /* Read min/max lat/lon */
  minlon = atof(argv[5]);
  maxlon = atof(argv[6]);
  minlat = atof(argv[7]);
  maxlat = atof(argv[8]);

  /* Read coord and data precision */
  coord_prec = atoi(argv[9]);
  data_prec = atoi(argv[10]);
  if (data_prec < 0) {
    fprintf(stderr,"%s: ERROR: Invalid data precision (%d)\n",argv[0],data_prec);
    exit(1);
  }

  /* Open the output file */
  if((fout=fopen(argv[11],"w"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[11]);
    exit(1);
  }

  /* Read in the header */
  fscanf(fin,"%*s %d",&ncols);
  fscanf(fin,"%*s %d",&nrows);
  fscanf(fin,"%*s %s",tmpstr2);
  xllcorner = atof(tmpstr2);
  fscanf(fin,"%*s %s",tmpstr2);
  yllcorner = atof(tmpstr2);
  fscanf(fin,"%*s %s",tmpstr2);
  cellsize = atof(tmpstr2);
  if (int_data) {
    fscanf(fin,"%*s %d",&tmp_int);
    nodata = (float)tmp_int;
  }
  else {
    fscanf(fin,"%*s %f",&nodata);
  }

  /* Compute number of output rows and cols */
  ncols_out = (int)((maxlon - minlon)/cellsize_out);
  if ((maxlon - minlon)/cellsize_out > (double)ncols_out) ncols_out++;
  nrows_out = (int)((maxlat - minlat)/cellsize_out);
  if ((maxlat - minlat)/cellsize_out > (double)nrows_out) nrows_out++;

  /* Allocate output data array */
  if ( (data_out = (float**)calloc(nrows_out,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data_out array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_out; i++) {
    if ( (data_out[i] = (float*)calloc(ncols_out,sizeof(float))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data_out array\n",argv[0]);
      exit(1);
    }
  }

  /* Allocate count array */
  if ( (count = (int**)calloc(nrows_out,sizeof(int*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for count array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_out; i++) {
    if ( (count[i] = (int*)calloc(ncols_out,sizeof(int))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for count array\n",argv[0]);
      exit(1);
    }
  }

  /* Initialize output data arrays */
  for(i=0;i<nrows_out;i++) {
    for(j=0;j<ncols_out;j++) {
      data_out[i][j] = 0;
      count[i][j] = 0;
    }
  }

  /* Read in the data */
  for(i=0;i<nrows;i++) {
    y = yllcorner + (nrows-0.5-i)*cellsize;
    for(j=0;j<ncols;j++) {
      x = xllcorner + (j+0.5)*cellsize;
      utm2latlon(x,y,zone,&lat,&lon);
      ii = nrows_out - 1 - (int)((lat-minlat)/cellsize_out);
      jj = (int)((lon-minlon)/cellsize_out);
      if (int_data) {
        fscanf(fin,"%d",&tmp_int);
        data = (float)tmp_int;
      }
      else {
        fscanf(fin,"%f",&data);
      }
      if (ii >= 0 && ii < nrows_out && jj >= 0 && jj < ncols_out
        && ((int_data && data != (int)nodata) || data != nodata)) {
        data_out[ii][jj] += data;
        count[ii][jj]++;
      }
    }
  }
  fclose(fin); 

  /* Finish averaging */
  for(i=0;i<nrows_out;i++) {
    for(j=0;j<ncols_out;j++) {
      if (count[i][j] > 0) {
        data_out[i][j] /= count[i][j];
      }
      else {
        data_out[i][j] = nodata;
      }
    }
  }

  /* Fill holes */
  for(i=0;i<nrows_out;i++) {
    for(j=0;j<ncols_out;j++) {
      if (data_out[i][j] == nodata) {
        data_out[i][j] = 0;
        tmp_count = 0;
        for (ii=-1; ii<=1; ii++) {
          for (jj=-1; jj<=1; jj++) {
            if (i+ii >= 0 && i+ii < nrows_out && j+jj >= 0 && j+jj < ncols_out && (ii != 0 || jj != 0) && data_out[i+ii][j+jj] != nodata) {
              tmp_sum += data_out[i+ii][j+jj];
              tmp_count++;
            }
          }
        }
        if (tmp_count >= 4) {
          data_out[i][j] /= tmp_count;
        }
        else {
          data_out[i][j] = nodata;
        }
      }
    }
  }

  /* Write the output header */
  fprintf(fout,"ncols %d\n",ncols_out);
  fprintf(fout,"nrows %d\n",nrows_out);
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",coord_prec);
  strcat(precstr,precstr2);
  sprintf(fmtstr,"xllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,minlon);
  fprintf(fout,"\n");
  sprintf(fmtstr,"yllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,minlat);
  fprintf(fout,"\n");
  sprintf(fmtstr,"cellsize ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,cellsize_out);
  fprintf(fout,"\n");
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  strcat(precstr,precstr2);
  sprintf(fmtstr,"nodata ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,nodata);
  fprintf(fout,"\n");

  /* Write the output data */
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  strcat(precstr,precstr2);
  strcpy(fmtstr,precstr);
  for (i=0; i<nrows_out; i++) {
    for (j=0; j<ncols_out; j++) {
      fprintf(fout,fmtstr,data_out[i][j]);
      if (j < ncols_out-1) {
        fprintf(fout," ");
      }
      else {
        fprintf(fout,"\n");
      }
    }
  }
  fclose(fout);

  return 0;
  
}

void utm2latlon(double x, double y, int zone, double *lat, double *lon) {

/******************************************************************************
* This function converts UTM Easting and Northing (with zone) to
* latitude and longitude (assuming WGS84 datum)
*
* Inputs
* x    : Easting (m)
* y    : Northing (m)
* zone : Longitudinal UTM zone
*
* Outputs
* lat  : latitude (deg)
* lon  : longitude (deg)
******************************************************************************/

  // Constants
  double PI = 3.1415927;

  // Datum-specific constants (assuming WGS84 datum)
  double x_offset = 500000; // False Easting (m)
  double k0 = 0.9996; // Scale factor
  double a = 6378137; // Equatorial radius of Earth (m)
  double b = 6356752.3142; // Polar radius of Earth (m)

  // Other terms
  double e; // eccentricity
  double eSq; // eccentricity squared
  double ePmSq; // e prime squared
  double e1;
  double J1;
  double J2;
  double J3;
  double J4;
  double lon0; // Central meridian of zone (deg)
  double M;    // Meridional arc
  double mu;
  double fp;   // Footprint latitude
  double C1;
  double tanfp;
  double T1;
  double tmp;
  double R1;   // Radius of curvature of Earth in meridional plane
  double N1;   // Radius of curvature of Earth perpendicular to meridional plane
  double D;
  double Q1;
  double Q2;
  double Q3;
  double Q4;
  double Q5;
  double Q6;
  double Q7;
  double phi;  // Output latitude in radians
  double dlon; // Offset between output longitude and centeral meridian

  // Compute derived constants
  lon0 = zone*6 - 183;
  e = sqrt(1-(b*b)/(a*a));

  // Some useful intermediate terms
  eSq = e*e;
  ePmSq = (eSq)/(1-eSq);

  // Compute footprint latitude
  e1 = (1 - sqrt(1-eSq))/(1 + sqrt(1-eSq));
  J1 = (3*e1/2 - 27*e1*e1*e1/32);
  J2 = (21*e1*e1/16 - 55*e1*e1*e1*e1/256);
  J3 = (151*e1*e1*e1/96);
  J4 = (1097*e1*e1*e1*e1/512);
  M = y/k0;
  mu = M/(a*(1-eSq/4 - 3*eSq*eSq/64 - 5*eSq*eSq*eSq/256));
  fp = mu + J1*sin(2*mu) + J2*sin(4*mu) + J3*sin(6*mu) + J4*sin(8*mu);

  // More intermediate terms
  C1 = ePmSq*cos(fp)*cos(fp);
  tanfp = sin(fp)/cos(fp);
  T1 = tanfp*tanfp;
  tmp = sqrt(1-eSq*sin(fp)*sin(fp));
  R1 = a*(1-eSq)/(tmp*tmp*tmp);
  N1 = a/tmp;
  D = (x-x_offset)/(N1*k0);
  Q1 = N1*tanfp/R1;
  Q2 = D*D/2;
  Q3 = (5 + 3*T1 + 10*C1 - 4*C1*C1 - 9*ePmSq)*D*D*D*D/24;
  Q4 = (61 + 90*T1 + 298*C1 + 45*T1*T1 - 3*C1*C1 - 252*ePmSq)*D*D*D*D*D*D/720;
  Q5 = D;
  Q6 = (1 + 2*T1 + C1)*D*D*D/6;
  Q7 = (5 - 2*C1 + 28*T1 - 3*C1*C1 + 8*ePmSq + 24*T1*T1)*D*D*D*D*D/120;

  // Finally, lat and lon
  phi = fp - Q1*(Q2 - Q3 + Q4);
  *lat = phi*180/PI;
  dlon = (Q5 - Q6 + Q7)/cos(fp);
  dlon *=180/PI;
  *lon = lon0 + dlon;

}
