#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265

void latlon2utm(double, double, int, double *, double *);

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:
  This program reprojects a geographic grid to UTM coordinates.  User must specify
  output domain via min/max Easting/Northing, and user must specify output cell size in
  meters.  All input pixels are convolved with a gaussian kernel of user-specified
  width; each output pixel value is the normalized sum (i.e. average) of all the
  convolved input kernels falling within the output pixel.

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
  double minx,maxx,miny,maxy;
  double sigma;
  float **data_out;
  float **count;
  float **kernel;
  float data;
  double x,y;
  double lat,lon;
  int tmp_int;
  int zone;
  int i_out,j_out;
  int ii_out,jj_out;
  int contribute;
  int width;
  int kernel_center;
  float kernel_a, kernel_b;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];

  /* Usage */
  if(argc!=13) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <zone> <resolution> <minx> <maxx> <miny> <maxy> <radius> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <zone>       UTM longitudinal zone\n");
    fprintf(stdout,"  <resolution> Output resolution (cellsize, in m)\n");
    fprintf(stdout,"  <minx>       Western boundary of output grid (m)\n");
    fprintf(stdout,"  <maxx>       Eastern boundary of output grid (m)\n");
    fprintf(stdout,"  <miny>       Southern boundary of output grid (m)\n");
    fprintf(stdout,"  <maxy>       Northern boundary of output grid (m)\n");
    fprintf(stdout,"  <radius>     Radius (i.e. sigma, or standard deviation) of gaussian kernel (pixels)\n");
    fprintf(stdout,"               NOTE: the gaussian kernel will act as a low-pass filter on the data,\n");
    fprintf(stdout,"               with cutoff wavelength equal to 2*PI*radius (pixels).\n");
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

  /* Read zone */
  zone = atoi(argv[3]);
  if (zone < 0) {
    fprintf(stderr,"%s: ERROR: Invalid UTM zone (%d)\n",argv[0],zone);
    exit(1);
  }

  /* Read output cellsize (resolution) */
  cellsize_out = atof(argv[4]);
  if (cellsize_out <= 0) {
    fprintf(stderr,"%s: ERROR: Invalid output cellsize (resolution) (%s)\n",argv[0],argv[4]);
    exit(1);
  }

  /* Read min/max lat/lon */
  minx = atof(argv[5]);
  maxx = atof(argv[6]);
  miny = atof(argv[7]);
  maxy = atof(argv[8]);

  /* Read sigma */
  sigma = atof(argv[9]);
  if (sigma <= 0) {
    fprintf(stderr,"%s: ERROR: Invalid kernel radius (%s)\n",argv[0],argv[9]);
    exit(1);
  }

  /* Read coord and data precision */
  coord_prec = atoi(argv[10]);
  data_prec = atoi(argv[11]);
  if (data_prec < 0) {
    fprintf(stderr,"%s: ERROR: Invalid data precision (%d)\n",argv[0],data_prec);
    exit(1);
  }

  /* Open the output file */
  if((fout=fopen(argv[12],"w"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[12]);
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

  /* Compute gaussian kernel */
  width = 2*((int)(3*sigma))+1;
  if ( (kernel = (float**)calloc(width,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for kernel array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<width; i++) {
    if ( (kernel[i] = (float*)calloc(width,sizeof(float))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for kernel array\n",argv[0]);
      exit(1);
    }
  }
  kernel_center = (width-1)/2;
  kernel_b = 1/(2*sigma*sigma);
  kernel_a = (1/PI)*kernel_b;
  for (ii=0; ii<width; ii++) {
    for (jj=0; jj<width; jj++) {
      kernel[ii][jj] = kernel_a * exp(-kernel_b * ( (ii-kernel_center)*(ii-kernel_center) + (jj-kernel_center)*(jj-kernel_center) ));
    }
  }

  /* Compute number of output rows and cols */
  ncols_out = (int)((maxx - minx)/cellsize_out);
  if ((maxx - minx)/cellsize_out > (double)ncols_out) ncols_out++;
  nrows_out = (int)((maxy - miny)/cellsize_out);
  if ((maxy - miny)/cellsize_out > (double)nrows_out) nrows_out++;

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
  if ( (count = (float**)calloc(nrows_out,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for count array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_out; i++) {
    if ( (count[i] = (float*)calloc(ncols_out,sizeof(float))) == NULL ) {
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
    lat = yllcorner + (nrows-0.5-i)*cellsize;
    for(j=0;j<ncols;j++) {
      lon = xllcorner + (j+0.5)*cellsize;
      latlon2utm(lat,lon,zone,&x,&y);
      i_out = nrows_out - 1 - (int)((y-miny)/cellsize_out);
      j_out = (int)((x-minx)/cellsize_out);
      if (int_data) {
        fscanf(fin,"%d",&tmp_int);
        data = (float)tmp_int;
      }
      else {
        fscanf(fin,"%f",&data);
      }
      contribute = 0;
      for (ii=0; ii < width; ii++) {
        ii_out = i_out+ii-kernel_center;
        for (jj=0; jj < width; jj++) {
          jj_out = j_out+jj-kernel_center;
          if (ii_out >= 0 && ii_out < nrows_out && jj_out >= 0 && jj_out < ncols_out) {
            contribute = 1;
          }
        }
      }
      if (contribute && ((int_data && data != (int)nodata) || data != nodata)) {
        for (ii=0; ii < width; ii++) {
          ii_out = i_out+ii-kernel_center;
          for (jj=0; jj < width; jj++) {
            jj_out = j_out+jj-kernel_center;
            if (ii_out >= 0 && ii_out < nrows_out && jj_out >= 0 && jj_out < ncols_out) {
              data_out[ii_out][jj_out] += data*kernel[ii][jj];
              count[ii_out][jj_out] += kernel[ii][jj];
            }
          }
        }
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

  /* Write the output header */
  fprintf(fout,"ncols %d\n",ncols_out);
  fprintf(fout,"nrows %d\n",nrows_out);
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",coord_prec);
  strcat(precstr,precstr2);
  sprintf(fmtstr,"xllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,minx);
  fprintf(fout,"\n");
  sprintf(fmtstr,"yllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,miny);
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

void latlon2utm(double lat, double lon, int zone, double *x, double *y) {

/******************************************************************************
* This function converts latitude and longitude (with UTM zone) to
* UTM Easting and Northing (assuming WGS84 datum)
*
* Inputs
* lat  : latitude (deg)
* lon  : longitude (deg)
* zone : Longitudinal UTM zone
*
* Outputs
* x    : Easting (m)
* y    : Northing (m)
******************************************************************************/

  // Datum-specific constants (assuming WGS84 datum)
  double x_offset = 500000; // False Easting (m)
  double k0 = 0.9996; // Scale factor
  double a = 6378137; // Equatorial radius of Earth (m)
  double b = 6356752.3142; // Polar radius of Earth (m)

  // Other terms
  double lon0; // Central meridian of zone (deg)
  double dlon; // Offset between input longitude and centeral meridian (1e4 sec)
  double phi;  // input latitude in radians
  double e; // eccentricity
  double eSq; // eccentricity squared
  double ePmSq; // e prime squared
  double n, n2, n3, n4, n5;
  double A0;
  double B0;
  double C0;
  double D0;
  double E0;
  double M;    // Meridional arc
  double tmp;
  double R1;   // Radius of curvature of Earth in meridional plane
  double N1;   // Radius of curvature of Earth perpendicular to meridional plane
  double tanphi;
  double sin1sec = 4.8481368e-6;
  double K1;
  double K2;
  double K3;
  double K4;
  double K5;

  // Compute derived constants
  lon0 = zone*6 - 183;
  dlon = (lon-lon0)*3600/1e4;
  phi = lat*PI/180;
  e = sqrt(1-(b*b)/(a*a));

  // Some useful intermediate terms
  eSq = e*e;
  ePmSq = (eSq)/(1-eSq);
  n = (a-b)/(a+b);
  n2 = n*n;
  n3 = n2*n;
  n4 = n3*n;
  n5 = n4*n;

  // Compute meridional arc
  A0 = a*(1 - n + (5/4)*(n2-n3) + (81/64)*(n4-n5));
  B0 = (3*a*n/2)*(1 - n + (7/8)*(n2-n3) + (55/64)*(n4-n5));
  C0 = (15*a*n2/16)*(1 - n + (3/4)*(n2-n3));
  D0 = (35*a*n3/48)*(1 - n + (11/16)*(n2-n3));
  E0 = (315*a*n4/51)*(1 - n);
  M = A0*phi - B0*sin(2*phi) + C0*sin(4*phi) - D0*sin(6*phi) + E0*sin(8*phi);

  // More intermediate terms
  tmp = sqrt(1-eSq*sin(phi)*sin(phi));
  R1 = a*(1-eSq)/(tmp*tmp*tmp);
  N1 = a/tmp;
  tanphi = sin(phi)/cos(phi);
//  tmp = (1/3600)*PI/180;
//  sin1sec = sin(tmp);
  K1 = M*k0;
  K2 = k0*sin1sec*sin1sec*N1*sin(phi)*cos(phi)*1e8/2;
  K3 = (k0*sin1sec*sin1sec*sin1sec*sin1sec*N1*sin(phi)*cos(phi)*cos(phi)*cos(phi)/24)*(5 - tanphi*tanphi - 9*ePmSq*cos(phi)*cos(phi) + 4*ePmSq*ePmSq*cos(phi)*cos(phi)*cos(phi)*cos(phi))*1e16;
  K4 = k0*sin1sec*N1*cos(phi)*1e4;
  K5 = (k0*sin1sec*sin1sec*sin1sec*N1*cos(phi)*cos(phi)*cos(phi)/6)*(1 - tanphi*tanphi + ePmSq*cos(phi)*cos(phi))*1e12;

  // Finally, x and y
  *y = K1 + K2*dlon*dlon + K3*dlon*dlon*dlon*dlon;
  *x = x_offset + K4*dlon + K5*dlon*dlon*dlon;

}
