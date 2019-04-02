#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define E_RADIUS 6371000.0 /* average radius of the earth in m*/    
#define PI 3.14159265

double  get_dist(double lat1, double long1, double lat2, double long2);

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:

  This program low-pass filters an arcinfo style ascii grid file by applying a
  moving window average to all the pixels.  There are multiple averaging methods,
  including a simple average of all values in the window, a truncated average,
  and a Gaussian smoother for which the user must specify the cutoff wavelength
  (in pixels).

  In the case of the Gaussian smoother, the window width is an automatic function
  if the cutoff wavelength.

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
  double xllcorner,yllcorner;
  double cellsize;
  float nodata;
  float **frac,**data,**data_out;
  int tmp_int;
  int count;
  int old_i,old_j;
  int width,width_j;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];
  float havg, hstd;
  float nodata_out = -9999;
  float data_mean, data_var, data_std;
  float frac_mean, frac_var, frac_std;
  float cov, corr;
  int window_count, count0, count1;
  int trunc, min_found, max_found;
  float data_mean_cover, total_frac;
  char method[10];
  int kernel_center_i,kernel_center_j;
  float sigma, kernel_a, kernel_b;
  float kernel_sum;
  float **kernel;
  float **window_data;
  int tmp_i, tmp_j;
  int window_min, window_max;
  int first;
  int geog;
  float length;
  double maxlat,dy,dx;
  float aspect_ratio;

  /* Usage */
  if(argc!=10) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <method> <length> <geog> <trunc> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <method>     Smoothing method (\"mean\",\"gauss\")\n");
    fprintf(stdout,"  <length>     For \"mean\", length = width (in pixels) of smoothing window\n");
    fprintf(stdout,"               For \"gauss\", length = cut-off wavelength (in pixels);\n");
    fprintf(stdout,"                 in this case, smoothing window width will be set to:\n");
    fprintf(stdout,"                   2*((int)(3*sigma))+1\n");
    fprintf(stdout,"                 where\n");
    fprintf(stdout,"                   sigma = length/(2*pi)\n");
    fprintf(stdout,"                         = the radius of the gaussian kernel's inflection point\n");
    fprintf(stdout,"                 so that the window contains 3*sigma on each side of the central pixel\n");
    fprintf(stdout,"  <geog>       1 = rows and columns are latitude and longitude, respectively (unequal in size), and the geographical aspect ratio will be taken into account; 0 = rows and columns are equal in size\n");
    fprintf(stdout,"  <trunc>      Number of extreme values in the window to ignore on each side of the distribution (generally set this to 0 unless you want a truncated mean)\n");
    fprintf(stdout,"  <coord_prec> Precision of coordinates, i.e. number of decimal places\n");
    fprintf(stdout,"  <data_prec>  Precision of data, i.e. number of decimal places (ignored for \"int\" data)\n");
    fprintf(stdout,"  <out_grid>   Output grid file name\n");
    exit(0);
  }

  /* Open the input grid file */
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

  /* Read method */
  strcpy(method,argv[3]);
  if (!strcmp(method,"mean") || !strcmp(method,"gauss")) {
    /* nothing to do */
  }
  else {
    fprintf(stderr,"%s: ERROR: method (%s) is not a supported method\n",argv[0],method);
    exit(1);
  }

  /* Read length */
  length = atof(argv[4]);

  /* Read geog flag */
  geog = atoi(argv[5]);
  if (geog != 0) geog = 1;

  /* Read number of pixels to truncate */
  trunc = atoi(argv[6]);
  if (trunc < 0) {
    fprintf(stderr,"%s: ERROR: number of values to ignore (%s) must not be negative\n",argv[0],argv[6]);
    exit(1);
  }

  /* Read coord and data precision */
  coord_prec = atoi(argv[7]);
  data_prec = atoi(argv[8]);
  if (data_prec < 0) {
    fprintf(stderr,"%s: ERROR: Invalid data precision (%d)\n",argv[0],data_prec);
    exit(1);
  }

  /* Open the output file */
  if((fout=fopen(argv[9],"w"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[9]);
    exit(1);
  }

  if (!strcmp(method,"mean")) {
    width = (int)length;
    if (width <= 0) {
      fprintf(stderr,"%s: ERROR: window width (%s) is not a positive integer\n",argv[0],argv[3]);
      exit(1);
    }
  }
  else if (!strcmp(method,"gauss")) {
    sigma = length/(2*PI);
    width = 2*((int)(3*sigma))+1;
  }

  /* Read the input file header */
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

  /* Compute maximum aspect ratio */
  if (geog) {
    maxlat = yllcorner + nrows*cellsize;
    if (fabs(maxlat) < fabs(yllcorner)) maxlat = yllcorner;
    dy = get_dist(maxlat,xllcorner,maxlat+cellsize,xllcorner);
    dx = get_dist(maxlat,xllcorner,maxlat,xllcorner+cellsize);
    aspect_ratio = dy/dx;
  }
  else {
    aspect_ratio = 1;
  }
  width_j = (int)((float)width*aspect_ratio);

  /* Allocate data arrays */
  if ( (data = (float**)calloc(nrows,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows; i++) {
    if ( (data[i] = (float*)calloc(ncols,sizeof(float))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data array\n",argv[0]);
      exit(1);
    }
  }
  if ( (data_out = (float**)calloc(nrows,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data_out array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows; i++) {
    if ( (data_out[i] = (float*)calloc(ncols,sizeof(float))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data_out array\n",argv[0]);
      exit(1);
    }
  }
  if ( (window_data = (float**)calloc(width,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for window_data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<width; i++) {
    if ( (window_data[i] = (float*)calloc(width_j,sizeof(float))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for window_data array\n",argv[0]);
      exit(1);
    }
  }

  /* Create filter kernel, if appropriate */
  kernel_center_i = (width-1)/2;
  kernel_center_j = (width_j-1)/2;
  if (!strcmp(method,"gauss")) {
    if ( (kernel = (float**)calloc(width,sizeof(float*))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for kernel array\n",argv[0]);
      exit(1);
    }
    for (i=0; i<width; i++) {
      if ( (kernel[i] = (float*)calloc(width_j,sizeof(float))) == NULL ) {
        fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for kernel array\n",argv[0]);
        exit(1);
      }
    }
    kernel_b = 1/(2*sigma*sigma);
    kernel_a = (1/PI)*kernel_b;
    if (!geog) {
      for (ii=0; ii<width; ii++) {
        for (jj=0; jj<width_j; jj++) {
          kernel[ii][jj] = kernel_a * exp(-kernel_b * ( (ii-kernel_center_i)*(ii-kernel_center_i) + (jj-kernel_center_j)*(jj-kernel_center_j)/(aspect_ratio*aspect_ratio) ));
        }
      }
    }
  }

  /* Read in the data */
  for(i=0;i<nrows;i++) {
    for(j=0;j<ncols;j++) {
      if (int_data) {
        fscanf(fin,"%d",&tmp_int);
        data[i][j] = (float)tmp_int;
      }
      else {
        fscanf(fin,"%f",&data[i][j]);
      }
    }
  }
  fclose(fin); 

  /* Apply filter to moving window */
  for(i=0;i<nrows;i++) {

    if (!strcmp(method,"gauss") && geog) {
      // Recompute kernel based on local aspect ratio
      maxlat = yllcorner + (nrows-1-i)*cellsize;
      dx = get_dist(maxlat,xllcorner,maxlat,xllcorner+cellsize);
      aspect_ratio = dy/dx;
      for (ii=0; ii<width; ii++) {
        for (jj=0; jj<width_j; jj++) {
          kernel[ii][jj] = kernel_a * exp(-kernel_b * ( (ii-kernel_center_i)*(ii-kernel_center_i) + (jj-kernel_center_j)*(jj-kernel_center_j)/(aspect_ratio*aspect_ratio) ));
        }
      }
    }

    for(j=0;j<ncols;j++) {

      /* Gather window data */
      for (ii=0; ii<width; ii++) {
        tmp_i = i-kernel_center_i+ii;
        for (jj=0; jj<width_j; jj++) {
          tmp_j = j-kernel_center_j+jj;
          if (tmp_i > 0 && tmp_j > 0 && tmp_i < nrows && tmp_j < ncols && data[tmp_i][tmp_j] != nodata) {
            window_data[ii][jj] = data[tmp_i][tmp_j];
          }
          else {
            window_data[ii][jj] = nodata;
          }
        }
      }

      /* Set extremes to nodata if appropriate */
      /* Note: for large windows, this might be done better by sorting the data in the window */
      if (trunc) {
        for (k=0; k<trunc; k++) {
          first = 1;
          for (ii=0; ii<width; ii++) {
            for (jj=0; jj<width_j; jj++) {
              if (window_data[ii][jj] != nodata) {
                if (first) {
                  window_min = window_data[ii][jj];
                  window_max = window_data[ii][jj];
                  first = 0;
                }
                else {
                  window_min = (window_min > window_data[ii][jj]) ? window_data[ii][jj] : window_min;
                  window_max = (window_max > window_data[ii][jj]) ? window_max : window_data[ii][jj];
                }
              }
            }
          }
          min_found = 0;
          max_found = 0;
          for (ii=0; ii<width; ii++) {
            for (jj=0; jj<width_j; jj++) {
              if (window_data[ii][jj] == window_min && !min_found) {
                window_data[ii][jj] = nodata;
                min_found = 1;
              }
              else if (window_data[ii][jj] == window_max && !max_found) {
                window_data[ii][jj] = nodata;
                max_found = 1;
              }
            }
          }
        }
      }

      /* Now apply filter to window */
      /* If we ever add a median method, we'll have to sort the data in the window */
      if (!strcmp(method,"mean")) {
        data_out[i][j] = 0;
        count = 0;
        for (ii=0; ii<width; ii++) {
          for (jj=0; jj<width_j; jj++) {
            if (window_data[ii][jj] != nodata) {
              data_out[i][j] += window_data[ii][jj];
              count++;
            }
          }
        }
        if (count > 0) {
          data_out[i][j] /= count;
        }
        else {
          /* If no valid data, we have to skip this window */
          data_out[i][j] = nodata_out;
          continue;
        }
      }

      else if (!strcmp(method,"gauss")) {
        data_out[i][j] = 0;
        kernel_sum = 0;
        for (ii=0; ii<width; ii++) {
          for (jj=0; jj<width_j; jj++) {
            if (window_data[ii][jj] != nodata) {
              data_out[i][j] += kernel[ii][jj]*window_data[ii][jj];
              kernel_sum += kernel[ii][jj];
            }
          }
        }
        if (kernel_sum > 0) {
          data_out[i][j] /= kernel_sum;
        }
        else {
          /* If no valid data, we have to skip this window */
          data_out[i][j] = nodata_out;
          continue;
        }
      }
    }
  }

  /* Write the output header */
  fprintf(fout,"ncols %d\n",ncols);
  fprintf(fout,"nrows %d\n",nrows);
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",coord_prec);
  strcat(precstr,precstr2);
  sprintf(fmtstr,"xllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,xllcorner);
  fprintf(fout,"\n");
  sprintf(fmtstr,"yllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,yllcorner);
  fprintf(fout,"\n");
  sprintf(fmtstr,"cellsize ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,cellsize);
  fprintf(fout,"\n");
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  strcat(precstr,precstr2);
  sprintf(fmtstr,"nodata ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,nodata_out);
  fprintf(fout,"\n");

  /* Write the output data */
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  strcat(precstr,precstr2);
  strcpy(fmtstr,precstr);
  for (i=0; i<nrows; i++) {
    for (j=0; j<ncols; j++) {
      fprintf(fout,fmtstr,data_out[i][j]);
      if (j < ncols-1) {
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

double  get_dist(double lat1, double long1, double lat2, double long2)
{
  double theta1;
  double phi1;
  double theta2;
  double phi2;
  double dtor;
  double term1;
  double term2;
  double term3;
  double temp;
  double dist;

  dtor = 2.0*PI/360.0;
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  temp = (double) (1.0 < temp) ? 1.0 : temp;
  dist = E_RADIUS*acos(temp);/*area between lat & lon*/

  return dist;
}  
