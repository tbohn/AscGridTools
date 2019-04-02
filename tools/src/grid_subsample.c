#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:
  This program sub-samples a grid to a higher resolution (smaller cellsize).
  The user can choose which method to use: nearest neighbor, moving average,
  or gaussian smoothing.  Output grid will cover the same domain as the input
  grid.

  Modifications:

******************************************************************************/
{
  FILE *fin,*fout;
  int i,j,k;
  int ii,jj;
  char argstr[10];
  int int_data;
  int int_data_out;
  int coord_prec;
  int data_prec;
  int nrows,ncols;
  int nrows_out,ncols_out;
  double xllcorner,yllcorner;
  double cellsize;
  float nodata;
  double cellsize_out;
  float length;
  float ratio;
  int res_ratio;
  int width;
  char method[10];
  int kernel_center;
  float distance;
  float data_tmp;
  float sigma, kernel_a, kernel_b;
  float **kernel;
  float **data;
  float **data_out;
  float **count;
  double x,y;
  int i_out,j_out;
  int ii_out,jj_out;
  int contribute;
  int tmp_int, tmp_count, tmp_sum;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];

  /* Usage */
  if(argc!=10) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <cellsize_out> <method> <length> <type_out> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <cellsize_out> Output cell width, in input units\n");
    fprintf(stdout,"  <method>     Sub-sampling method; \"nn\" = nearest neighbor; \"mean\" = moving average; \"bilin\" = bi-linear interpolation, \"gauss\" = gaussian smoothing.\n");
    fprintf(stdout,"  <length>     Characteristic length (in output pixels) of averaging method:\n");
    fprintf(stdout,"                 For method=\"nn\" or \"bilin\", length is ignored (can be 0)\n");
    fprintf(stdout,"                 For method=\"mean\", length=width of averaging window\n");
    fprintf(stdout,"                 For method=\"gauss\", length=low-pass filter cut-off wavelength;\n");
    fprintf(stdout,"                   in this case, sigma, the radius of the gaussian inflection point, is:\n");
    fprintf(stdout,"                     sigma = length/(2*PI)\n");
    fprintf(stdout,"                   and the smoothing window width will be set to:\n");
    fprintf(stdout,"                     2*((int)(3*sigma))+1\n");
    fprintf(stdout,"                   so that the window contains 3*sigma on each side of the central pixel\n");
    fprintf(stdout,"  <type_out>   Output data type (\"int\" or \"float\")\n");
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

  /* Read output cellsize */
  cellsize_out = atof(argv[3]);
  if (cellsize_out <= 0) {
    fprintf(stderr,"%s: ERROR: Invalid output cellsize (%s)\n",argv[0],argv[3]);
    exit(1);
  }

  /* Read method */
  strcpy(method,argv[4]);
  if (!strcmp(method,"nn") || !strcmp(method,"mean") || !strcmp(method,"bilin") || !strcmp(method,"gauss")) {
    /* nothing to do */
  }
  else {
    fprintf(stderr,"%s: ERROR: method (%s) is not a supported method\n",argv[0],method);
    exit(1);
  }

  /* Read length */
  length = atof(argv[5]);
  if (!strcmp(method,"nn")) {
    length = 1;
  }

  /* Read output data type */
  strcpy(argstr,argv[6]);
  int_data_out = 0;
  if(!strcasecmp(argstr,"int"))   {
    int_data_out = 1;
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
  ratio = cellsize/cellsize_out;
  nrows_out = (int)((float)nrows*ratio);
  ncols_out = (int)((float)ncols*ratio);
  if (nrows/ratio > (double)nrows_out) nrows_out++;
  if (ncols/ratio > (double)ncols_out) ncols_out++;

  /* Allocate data array */
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

  /* Create filter kernel */
  if (!strcmp(method,"gauss")) {
    sigma = length/(2*PI);
    width = 2*((int)(3*sigma))+1;
  }
  else if (!strcmp(method,"bilin")) {
    width = 2*((int)(cellsize/cellsize_out + 0.5))+1;
  }
  else {
    width = length;
  }
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
  if (!strcmp(method,"gauss")) {
    kernel_center = (width-1)/2;
    kernel_b = 1/(2*sigma*sigma);
    kernel_a = (1/PI)*kernel_b;
    for (ii=0; ii<width; ii++) {
      for (jj=0; jj<width; jj++) {
        distance = sqrt((ii-kernel_center)*(ii-kernel_center) + (jj-kernel_center)*(jj-kernel_center));
        kernel[ii][jj] = kernel_a * exp(-kernel_b * ( distance * distance ));
      }
    }
  }
  else if (!strcmp(method,"bilin")) {
    kernel_center = (width-1)/2;
    for (ii=0; ii<width; ii++) {
      for (jj=0; jj<width; jj++) {
        distance = sqrt((ii-kernel_center)*(ii-kernel_center) + (jj-kernel_center)*(jj-kernel_center));
        kernel[ii][jj] = 1 - distance/kernel_center;
        if (kernel[ii][jj] < 0) {
          kernel[ii][jj] = 0;
        }
      }
    }
  }
  else { /* mean */
    kernel_center = (width-1)/2;
    for (ii=0; ii<width; ii++) {
      for (jj=0; jj<width; jj++) {
        kernel[ii][jj] = 1;
      }
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
  if (!strcmp(method,"gauss") || !strcmp(method,"mean") || !strcmp(method,"bilin")) {
    for(i=0;i<nrows;i++) {
      y = yllcorner + (nrows-0.5-i)*cellsize;
      for(j=0;j<ncols;j++) {
        x = xllcorner + (j+0.5)*cellsize;
        i_out = nrows_out - 1 - (int)((y-yllcorner)/cellsize_out);
        j_out = (int)((x-xllcorner)/cellsize_out);
        if (int_data) {
          fscanf(fin,"%d",&tmp_int);
          data_tmp = (float)tmp_int;
        }
        else {
          fscanf(fin,"%f",&data_tmp);
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
        if (contribute && ((int_data && data_tmp != (int)nodata) || data_tmp != nodata)) {
          for (ii=0; ii < width; ii++) {
            ii_out = i_out+ii-kernel_center;
            for (jj=0; jj < width; jj++) {
              jj_out = j_out+jj-kernel_center;
              if (ii_out >= 0 && ii_out < nrows_out && jj_out >= 0 && jj_out < ncols_out) {
                data_out[ii_out][jj_out] += data_tmp*kernel[ii][jj];
                count[ii_out][jj_out] += kernel[ii][jj];
              }
            }
          }
        }
      }
    }

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

  }
  else {

    for(i=0;i<nrows;i++) {
      for(j=0;j<ncols;j++) {
        if (int_data) {
          fscanf(fin,"%d",&tmp_int);
          data[i][j] = (float)tmp_int;
        }
        else {
          fscanf(fin,"%f",&(data[i][j]));
        }
      }
    }
    for(i_out=0;i_out<nrows_out;i_out++) {
      y = yllcorner + (nrows_out-0.5-i_out)*cellsize_out;
      i = nrows-1 - (int)((y-yllcorner)/cellsize);
      for(j_out=0;j_out<ncols_out;j_out++) {
        x = xllcorner + j_out*cellsize_out;
        j = (int)((x-xllcorner)/cellsize);
        data_out[i_out][j_out] = data[i][j];
      }
    }

  }
  fclose(fin);

  if (length > 1 && width > 1) {

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

  }

  /* Write the output header */
  fprintf(fout,"ncols %d\n",ncols_out);
  fprintf(fout,"nrows %d\n",nrows_out);
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
  fprintf(fout,fmtstr,cellsize_out);
  fprintf(fout,"\n");
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  if (int_data_out) {
    fprintf(fout,"nodata %d\n",(int)nodata);
  }
  else {
    strcat(precstr,precstr2);
    sprintf(fmtstr,"nodata ");
    strcat(fmtstr,precstr);
    fprintf(fout,fmtstr,nodata);
    fprintf(fout,"\n");
  }

  /* Write the output data */
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  strcat(precstr,precstr2);
  strcpy(fmtstr,precstr);
  for (i=0; i<nrows_out; i++) {
    for (j=0; j<ncols_out; j++) {
      if (int_data_out) {
        fprintf(fout,"%d",(int)data_out[i][j]);
      }
      else {
        fprintf(fout,fmtstr,data_out[i][j]);
      }
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

