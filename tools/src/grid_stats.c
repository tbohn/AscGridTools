#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:

  This program computes statistics (e.g., mean, standard deviation, etc.) of an
  arcinfo style ascii grid file, at a coarser resolution.

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
  float **data,**data_out;
  float mean, var, min, max, sum;
  int tmp_int;
  int count;
  int old_i,old_j;
  int width;
  char statstr[5];
  char precstr[10];
  char tmpstr[500];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];

  /* Usage */
  if(argc!=8) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <stat> <width> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <stat>       Statistic to compute (\"mean\", \"var\", \"std\", \"min\", \"max\", \"sum\")\n");
    fprintf(stdout,"  <width>      Width of analysis window; this is the resolution at which statistics will be output; a value of 0 == window encompasses entire grid\n");
    fprintf(stdout,"  <coord_prec> Precision of coordinates, i.e. number of decimal places\n");
    fprintf(stdout,"  <data_prec>  Precision of data, i.e. number of decimal places (ignored for \"int\" data)\n");
    fprintf(stdout,"  <out_grid>   Output grid file prefix; \".mean.asc\", \".std.asc\", etc. will be appended to this prefix to build the various output filenames\n");
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

  /* Read stat type */
  strcpy(statstr,argv[3]);

  /* Read width */
  width = atoi(argv[4]);
  if (width < 0) {
    fprintf(stderr,"%s: ERROR: window width (%s) must be a positive integer or 0\n",argv[0],argv[4]);
    exit(1);
  }

  /* Read coord and data precision */
  coord_prec = atoi(argv[5]);
  data_prec = atoi(argv[6]);
  if (data_prec < 0) {
    fprintf(stderr,"%s: ERROR: Invalid data precision (%d)\n",argv[0],data_prec);
    exit(1);
  }

  /* Open the output files */
  if((fout=fopen(argv[7],"w"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[7]);
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
  if (width != 0) {
    nrows_out = (int)(nrows/width);
    ncols_out = (int)(ncols/width);
    if (nrows/width > (double)nrows_out) nrows_out++;
    if (ncols/width > (double)ncols_out) ncols_out++;
  }
  else {
    nrows_out = 1;
    ncols_out = 1;
    if (nrows > ncols) {
      width = nrows;
    }
    else {
      width = ncols;
    }
  }
  cellsize_out = cellsize * width;

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

  /* Allocate output data arrays */
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

  /* Compute statistics */
  for(i=0;i<nrows_out;i++) {
    for(j=0;j<ncols_out;j++) {
      mean = 0;
      count = 0;
      for (ii=0; ii<width; ii++) {
        old_i = i*width + ii;
        for (jj=0; jj<width; jj++) {
          old_j = j*width + jj;
          if ( old_i < nrows && old_j < ncols
               && (int_data && (int)data[old_i][old_j] != (int)nodata
                   || data[old_i][old_j] != nodata) ) {
            if (!strcasecmp(statstr,"min")) {
              if (count == 0) {
                data_out[i][j] = data[old_i][old_j];
              }
              else if (data[old_i][old_j] < data_out[i][j]) {
                data_out[i][j] = data[old_i][old_j];
              }
            }
            else if (!strcasecmp(statstr,"max")) {
              if (count == 0) {
                data_out[i][j] = data[old_i][old_j];
              }
              else if (data[old_i][old_j] > data_out[i][j]) {
                data_out[i][j] = data[old_i][old_j];
              }
            }
            else {
              mean += data[old_i][old_j];
            }
            count++;
          }
        }
      }
      sum = mean;
      if (count > 0) {
        mean /= count;
      }
      else {
        mean = nodata;
      }
      if (!strcasecmp(statstr,"sum")) {
        data_out[i][j] = sum;
      }
      else if (!strcasecmp(statstr,"mean")) {
        data_out[i][j] = mean;
      }
      else {
        var = 0;
        for (ii=0; ii<width; ii++) {
          old_i = i*width + ii;
          for (jj=0; jj<width; jj++) {
            old_j = j*width + jj;
            if ( old_i < nrows && old_j < ncols
                 && (int_data && (int)data[old_i][old_j] != (int)nodata
                     || data[old_i][old_j] != nodata) ) {
              var += (data[old_i][old_j]-mean)*(data[old_i][old_j]-mean);
            }
          }
        }
        if (count > 1) {
          var /= count-1;
        }
        else {
          var = nodata;
        }
        if (!strcasecmp(statstr,"var")) {
          data_out[i][j] = var;
        }
        else if (!strcasecmp(statstr,"std")) {
          if (count > 1) {
            data_out[i][j] = sqrt(var);
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

