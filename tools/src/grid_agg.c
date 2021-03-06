#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:
  This program aggregates a grid to a coarser resolution.

  NOTE: This program currently only aggregates by averaging, and only outputs
  float/double data.

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
  double xllcorner,yllcorner,yllcorner_out;
  double cellsize;
  float nodata;
  double cellsize_out;
  float **data,**data_out;
  int tmp_int;
  int count;
  int old_i,old_j;
  int res_ratio;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];

  /* Usage */
  if(argc!=7) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <res_ratio> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <res_ratio>  Ratio of output/input resolution (or ratio of output/input cellsize) - must be integer\n");
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

  /* Read resolution ratio */
  res_ratio = atoi(argv[3]);
  if (res_ratio <= 0) {
    fprintf(stderr,"%s: ERROR: ratio of output to input resolution (%s) is not a positive integer\n",argv[0],argv[3]);
    exit(1);
  }

  /* Read coord and data precision */
  coord_prec = atoi(argv[4]);
  data_prec = atoi(argv[5]);
  if (data_prec < 0) {
    fprintf(stderr,"%s: ERROR: Invalid data precision (%d)\n",argv[0],data_prec);
    exit(1);
  }

  /* Open the output file */
  if((fout=fopen(argv[6],"w"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[6]);
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
  ncols_out = (int)(ncols/res_ratio);
  if (ncols/res_ratio > (double)ncols_out) ncols_out++;
  nrows_out = (int)(nrows/res_ratio);
  if (nrows/res_ratio > (double)nrows_out) nrows_out++;
  cellsize_out = cellsize * res_ratio;
//  yllcorner_out = yllcorner + nrows*cellsize - nrows_out*cellsize_out;
  yllcorner_out = yllcorner;

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

  /* Aggregate to new resolution */
  for(i=nrows_out-1;i>=0;i--) {
    for(j=0;j<ncols_out;j++) {
      data_out[i][j] = 0;
      count = 0;
      for (ii=0; ii<res_ratio; ii++) {
        old_i = nrows-1 - ((nrows_out-1-i)*res_ratio + ii);
        for (jj=0; jj<res_ratio; jj++) {
          old_j = j*res_ratio + jj;
          if ( old_i >= 0 && old_i < nrows && old_j >= 0 && old_j < ncols
               && (int_data && data[old_i][old_j] != (int)nodata
                   || data[old_i][old_j] != nodata) ) {
            data_out[i][j] += data[old_i][old_j];
            count++;
          }
        }
      }
      if (count > 0) {
        data_out[i][j] /= count;
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
  fprintf(fout,fmtstr,xllcorner);
  fprintf(fout,"\n");
  sprintf(fmtstr,"yllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,yllcorner_out);
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

