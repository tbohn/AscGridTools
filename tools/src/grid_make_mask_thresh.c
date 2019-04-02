#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:
  This program makes a mask from a grid file.  The user may specify a threshold
  value and a condition, such as "lt" for "less than", "ge" for "greater than or
  equal to", "eq" for "equal to", etc.  All pixels satisfying the condition in
  relation to the threshold value are set to 1 in the output mask.  All other pixels
  are set to 0 in the output mask.  The output mask's "NODATA" value is set to 0.

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
  int tmp_int;
  int count;
  int old_i,old_j;
  int res_ratio;
  float nrows_out_tmp, ncols_out_tmp;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];
  char condition[3];
  float threshold;

  /* Usage */
  if(argc!=7) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <condition> <threshold> <coord_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Input data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <condition>  lt, le, eq, ge, gt\n");
    fprintf(stdout,"  <threshold>  Value to compare grid values to; if condition is satisfied, returns 1, else 0\n");
    fprintf(stdout,"  <coord_prec> Precision of coordinates, i.e. number of decimal places\n");
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

  /* Read condition */
  strcpy(condition,argv[3]);

  /* Read threshold */
  threshold = atof(argv[4]);

  /* Read coord precision */
  coord_prec = atoi(argv[5]);

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
      if (data[i][j] != nodata) {
        if (!strcmp(condition,"lt") && data[i][j] < threshold) {
          data[i][j] = 1;
        }
        else if (!strcmp(condition,"le") && data[i][j] <= threshold) {
          data[i][j] = 1;
        }
        else if (!strcmp(condition,"eq") && data[i][j] == threshold) {
          data[i][j] = 1;
        }
        else if (!strcmp(condition,"ge") && data[i][j] >= threshold) {
          data[i][j] = 1;
        }
        else if (!strcmp(condition,"gt") && data[i][j] > threshold) {
          data[i][j] = 1;
        }
        else {
          data[i][j] = 0;
        }
      }
      else {
        data[i][j] = 0;
      }
    }
  }
  fclose(fin); 

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
  fprintf(fout,"nodata 0\n");

  /* Write the output data */
  for (i=0; i<nrows; i++) {
    for (j=0; j<ncols; j++) {
      fprintf(fout,"%d",(int)data[i][j]);
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

