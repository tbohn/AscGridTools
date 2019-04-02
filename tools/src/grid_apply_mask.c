#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2008-Jun-20

  Description:
  This program applies a mask (in Arc/Info ascii grid format) to a grid (also
  in Arc/Info ascii grid format).  Pixels in the input grid that correspond to
  mask pixels having "NODATA" values are set to the input grid's "NODATA" value.

  NOTE: This program currently only handles masks of type integer.  Also, mask
  has to have the same grid resolution and boundaries as the input data grid.

  Modifications:

******************************************************************************/
{
  FILE *fin,*fmask,*fout;
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
  float **data;
  int tmp_int;
  int count;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];
  int nodata_mask;
  int **mask;

  /* Usage */
  if(argc!=7) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <mask_file> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <mask_file>  Mask file\n");
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

  /* Open mask file */
  if((fmask=fopen(argv[3],"r"))==NULL)   {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[3]);
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

  /* Read in the header - data file */
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

  /* Read in the header - mask file */
  fscanf(fmask,"%*s %d",&ncols);
  fscanf(fmask,"%*s %d",&nrows);
  fscanf(fmask,"%*s %s",tmpstr2);
  xllcorner = atof(tmpstr2);
  fscanf(fmask,"%*s %s",tmpstr2);
  yllcorner = atof(tmpstr2);
  fscanf(fmask,"%*s %s",tmpstr2);
  cellsize = atof(tmpstr2);
  fscanf(fmask,"%*s %d",&nodata_mask);

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

  /* Allocate mask array */
  if ( (mask = (int**)calloc(nrows,sizeof(int*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for mask array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows; i++) {
    if ( (mask[i] = (int*)calloc(ncols,sizeof(int))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for mask array\n",argv[0]);
      exit(1);
    }
  }

  /* Read in the data & mask, and apply mask to data */
  for(i=0;i<nrows;i++) {
    for(j=0;j<ncols;j++) {
      fscanf(fmask,"%d",&mask[i][j]);
      if (int_data) {
        fscanf(fin,"%d",&tmp_int);
        data[i][j] = (float)tmp_int;
      }
      else {
        fscanf(fin,"%f",&data[i][j]);
      }
      if (mask[i][j] == nodata_mask && data[i][j] != nodata) {
        data[i][j] = nodata;
      }
    }
  }
  fclose(fin); 
  fclose(fmask); 

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
  fprintf(fout,fmtstr,nodata);
  fprintf(fout,"\n");

  /* Write the output data */
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",data_prec);
  strcat(precstr,precstr2);
  strcpy(fmtstr,precstr);
  for (i=0; i<nrows; i++) {
    for (j=0; j<ncols; j++) {
      fprintf(fout,fmtstr,data[i][j]);
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

