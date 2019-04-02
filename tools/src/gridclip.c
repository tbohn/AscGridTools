#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
  
int main (int argc, char *argv[])
/******************************************************************************
  This program takes an arcinfo ascii grid file and "clips" it
  to specified boundaries, i.e. creates a new grid file containing
  just the subset of the original grid that falls within the
  boundaries.

  Author: Ted Bohn 2007-Dec-20

  Modifications:

******************************************************************************/
{
  FILE *fin,*fout;
  double xmin,xmax;
  double ymin,ymax;
  int nrows_out,ncols_out;
  char argstr[10];
  int int_data;
  int coord_prec;
  int data_prec;
  int nrows,ncols;
  double xllcorner,yllcorner;
  double cellsize;
  int tmp_int;
  float nodata;
  double y_center,x_center;
  float tmp_data;
  int i,j;
  int ii,jj;
  float **data;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];
  
  /* Usage */
  if(argc!=10) {
    fprintf(stdout,"Usage: %s <in_grid> <type> <xmin> <xmax> <ymin> <ymax> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid>    Input grid file name\n");
    fprintf(stdout,"  <type>       Data type (\"int\" or \"float\")\n");
    fprintf(stdout,"  <xmin>       Minimum x coordinate (western boundary)\n");
    fprintf(stdout,"  <xmax>       Maximum x coordinate (eastern boundary)\n");
    fprintf(stdout,"  <ymin>       Minimum y coordinate (southern boundary)\n");
    fprintf(stdout,"  <ymax>       Minimum y coordinate (northern boundary)\n");
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

  /* Read the data type */
  strcpy(argstr,argv[2]);
  int_data = 0;
  if(!strcasecmp(argstr,"int"))   {
    int_data = 1;
  }

  /* Read min and max x and y values */
  xmin = atof(argv[3]);
  xmax = atof(argv[4]);
  ymin = atof(argv[5]);
  ymax = atof(argv[6]);

  /* Read coordinate and data precisions */
  coord_prec = atoi(argv[7]);
  data_prec = atoi(argv[8]);

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

  /* Compute dimensions of output data array */
  ncols_out = (int)((xmax - xmin)/cellsize);
  if ((xmax - xmin)/cellsize > (double)ncols_out) ncols_out++;
  nrows_out = (int)((ymax - ymin)/cellsize);
  if ((ymax - ymin)/cellsize > (double)nrows_out) nrows_out++;
  xmax = xmin + ncols_out*cellsize;
  ymax = ymin + nrows_out*cellsize;

  /* Allocate data array */
  if ( (data = (float**)calloc(nrows_out,sizeof(float*))) == NULL ) {
    fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data array\n",argv[0]);
    exit(1);
  }
  for (i=0; i<nrows_out; i++) {
    if ( (data[i] = (float*)calloc(ncols_out,sizeof(float))) == NULL ) {
      fprintf(stderr,"%s: ERROR: cannot allocate sufficient memory for data array\n",argv[0]);
      exit(1);
    }
  }

  /* Initialize data array */
  for(i=0;i<nrows_out;i++) {
    for(j=0;j<ncols_out;j++) {
      data[i][j] = nodata;
    }
  }

  /* Initialize y_center and x_center */
  x_center = xllcorner + cellsize*0.5;
  y_center = yllcorner + cellsize*(nrows-0.5);

  /* Read in the data */
  for(i=0;i<nrows;i++) {
    for(j=0;j<ncols;j++) {
      if (int_data) {
        fscanf(fin,"%d",&tmp_int);
      }
      else {
        fscanf(fin,"%f",&tmp_data);
      }
      if (y_center >= ymin && y_center < ymax) {
        if (x_center >= xmin && x_center < xmax) {
          ii = nrows_out - 1 - (int)((y_center-ymin)/cellsize);
          jj = (int)((x_center-xmin)/cellsize);
          if (int_data) {
            data[ii][jj] = (float)tmp_int;
          }
          else {
            data[ii][jj] = tmp_data;
          }
        }
      }
      x_center += cellsize;
    }
    y_center -= cellsize;
    x_center = xllcorner + cellsize*0.5;
  }
  fclose(fin); 


  /* Write the output header */
  fprintf(fout,"ncols %d\n",ncols_out);
  fprintf(fout,"nrows %d\n",nrows_out);
  sprintf(precstr,"%%");
  sprintf(precstr2,".%df",coord_prec);
  strcat(precstr,precstr2);
  sprintf(fmtstr,"xllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,xmin);
  fprintf(fout,"\n");
  sprintf(fmtstr,"yllcorner ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,ymin);
  fprintf(fout,"\n");
  sprintf(fmtstr,"cellsize ");
  strcat(fmtstr,precstr);
  fprintf(fout,fmtstr,cellsize);
  fprintf(fout,"\n");
  if (int_data) {
    tmp_int = (int)nodata;
    fprintf(fout,"nodata %d\n",tmp_int);
  }
  else {
    sprintf(precstr,"%%");
    sprintf(precstr2,".%df",data_prec);
    strcat(precstr,precstr2);
    sprintf(fmtstr,"nodata ");
    strcat(fmtstr,precstr);
    fprintf(fout,fmtstr,nodata);
    fprintf(fout,"\n");
  }

  /* Write the output data */
  if (int_data) {
    sprintf(fmtstr,"%%d");
  }
  else {
    sprintf(precstr,"%%");
    sprintf(precstr2,".%df",data_prec);
    strcat(precstr,precstr2);
    strcpy(fmtstr,precstr);
  }
  for (ii=0; ii<nrows_out; ii++) {
    for (jj=0; jj<ncols_out; jj++) {
      if (int_data) {
        fprintf(fout,fmtstr,(int)data[ii][jj]);
      }
      else {
        fprintf(fout,fmtstr,data[ii][jj]);
      }
      if (jj < ncols_out-1) {
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

