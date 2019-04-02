#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2007-Dec-20

  Description:

  This program "overlays" one arcinfo ascii grid file on top of another, i.e.
  the data in the output grid equal the data in the first input grid, with gaps
  in the data from the first input grid filled by data from the second input grid.

  In other words, the grid cells (pixels) in the output grid are assigned as follows:

    if (data_1 != NODATA_1)
      data_out = data_1
    else if (data_2 != NODATA_2)
      data_out = data_2
    else
      data_out = NODATA_OUT


  NOTE: both grids must be of same dimensions and resolution.

  Modifications:

******************************************************************************/
{
  FILE *fin1,*fin2,*fout;
  int nrows,ncols;
  char argstr[10];
  int int_data1, int_data2, int_data_out;
  int coord_prec;
  int data_prec;
  double xllcorner,yllcorner;
  double cellsize;
  int tmp_int1,tmp_int2;
  float nodata1, nodata2, nodata_out;
  float tmp_float1,tmp_float2;
  int i,j;
  float **data;
  char precstr[10];
  char tmpstr[10];
  char precstr2[10];
  char fmtstr[30];
  char tmpstr2[30];
  char operation[10];
  int union_switch;
  int const1,const2;
  
  /* Usage */
  if(argc!=10) {
    fprintf(stdout,"Usage: %s <in_grid_1> <type1> <in_grid_2> <type2> <nodata_out> <type_out> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid_1>  First input grid file name, or a numerical value (when type1 = const)\n");
    fprintf(stdout,"  <type1>      Data type of first file (\"int\" or \"float\" or \"const\")\n");
    fprintf(stdout,"  <in_grid_2>  Second input grid file name, or a numerical value (when type2 = const)\n");
    fprintf(stdout,"  <type2>      Data type of second file (\"int\" or \"float\" or \"const\")\n");
    fprintf(stdout,"  <nodata_out> Nodata value for the output file\n");
    fprintf(stdout,"  <type_out>   Data type for the output file (\"int\" or \"float\")\n");
    fprintf(stdout,"  <coord_prec> Precision of coordinates, i.e. number of decimal places\n");
    fprintf(stdout,"  <data_prec>  Precision of data, i.e. number of decimal places (ignored for \"int\" data)\n");
    fprintf(stdout,"  <out_grid>   Output grid file name\n");
    exit(0);
  }

  /* Data type - file 1 */
  strcpy(argstr,argv[2]);
  int_data1 = 0;
  const1 = 0;
  if(!strcasecmp(argstr,"int")) {
    int_data1 = 1;
  }
  else if (!strcasecmp(argstr,"const")) {
    const1 = 1;
  }

  /* Open input file 1 */
  if (const1) {
    tmp_float1 = atof(argv[1]);
  }
  else {
    if((fin1=fopen(argv[1],"r"))==NULL) {
      fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[1]);
      exit(1);
    }
  }

  /* Data type - file 2 */
  strcpy(argstr,argv[4]);
  int_data2 = 0;
  const2 = 0;
  if(!strcasecmp(argstr,"int")) {
    int_data2 = 1;
  }
  else if (!strcasecmp(argstr,"const")) {
    const2 = 1;
  }

  /* Open input file 2 */
  if (const2) {
    tmp_float2 = atof(argv[3]);
  }
  else {
    if((fin2=fopen(argv[3],"r"))==NULL) {
      fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[3]);
      exit(1);
    }
  }

  if (const1 && const2) {
    fprintf(stderr,"%s: ERROR: both inputs cannot be constants\n",argv[0]);
    exit(1);
  }

  /* Nodata_out */
  nodata_out = atof(argv[5]);

  /* Data type - output file */
  strcpy(argstr,argv[6]);
  int_data_out = 0;
  if(!strcasecmp(argstr,"int")) {
    int_data_out = 1;
  }

  /* Coordinate and data precision */
  coord_prec = atoi(argv[7]);
  data_prec = atoi(argv[8]);

  /* Open the output file */
  if((fout=fopen(argv[9],"w"))==NULL) {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[9]);
    exit(1);
  }

  if (!const1) {
    /* Read in the header - file 1 */
    fscanf(fin1,"%*s %d",&ncols);
    fscanf(fin1,"%*s %d",&nrows);
    fscanf(fin1,"%*s %s",tmpstr2);
    xllcorner = atof(tmpstr2);
    fscanf(fin1,"%*s %s",tmpstr2);
    yllcorner = atof(tmpstr2);
    fscanf(fin1,"%*s %s",tmpstr2);
    cellsize = atof(tmpstr2);
    if (int_data1) {
      fscanf(fin1,"%*s %d",&tmp_int1);
      nodata1 = (float)tmp_int1;
    }
    else {
      fscanf(fin1,"%*s %f",&nodata1);
    }
  }
  else {
    nodata1 = tmp_float1+1;
  }

  if (!const2) {
    /* Read in the header - file 2 */
    fscanf(fin2,"%*s %d",&ncols);
    fscanf(fin2,"%*s %d",&nrows);
    fscanf(fin2,"%*s %s",tmpstr2);
    xllcorner = atof(tmpstr2);
    fscanf(fin2,"%*s %s",tmpstr2);
    yllcorner = atof(tmpstr2);
    fscanf(fin2,"%*s %s",tmpstr2);
    cellsize = atof(tmpstr2);
    if (int_data2) {
      fscanf(fin2,"%*s %d",&tmp_int1);
      nodata2 = (float)tmp_int1;
    }
    else {
      fscanf(fin2,"%*s %f",&nodata2);
    }
  }
  else {
    nodata2 = tmp_float2+1;
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

  /* Read in the data and perform operation on it */
  for(i=0;i<nrows;i++) {
    for(j=0;j<ncols;j++) {
      if (!const1) {
        if (int_data1) {
          fscanf(fin1,"%d",&tmp_int1);
          tmp_float1 = (float)tmp_int1;
        }
        else {
          fscanf(fin1,"%f",&tmp_float1);
        }
      }
      if (!const2) {
        if (int_data2) {
          fscanf(fin2,"%d",&tmp_int2);
          tmp_float2 = (float)tmp_int2;
        }
        else {
          fscanf(fin2,"%f",&tmp_float2);
        }
      }
      if (tmp_float1 != nodata1) {
        data[i][j] = tmp_float1;
      }
      else if (tmp_float2 != nodata2) {
        data[i][j] = tmp_float2;
      }
      else {
        data[i][j] = nodata_out;
      }
    }
  }
  if (!const1) {
    fclose(fin1); 
  }
  if (!const2) {
    fclose(fin2); 
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
  if (int_data_out) {
    tmp_int1 = (int)nodata_out;
    fprintf(fout,"nodata %d\n",tmp_int1);
  }
  else {
    sprintf(precstr,"%%");
    sprintf(precstr2,".%df",data_prec);
    strcat(precstr,precstr2);
    sprintf(fmtstr,"nodata ");
    strcat(fmtstr,precstr);
    fprintf(fout,fmtstr,nodata_out);
    fprintf(fout,"\n");
  }

  /* Write the output data */
  if (int_data_out) {
    sprintf(fmtstr,"%%d");
  }
  else {
    sprintf(precstr,"%%");
    sprintf(precstr2,".%df",data_prec);
    strcat(precstr,precstr2);
    strcpy(fmtstr,precstr);
  }
  for (i=0; i<nrows; i++) {
    for (j=0; j<ncols; j++) {
      if (int_data_out) {
        fprintf(fout,fmtstr,(int)data[i][j]);
      }
      else {
        fprintf(fout,fmtstr,data[i][j]);
      }
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

