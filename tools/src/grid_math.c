#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char *argv[])
/******************************************************************************
  Author: Ted Bohn 2007-Dec-20

  Description:
  This program takes two arcinfo ascii grid files and performs a mathematical
  operation on them (see usage for the list of choices).

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
  double nodata1, nodata2, nodata_out;
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
  if(argc!=12) {
    fprintf(stdout,"Usage: %s <in_grid_1> <type1> <in_grid_2> <type2> <operation> <nodata_out> <type_out> <union> <coord_prec> <data_prec> <out_grid>\n",argv[0]);
    fprintf(stdout,"  <in_grid_1>  First input grid file name, or a numerical value (when type1 = const)\n");
    fprintf(stdout,"  <type1>      Data type of first file (\"int\" or \"float\" or \"const\")\n");
    fprintf(stdout,"  <in_grid_2>  Second input grid file name, or a numerical value (when type2 = const)\n");
    fprintf(stdout,"  <type2>      Data type of second file (\"int\" or \"float\" or \"const\")\n");
    fprintf(stdout,"  <operation>  Mathematical operation (\"+\",\"-\",\"*\",\"/\",\"min\",\"max\",\"avg\",\"gt\",\"ge\",\"eq\",\"le\",\"lt\")\n");
    fprintf(stdout,"  <nodata_out> Nodata value for the output file\n");
    fprintf(stdout,"  <type_out>   Data type for the output file (\"int\" or \"float\")\n");
    fprintf(stdout,"  <union>      What to do for cells that are nodata in one of the files; 0 = set to nodata_out; 1 = take value from other file\n");
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

  /* Operation */
  strcpy(operation,argv[5]);

  /* Nodata_out */
  nodata_out = atof(argv[6]);

  /* Data type - output file */
  strcpy(argstr,argv[7]);
  int_data_out = 0;
  if(!strcasecmp(argstr,"int")) {
    int_data_out = 1;
  }

  /* Union switch */
  union_switch = atoi(argv[8]);

  /* Coordinate and data precision */
  coord_prec = atoi(argv[9]);
  data_prec = atoi(argv[10]);

  /* Open the output file */
  if((fout=fopen(argv[11],"w"))==NULL) {
    fprintf(stderr,"%s: ERROR: Unable to open %s\n",argv[0],argv[11]);
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
      nodata1 = (double)tmp_int1;
    }
    else {
      fscanf(fin1,"%*s %s",tmpstr2);
      nodata1 = atof(tmpstr2);
    }
  }
  else {
    nodata1 = (double)(tmp_float1+1);
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
      nodata2 = (double)tmp_int1;
    }
    else {
      fscanf(fin2,"%*s %s",tmpstr2);
      nodata2 = atof(tmpstr2);
    }
  }
  else {
    nodata2 = (double)(tmp_float2+1);
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
      if (tmp_float1 != nodata1 && tmp_float2 != nodata2) {
        if (!strcmp(operation,"+")) {
          data[i][j] = tmp_float1 + tmp_float2;
        }
        else if (!strcmp(operation,"-")) {
          data[i][j] = tmp_float1 - tmp_float2;
        }
        else if (!strcmp(operation,"*")) {
          data[i][j] = tmp_float1 * tmp_float2;
        }
        else if (!strcmp(operation,"/")) {
          if (tmp_float2 > 0) {
            data[i][j] = tmp_float1 / tmp_float2;
          }
          else {
            data[i][j] = nodata_out;
          }
        }
        else if (!strcmp(operation,"min")) {
          data[i][j] = (tmp_float1 > tmp_float2) ? tmp_float2 : tmp_float1;
        }
        else if (!strcmp(operation,"max")) {
          data[i][j] = (tmp_float1 > tmp_float2) ? tmp_float1 : tmp_float2;
        }
        else if (!strcmp(operation,"avg")) {
          data[i][j] = 0.5 * (tmp_float1 + tmp_float2);
        }
        else if (!strcmp(operation,"gt")) {
          data[i][j] = (tmp_float1 > tmp_float2) ? 1 : 0;
        }
        else if (!strcmp(operation,"ge")) {
          data[i][j] = (tmp_float1 >= tmp_float2) ? 1 : 0;
        }
        else if (!strcmp(operation,"eq")) {
          data[i][j] = (tmp_float1 == tmp_float2) ? 1 : 0;
        }
        else if (!strcmp(operation,"le")) {
          data[i][j] = (tmp_float1 <= tmp_float2) ? 1 : 0;
        }
        else if (!strcmp(operation,"lt")) {
          data[i][j] = (tmp_float1 < tmp_float2) ? 1 : 0;
        }
        else {
          data[i][j] = nodata_out;
        }
      }
      else {
        if (union_switch) {
          if (tmp_float1 == nodata1 && tmp_float2 != nodata2) {
            data[i][j] = tmp_float2;
          }
          else if (tmp_float1 != nodata1 && tmp_float2 == nodata2) {
            data[i][j] = tmp_float1;
          }
          else {
            data[i][j] = nodata_out;
          }
        }
        else {
          data[i][j] = nodata_out;
        }
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
  sprintf(precstr2,".%dg",coord_prec);
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

