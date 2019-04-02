/*
 * SUMMARY:      grid_topo_index_dinf.c - Calculate the topographic wetness index using D-infinity algorithm
 * AUTHOR:       Ted Bohn
 * ORG:          University of Washington, Department of Civil Engineering
 * DESCRIPTION:  Calculate the topographic wetness index of a DEM (ascii arc/info grid format) using D-infinity algorithm
 * DESCRIP-END.
 * Modified:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUMBANDS 10.
#define MAXSTRING 500
#define NNEIGHBORS   8
#define OUTSIDEBASIN -99
#define E_RADIUS 6371000.0 /* average radius of the earth in m*/    
#define PI 3.1415
typedef struct {
  int  Rank;
  int x;
  int y;
  float TopoIndex;
  float **dem;
} ITEM;

typedef struct {
  float Dem;
  char Mask;
  float TopoIndex;
} FINEPIX;

typedef struct {
  float Dem;
  char Mask;
} TOPOPIX;

void quick(ITEM *item,int count);
void qs(ITEM *item ,int left,int right);
double cell_area(double cent_long, double cent_lat, double size);
double  get_dist(double lat1, double long1, double lat2, double long2);

int xneighbor[NNEIGHBORS] = {
#if NNEIGHBORS == 4
 0, 1, 0, -1
#elif NNEIGHBORS == 8
  -1, 0, 1, 1, 1, 0, -1, -1
#endif
  };
 int yneighbor[NNEIGHBORS] = {
#if NNEIGHBORS == 4
  -1, 0, 1, 0
#elif NNEIGHBORS == 8
 1, 1, 1, 0, -1, -1, -1, 0
#endif
    };


/*****************************************************************************
  Function name: CalcTopoIndex()

  Purpose      : Calculate the topographic wetness index for the generation of wetness index fractional area curve.  
                 
  Returns      : float, topographic wetness index

  Modifies     : none
   
  Comments     : Based on the topographic index, ln (a/tan beta) from TOPMODEL
                 (Beven & Kirkby 1979). Calculated according to Wolock 1995.

                 (modified to use Tarboton's D-infinity algorithm)

  The surrounding grid cells are numbered in the following way

                |-----| DX

          0-----1-----2  ---
          |\    |    /|   |
          | \   |   / |   |
          |  \  |  /  |   | DY
          |   \ | /   |   |
          |    \|/    |   |
          7-----*-----3  ---
          |    /|\    |
          |   / | \   |
          |  /  |  \  |
          | /   |   \ |
          |/    |    \|
          6-----5-----4

  For the current implementation it is assumed that the resolution is the 
  same in both the X and the Y direction.  

  Source       : Beven K.J. and M.J. Kirkby, 1979, A physically based, variable 
                 contributing area model of basin hydrology, Hydrol Sci Bull 24, 
                 43-69.
        
                 Wolock, David M. and Gregory J. McCabe, Jr., 1995, Comparison
                 of single and multiple flow direction algorithms for computing
                 topographic parameters in TOPMODEL, Water Resources Research,
                 31 (5), 1315-1324.

                 Tarboton, 1997, A new method for the determination of flow
                 directions and upslope areas in grid digital elevation models,
                 Water Resources Research, 33 (2), 309-319.

*****************************************************************************/

void main(int argc ,char *argv[], FINEPIX ***FineMap)
  
{
  FILE *fw,*fa,*fb,*fdem,*fs, *ff;
  char demfile[100],wifile[100],afile[100],bfile[100],fracfile[100],soilfile[100];
  char topoindexmap[100];
  char tempstr[MAXSTRING];
  int printmap,flag;
  int i, j, k, x, y, n, lower;  /* counters */
  float dx, dy;
  int  i_corner,j_corner;
  int cols, rows;
  double xorig, yorig, size;
  float nodata;
  char tmpstr[30];
  char units[30];
  // int  nodata;
  float resolution;
  float celev;
  int count;
  float Xmax,Ymax,max,min;
  float neighbor_elev[NNEIGHBORS], neighbor_lat[NNEIGHBORS], neighbor_long[NNEIGHBORS];
  float temp_slope[NNEIGHBORS];
  float new_temp_slope[NNEIGHBORS];
  float delta_a[NNEIGHBORS];
  float length_diagonal;
  float **a;                    /* Upslope contributing area per unit contour (m2)(dem) */
  float **a_pixel;                    /* Local area of pixel (m2) */
  //float **wetnessindex;   
  double **wetnessindex;                    
  float **tanbeta;
  float **contour_length;
  float **dem;
//  float corner_lat, corner_lon;
  float lat,lon,lat1,lat2,long1,long2,cent_lat,cent_long;
  int cellnum,numcells;
  int coarsei, coarsej;
  int first;
  float VERTRES;
  //  float *OrderedCellsfine;
    ITEM *OrderedCellsfine;
    //  FINEPIX ***FineMap;
    // float *dist;
  TOPOPIX  **TopoMap;
  
  if (argc != 8)
    {
      printf("Usage : %s <demfile> <units> <vertres> <wifile> <afile> <bfile> <fracfile>\n", argv[0]);
      printf("\t\tdemfile : DEM with arcinfo header\n");
      printf("\t\tunits   : Length units for DEM lateral dimensions (\"deg\" or \"m\")\n");
      printf("\t\tvertres : Vertical resolution of DEM\n");
      printf("\t\twifile  : (output) arcinfo map of topographic wetness index values\n");
      printf("\t\tafile   : (output) arcinfo map of contributing area values\n");
      printf("\t\tbfile   : (output) arcinfo map of tanbeta values\n");
      printf("\t\tfracfile: (output) file listing the wetness index values and their fractional areas\n");
      exit(0);
    }
  strcpy(demfile,argv[1]);
  strcpy(units,argv[2]);
  VERTRES = atof(argv[3]);
  strcpy(wifile,argv[4]);
  strcpy(afile,argv[5]);
  strcpy(bfile,argv[6]);
  strcpy(fracfile,argv[7]);
		       /*------------------------------------------------------*/
		       /*	 OPEN FILES*/
		       /*------------------------------------------------------*/
		       
 if((fdem=fopen(demfile,"r"))==NULL)
   { 
     printf("cannot open/read dem file,%s\n",demfile);
     exit(1);
   }
 if((fw=fopen(wifile,"w"))==NULL)
   { printf("cannot open wifile,%s\n",wifile);
     exit(1);
   }
 if((fa=fopen(afile,"w"))==NULL)
   { printf("cannot open afile,%s\n",afile);
     exit(1);
   }
 if((fb=fopen(bfile,"w"))==NULL)
   { printf("cannot open bfile,%s\n",bfile);
     exit(1);
   }
 if((ff=fopen(fracfile,"w"))==NULL)
   { printf("cannot open fracfile,%s\n",fracfile);
     exit(1);
   }
 
 printf("Files sucessfully opened...\n\n");
		       
		       /*-------------------------------------------------*/
		       /*Scan DEM header*/
		       /*----------------------------------------------*/
		       fscanf(fdem,"%s %d",tempstr,&cols);
		       fscanf(fdem,"%s %d",tempstr,&rows);
                       fscanf(fdem,"%*s %s",tmpstr);
                       xorig = atof(tmpstr);
                       fscanf(fdem,"%*s %s",tmpstr);
                       yorig = atof(tmpstr);
                       fscanf(fdem,"%*s %s",tmpstr);
                       size = atof(tmpstr);
		       fscanf(fdem,"%s %f",tempstr,&nodata);
		       printf("%d %d %.16f %.16f %.16f \n",cols,rows,xorig,yorig,size);
		       
		       
		       /**************************************************************************/
                        /* Allocate memory to arrays*/
		       /***************************************************************************/
		       
		       if(!(dem = (float**) calloc(rows,sizeof(float*))))
			 {
			   printf("Cannot allocate memory to first record: BASIN\n");
			   exit(8); 
			 }
		       for(i=0; i<rows;i++)
			 {
			   if(!(dem[i] = (float*) calloc(cols,sizeof(float)))) {
			   printf("Cannot allocate memory to first record: BASIN\n");
			   exit(8); }
			 }

		       if(!(wetnessindex = (double**) calloc(rows,sizeof(double*))))
			 {
			   printf("Cannot allocate memory to first record: BASIN\n");
			   exit(8); 
			 }
		       for(i=0; i<rows;i++)
			 {
			   if(!(wetnessindex[i] = (double*) calloc(cols,sizeof(double)))) {
			   printf("Cannot allocate memory to first record: BASIN\n");
			   exit(8); }
			 }

		       if (!(a = (float **)calloc(rows, sizeof(float *))))
			 {
			   printf("cannot allocate memory in the basin \n");
			   exit(1);
			 }
		       for(i=0; i<rows; i++) 
			 {
			   if (!(a[i] = (float *)calloc(cols, sizeof(float))))
			     {
			       printf("cannot allocate memory in the basin  \n");
			       exit(1);
			     }
			 }
		       if (!(a_pixel = (float **)calloc(rows, sizeof(float *))))
			 {
			   printf("cannot allocate memory in the basin \n");
			   exit(1);
			 }
		       for(i=0; i<rows; i++) 
			 {
			   if (!(a_pixel[i] = (float *)calloc(cols, sizeof(float))))
			     {
			       printf("cannot allocate memory in the basin  \n");
			       exit(1);
			     }
			 }
		       if (!(tanbeta = (float **)calloc(rows, sizeof(float *))))
			 {
			   printf("cannot allocate memory in the basin \n");
			   exit(1);
			 }
		       for(i=0; i<rows; i++) 
			 {
			   if (!(tanbeta[i] = (float *)calloc(cols, sizeof(float))))
			     {
			       printf("cannot allocate memory in the basin \n");
			       exit(1);
			     }
			 }
		       if (!(contour_length = (float **)calloc(rows, sizeof(float *))))
			 {
			   printf("cannot allocate memory in the basin \n");
			   exit(1);
			 }
		       for(i=0; i<rows; i++) 
			 {
			   if (!(contour_length[i] = (float *)calloc(cols, sizeof(float))))
			     {
			       printf("cannot allocate memory in the basin \n");
			       exit(1);  
			     }
			 }

		       /*--------------------------------------------*/
		       /* READ DEM FILES                         */
		       /*--------------------------------------------*/
		       printf("Reading in files\n");
                   first = 1;
		   for(i=0; i<rows;i++) 
			 {
			   for(j=0; j<cols; j++)
			     {
			       fscanf(fdem,"%f",&dem[i][j]);
                               if (dem[i][j] != nodata) {
                                 if (first) {
                                   min = dem[i][j];
                                   max = dem[i][j];
                                   first = 0;
                                 }
                                 else {
			           if(dem[i][j]>max)
				     max=dem[i][j];
			           if(dem[i][j]<min)
				     min = dem[i][j];
                                 }
                               }
			     }				 
			   //  printf("%d \n",dem[i][j]);
			 }

		    printf("max = %f min =%f \n",max,min);

		       fclose(fdem);
		       printf("done reading DEM files\n"); 

                       numcells = rows*cols;
		       if(!(OrderedCellsfine=(ITEM*) calloc(numcells,sizeof(ITEM)))) {
			 printf("Cannot allocate memory to first record: BASIN\n");
			 exit(1); }
		       
			   /* Calculate position in dem file. */
                           i_corner = 0;
                           j_corner = 0;
			   count =0;
			   Ymax=0;
			   Xmax=0;
			   // printf("i_corner =%d j_corner=%d",i_corner,j_corner);
			 
			    for(i=i_corner; i<(i_corner+rows);i++) 
			     {
			      for(j=j_corner; j<(j_corner+cols);j++)
				{
				   if(dem[i][j]!= nodata){
				     OrderedCellsfine[count].Rank = dem[i][j];
				     OrderedCellsfine[count].y = i;
				     OrderedCellsfine[count].x = j;
				     count++;
				     /*Find Xmax and Ymax*/
				     for (y = 0; y < rows; y++) {
				       if(y>Ymax)
					 Ymax=y;
				     }
				     for (x = 0; x < cols; x++){
				       if(x>Xmax)
					 Xmax=x;
				     }
				   }
				   
				 }
			     }
			   
			   if(count==0)
			     {
			       printf("ERROR: No data in current cell. \n");
			     }
			   else
			     {
			       /* Sort OrderedCellsfine/dems into ascending order */
			       quick(OrderedCellsfine, count);
			     }
			   /* Loop through all cells in descending order of elevation */
			   for (k = count-1; k >-1; k--) 
			     { 
                             if (OrderedCellsfine[k].Rank != nodata) {
			       y = OrderedCellsfine[k].y;/* y is rows*/
			       x = OrderedCellsfine[k].x;/*x is cols*/
			       cent_long = xorig + x*size+size/2;/*centre longtude */
			       cent_lat = yorig + y*size+size/2;/*centre latitude */
			       long1=xorig+x*size;/*current long*/ 
			       lat1=yorig + y*size;/*current lat of dem*/
			       long2=xorig;/*initial lat &long of dem*/
			       lat2=yorig;
                               if (!strcasecmp(units,"deg")) {
			         a_pixel[y][x]=cell_area(cent_long,cent_lat,size);
                               }
                               else {
			         a_pixel[y][x]=size*size;
                               }
			       a[y][x]+=a_pixel[y][x];
                               if (!strcasecmp(units,"deg")) {
			         dx = get_dist( cent_lat,xorig+x*size, cent_lat,xorig+x*size+size);
			         dy = get_dist(yorig + y*size, cent_long,yorig+y*size+size, cent_long);
                               }
                               else {
                                 dx = size;
                                 dy = size;
                               }
			       // printf("dx=%f ,dy=%f , a = %f \n",dx,dy, a[y][x]);
			       length_diagonal = sqrt((pow(dx, 2)) + (pow(dy, 2)));
			       /* fill neighbor array?*/
			       for (n = 0; n < NNEIGHBORS; n++) 
				 {
				   
				   int xn = x + xneighbor[n];
				   int yn = y + yneighbor[n];
				   
				   // Initialize neighbor_elev
				   neighbor_elev[n] = (float) OUTSIDEBASIN;
				   /*check to see if xn and yn are with in dem boundries*/
				   if(xn>=0 && yn>=0 && xn<cols && yn<rows)
				     {
				       neighbor_elev[n] = ((dem[yn][xn]!=nodata) ?   dem[yn][xn] :(float) OUTSIDEBASIN);
				       neighbor_long[n]= xorig+xn*size;
				       neighbor_lat[n]= yorig + yn*size;
				     }
				 }
			       celev = dem[y][x];
			       switch (NNEIGHBORS) { 
			       case 8:
				 lower = 0;
				 for (n = 0; n < NNEIGHBORS; n++) {
				   if(neighbor_elev[n] == OUTSIDEBASIN) {
				     neighbor_elev[n] = celev;
				   }
				   
				   /* Calculating tanbeta as tanbeta * length of cell boundary between
				      the cell of interest and downsloping neighbor. */
				   if(neighbor_elev[n] < celev){
				     if(n==0 || n==2 || n==4 || n==6){
				       temp_slope[n] = (celev - neighbor_elev[n])/length_diagonal;
				       //  new_temp_slope[n] = (E_RADIUS * PI *  temp_slope[n] / 180.0)/length_diagonal;
				       // new_temp_slope[n] = get_dist(yorig + y*size,xorig + x*size, neighbor_lat[n], neighbor_long[n])/length_diagonal;
				       // contour_length[y][x] += 0.4*get_dist(yorig + y*size,xorig + x*size,xorig,yorig);
				       contour_length[y][x] += 0.2*dx+0.2*dy;
				       tanbeta[y][x] += temp_slope[n]*(0.2*dx+0.2*dy);
				       delta_a[n] = a[y][x] *temp_slope[n]*(0.2*dx+0.2*dy);
				     }
					 
				     else {
				       if(n==1||n==5){
					 temp_slope[n] = (celev - neighbor_elev[n])/dy;
					 // new_temp_slope[n] = get_dist(yorig + y*size,xorig + x*size, neighbor_lat[n], neighbor_long[n])/length_diagonal;
					 //  contour_length[y][x] += 0.6*get_dist(yorig + y*size,xorig + x*size,yorig,xorig);
					 contour_length[y][x] += 0.6*dx ;
					 tanbeta[y][x] += temp_slope[n]*0.6*dx;
					 delta_a[n] = a[y][x] *temp_slope[n]*0.6*dx;
				       }
				       else if(n==3||n==7){
					 temp_slope[n]=(celev - neighbor_elev[n])/dx;
					 contour_length[y][x] += 0.6*dy;
					 tanbeta[y][x] += temp_slope[n]*0.6*dy;
					 delta_a[n] = a[y][x] *temp_slope[n]*0.6*dy;
					 
				       }
				     }
				   }
				   else
				     lower++;
				 } /* end for (n = 0; n < NNEIGHBORS; n++) { */
				     
				 if (lower == 8){
				   /* if this is a flat area then tanbeta = sum of (0.5 * vertical resolution of elevation 
				      data)/ horizontal distance between centers of neighboring grid cells */
				   tanbeta[y][x] = ((NNEIGHBORS/2)*((0.5 * VERTRES)/length_diagonal)) +
				     ((NNEIGHBORS/2)*((0.5 * VERTRES)/dx));
				 }
				     
				 /* Distributing total upslope area to downslope neighboring cells */
				 for (n = 0; n < NNEIGHBORS; n++) {
				   if(neighbor_elev[n]<celev){  
				     switch (n) {
				     case 0:
				       a[y+1][x-1] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 1:
				       a[y+1][x] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 2:
				       a[y+1][x+1] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 3:
				       a[y][x+1] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 4:
				       a[y-1][x+1] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 5:
				       a[y-1][x] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 6:
				       a[y-1][x-1] += delta_a[n]/tanbeta[y][x];
				       break;
				     case 7:
				       a[y][x-1] += delta_a[n]/tanbeta[y][x];
				       break;
				     default:
				       printf("VICTopoIndex error has occured\n");
				       // assert(0);
				     } /* end switch (n) {*/
				   }/* end 	 if(neighbor_elev[n]<celev){ */	
				 }/*  for (n = 0; n < NNEIGHBORS; n++) { */
				 break; /*end case 8: */
    
			       case 4:
				 printf("VICTopoIndex error has occured\n");
				 //assert(0);                       /* not set up to do this */
				 break;
			       default:
				 printf("VICTopoIndex error has occured\n");
				 //	assert(0);			/* other cases don't work either */
			       } /* end  switch (NNEIGHBORS) {  */
                             } /* end if != nodata */
			     } /* end  for (k = 0; k < count-1; k++) { */

  float sum = 0;
  for(i=i_corner; i<(i_corner+rows);i++) 
    {
      for(j=j_corner; j<(j_corner+cols);j++)   
	{
  	  /*total area*/
//	  sum+=a[i][j];
	  sum+=a_pixel[i][j];
	}
    }
  // printf("cellnum %d\n",cellnum);
  float area=0;
 float prev =0;
 float Frac_area=0;
 for (k = count-1; k >-1; k--) 
   {
     y = OrderedCellsfine[k].y;
     x = OrderedCellsfine[k].x;
   if (OrderedCellsfine[k].Rank != nodata) {
     prev=area;
     // OrderedCellsfine[k].TopoIndex = log(a[y][x]/tanbeta[y][x]);
     wetnessindex[y][x] = log ((double)(a[y][x]/tanbeta[y][x]));
     Frac_area =a_pixel[y][x]/sum;
     area=Frac_area+prev;
//	  printf("wetness=%f area = %f contrib area %.4f\n",wetnessindex[y][x],a_pixel[y][x],a[y][x]);
    
   }
   else {
     wetnessindex[y][x] = nodata;
   }
   }
 
 
  /**************************************/
  /*Start of wetness index profile*/
  /**************************************/
  float  MIN_wetness = 9999;
  float  MAX_wetness = 0;
  float total_area = 0;
  float sum_wetness =0;
  float  lamda  =0;
  int gen3=0;
  first=1;
  for(i=0; i<rows;i++) 
    {
      for(j=0; j<cols;j++)   
	{
        if (wetnessindex[i][j] != nodata) {
	  /*find min and max of wetness index*/
          if (first) {
            MIN_wetness = wetnessindex[i][j];
            MAX_wetness = wetnessindex[i][j];
            first = 0;
          }
          else {
	    if(wetnessindex[i][j]<MIN_wetness)
	      MIN_wetness = wetnessindex[i][j];
	    if(wetnessindex[i][j] >MAX_wetness)
	      MAX_wetness = wetnessindex[i][j];
          }
  	  /*find the total area of DEM*/
	  total_area+=a_pixel[i][j];
	  sum_wetness += wetnessindex[i][j]*a_pixel[i][j];
	  gen3++;
	}
	}
    }
  MIN_wetness = (float)(int)MIN_wetness;
  MAX_wetness = (float)(int)MAX_wetness+1;

  lamda = sum_wetness /total_area;
  printf("%f lamda \n", lamda);	     
  int  num_bins;
  float  bin_width;
  bin_width = 0.1;
  num_bins  = (int)((MAX_wetness-MIN_wetness)/bin_width)+1;
  float  **d;  
  if(!(d = (float**) calloc(num_bins,sizeof(float*))))
    {
      printf("Cannot allocate memory to first record: BASIN\n");
      exit(8); 
   }
  for(i=0; i<num_bins;i++)
    {
      if(!(d[i] = (float*) calloc(3,sizeof(float)))) {
	printf("Cannot allocate memory to first record: BASIN\n");
	exit(8); }
    } 
  
  float sum1 = 0;
  float sum2 = 0;
  float frac_area = 0;
  float  avg_gradient=0;
  int gen1=0;
  float minWI;
  float maxWI;
  int counter;
  for (counter = 0; counter < num_bins; counter++)
    {
      minWI = (int)MIN_wetness + counter*bin_width;
      maxWI = minWI + bin_width;
      sum1 = 0;
      sum2 = 0;
      gen1=0;
      for(i=0; i<rows;i++) 
	{
	  for(j=0; j<cols; j++)
	    {
	      if(dem[i][j]!=0 && wetnessindex[i][j] != nodata)
		{
	      if(minWI<=wetnessindex[i][j] && wetnessindex[i][j]<maxWI){
		sum1+=a_pixel[i][j];
		sum2+=a_pixel[i][j]*tanbeta[i][j];
		gen1+=1;
	      }
		}
	    }
	}
      
      frac_area=sum1/total_area;
      if(gen1==0)
	{
	  avg_gradient =0;
	}
      else {
        avg_gradient= sum2/sum1;
      }
      fprintf(ff,"%d %f %f\n",counter,maxWI,frac_area);
    }

  /*************************************************************************/
    /* Create output files...creates the out in the output file in the ASCII format along with  TopoIndex generatig Wetnes index Map */
    /*************************************************************************/
    //    printf("Output.....\n");
    fprintf(fw,"ncols %d\n",cols);
    fprintf(fw,"nrows %d\n",rows);
    fprintf(fw,"xllcorner %.16f\n",xorig);
    fprintf(fw,"yllcorner %.16f\n",yorig );
    fprintf(fw,"cellsize %.16f\n",size);
    fprintf(fw,"NODATA_value %f\n",nodata);
    // k=0;
    for(i=0; i<rows;i++) 
      {
	for(j=0; j<cols; j++)
	  {
	    if(dem[i][j]!=0)
	      {
		fprintf(fw, "%lf  ",wetnessindex[i][j]);
	      }
	    else
	      fprintf(fw, "0. "); 
	    // k++;
	  }
	
	fprintf(fw,"\n");
      }
    
    fclose(fw);

    fprintf(fa,"ncols %d\n",cols);
    fprintf(fa,"nrows %d\n",rows);
    fprintf(fa,"xllcorner %.16f\n",xorig);
    fprintf(fa,"yllcorner %.16f\n",yorig );
    fprintf(fa,"cellsize %.16f\n",size);
    fprintf(fa,"NODATA_value %f\n",nodata);
    // k=0;
    for(i=0; i<rows;i++) 
      {
	for(j=0; j<cols; j++)
	  {
	    if(dem[i][j]!=0)
	      {
		fprintf(fa, "%lf  ",a[i][j]);
	      }
	    else
	      fprintf(fa, "0. "); 
	    // k++;
	  }
	
	fprintf(fa,"\n");
      }
    
    fclose(fa);

    fprintf(fb,"ncols %d\n",cols);
    fprintf(fb,"nrows %d\n",rows);
    fprintf(fb,"xllcorner %.16f\n",xorig);
    fprintf(fb,"yllcorner %.16f\n",yorig );
    fprintf(fb,"cellsize %.16f\n",size);
    fprintf(fb,"NODATA_value %f\n",nodata);
    // k=0;
    for(i=0; i<rows;i++) 
      {
	for(j=0; j<cols; j++)
	  {
	    if(dem[i][j]!=0)
	      {
		fprintf(fb, "%lf  ",tanbeta[i][j]);
	      }
	    else
	      fprintf(fb, "0. "); 
	    // k++;
	  }
	
	fprintf(fb,"\n");
      }
    
    fclose(fb);

    for(i=0; i<rows; i++) { 
      free(a[i]);
      free(a_pixel[i]);
      free(tanbeta[i]);
      free(contour_length[i]);
      free(wetnessindex[i]);
    }
   
    free(a);
    free(a_pixel);
    free(tanbeta);
    free(contour_length);
    free(OrderedCellsfine);
    free(wetnessindex);
    
}
     
   /* -------------------------------------------------------------
   QuickSort
   ------------------------------------------------------------- */

/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/

void quick(ITEM *item,int count)
  
    
{
  qs(item,0,count-1);
}
		       
void qs(ITEM *item, int left,  int right)
     
     /**********************************************************************
 this is the quick sort subroutine - it returns the values in
 an array from high to low.
     **********************************************************************/
{
  register int i,j;
  ITEM x,y;
  
  i=left;
  j=right;
  x=item[(left+right)/2];
			 
  do {
    while(item[i].Rank<x.Rank&& i<right) i++;
    while(x.Rank<item[j].Rank && j>left) j--;
    
    if (i<=j) {
      y=item[i];
      item[i]=item[j];
      item[j]=y;
      i++;
      j--;
    }
  } while (i<=j);
  
  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right);
  
}

double cell_area(double cent_long, double cent_lat, double size) {
/***********************************************************************
  cell_area.c           Keith Cherkauer             

  This program computes the area being routed using the basin, or
  sub-basin flow direction file, and the drainage fraction file.  

***********************************************************************/

  
  int    i, row, col, nrows, ncols;
  int    direc;
  double cell_lat1, cell_lng;
  double ll_lat, ll_lng;
  double dlat, tmplat;
  // double cellsize;
  double fract;
  double tmpsum, areasum;
  double routarea, tmprout;
  double NODATA;
  double ddist;

  // if(argc!=4) {
  //  fprintf(stderr,"Usage: %s <longitude> <latitude> <res>\n",argv[0]);
  //  exit(0);
  // }

  // cellsize = atof(argv[3]);

  /** Process Direction Number **/
  dlat  = size / 10.;
  ddist = get_dist(0.,0.,dlat,0.);

  areasum  = 0.;

  cell_lat1 = cent_lat;
  cell_lng = cent_long;
  tmpsum = 0.0;
  tmplat = cell_lat1-0.5*size + 0.5*dlat;
  for(i=0;i<10;i++) {
    tmpsum += get_dist(tmplat,cell_lng-0.5*size,
			     tmplat,cell_lng+0.5*size) 
      * ddist;
    tmplat += dlat;
  }
  areasum = tmpsum;

  return areasum;

  // fprintf(stdout,"%lf\n",areasum);

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

/*******************************************************************************************
  facetFlow: computes direction of steepest descent for a north-east pixel facet.
             e0 = elevation of center pixel
             e1 = elevation of east neighbor
             e2 = elevation of north-east neighbor
             d1 = east-west pixel width
             d2 = north-south pixel width
             r = direction of steepest descent, in radians
             s = slope in direction of steepest descent
*******************************************************************************************/
void facetFlow(double e0, double e1, double e2, double d1, double d2, double *r, double *s)
{

  double s1, s2;
  double diag_angle, diag_dist;

  s1 = (e0 - e1)/d1;
  s2 = (e1 - e2)/d2;

  *r = atan2(s2,s1);
  *s = sqrt(s1*s1+s2*s2);

  diag_angle = atan2(d2,d1);
  diag_dist = sqrt(d1*d1+d2*d2);

  if (*r < 0) {
    *r = 0;
    *s = s1;
  }
  else if (*r > diag_angle) {
    *r = diag_angle;
    *s = (e0 - e2)/diag_dist;
  }

}
