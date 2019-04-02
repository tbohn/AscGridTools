# Descriptions and Usage of Grid Tools  
  
 - `grid_agg` - aggregates a grid to a coarser resolution via spatial averaging.  
  
   Usage: grid_agg <in_grid> <type> <res_ratio> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <res_ratio>  Ratio of output/input resolution (or ratio of output/input cellsize) - must be integer  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_agg_class` - aggregates a grid of land cover class pixels to a coarser resolution. The output is a set of several grid files, one per each class present in the input grid. Each output grid file contains a map of the fractional area covered by that class in each output grid cell.  
  
   Usage: grid_agg_class <in_grid> <res_ratio> <classID> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <res_ratio>  Ratio of output/input resolution (or ratio of output/input cellsize) - must be integer  
     <classID>    ID number of the class to compute fractions of  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_apply_mask` - applies the given mask to the input grid; all pixels outside the "valid" regions of the mask will be set to nodata values.  
  
   Usage: grid_apply_mask <in_grid> <type> <mask_file> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <mask_file>  Mask file  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `gridclip` - clips an input grid to a smaller (or larger) box bounded by the given min/max x/y coordinates. Where the output box extends beyond the limits of the input grid, the grid cells are assigned nodata values.  
  
   Usage: gridclip <in_grid> <type> <xmin> <xmax> <ymin> <ymax> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <xmin>       Minimum x coordinate (western boundary)  
     <xmax>       Maximum x coordinate (eastern boundary)  
     <ymin>       Minimum y coordinate (southern boundary)  
     <ymax>       Minimum y coordinate (northern boundary)  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_latlon2utm` - reprojects a grid from geographic projection to Universal Transverse Mercator (UTM) projection.  
  
   Usage: grid_latlon2utm <in_grid> <type> <zone> <resolution> <minx> <maxx> <miny> <maxy> <radius> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <zone>       UTM longitudinal zone  
     <resolution> Output resolution (cellsize, in m)  
     <minx>       Western boundary of output grid (m)  
     <maxx>       Eastern boundary of output grid (m)  
     <miny>       Southern boundary of output grid (m)  
     <maxy>       Northern boundary of output grid (m)  
     <radius>     Radius (i.e. sigma, or standard deviation) of gaussian kernel (pixels)  
                  NOTE: the gaussian kernel will act as a low-pass filter on the data,  
                  with cutoff wavelength equal to 2*PI*radius (pixels).  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_make_mask_thresh` - xxx  
  
   Usage: grid_make_mask_thresh <in_grid> <type> <condition> <threshold> <coord_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Input data type ("int" or "float")  
     <condition>  lt, le, eq, ge, gt  
     <threshold>  Value to compare grid values to; if condition is satisfied, returns 1, else 0  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <out_grid>   Output grid file name  
  
 - `grid_math` - xxx  
  
   Usage: grid_math <in_grid_1> <type1> <in_grid_2> <type2> <operation> <nodata_out> <type_out> <union> <coord_prec> <data_prec> <out_grid>  
     <in_grid_1>  First input grid file name, or a numerical value (when type1 = const)  
     <type1>      Data type of first file ("int" or "float" or "const")  
     <in_grid_2>  Second input grid file name, or a numerical value (when type2 = const)  
     <type2>      Data type of second file ("int" or "float" or "const")  
     <operation>  Mathematical operation ("+","-","*","/","min","max","avg","gt","ge","eq","le","lt")  
     <nodata_out> Nodata value for the output file  
     <type_out>   Data type for the output file ("int" or "float")  
     <union>      What to do for cells that are nodata in one of the files; 0 = set to nodata_out; 1 = take value from other file  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_overlay` - xxx  
  
   Usage: grid_overlay <in_grid_1> <type1> <in_grid_2> <type2> <nodata_out> <type_out> <coord_prec> <data_prec> <out_grid>  
     <in_grid_1>  First input grid file name, or a numerical value (when type1 = const)  
     <type1>      Data type of first file ("int" or "float" or "const")  
     <in_grid_2>  Second input grid file name, or a numerical value (when type2 = const)  
     <type2>      Data type of second file ("int" or "float" or "const")  
     <nodata_out> Nodata value for the output file  
     <type_out>   Data type for the output file ("int" or "float")  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_smooth` - xxx  
  
   Usage: grid_smooth <in_grid> <type> <method> <length> <geog> <trunc> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <method>     Smoothing method ("mean","gauss")  
     <length>     For "mean", length = width (in pixels) of smoothing window  
                  For "gauss", length = cut-off wavelength (in pixels);  
                    in this case, smoothing window width will be set to:  
                      2*((int)(3*sigma))+1  
                    where  
                      sigma = length/(2*pi)  
                            = the radius of the gaussian kernel's inflection point  
                    so that the window contains 3*sigma on each side of the central pixel  
     <geog>       1 = rows and columns are latitude and longitude, respectively (unequal in size), and the geographical aspect ratio will be taken into account; 0 = rows and columns are equal in size  
     <trunc>      Number of extreme values in the window to ignore on each side of the distribution (generally set this to 0 unless you want a truncated mean)  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_stats` - xxx  
  
   Usage: grid_stats <in_grid> <type> <stat> <width> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <stat>       Statistic to compute ("mean", "var", "std", "min", "max", "sum")  
     <width>      Width of analysis window; this is the resolution at which statistics will be output; a value of 0 == window encompasses entire grid  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file prefix; ".mean.asc", ".std.asc", etc. will be appended to this prefix to build the various output filenames  
  
 - `grid_subsample` - xxx  
  
   Usage: grid_subsample <in_grid> <type> <cellsize_out> <method> <length> <type_out> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <cellsize_out> Output cell width, in input units  
     <method>     Sub-sampling method; "nn" = nearest neighbor; "mean" = moving average; "bilin" = bi-linear interpolation, "gauss" = gaussian smoothing.  
     <length>     Characteristic length (in output pixels) of averaging method:  
                    For method="nn" or "bilin", length is ignored (can be 0)  
                    For method="mean", length=width of averaging window  
                    For method="gauss", length=low-pass filter cut-off wavelength;  
                      in this case, sigma, the radius of the gaussian inflection point, is:  
                        sigma = length/(2*PI)  
                      and the smoothing window width will be set to:  
                        2*((int)(3*sigma))+1  
                      so that the window contains 3*sigma on each side of the central pixel  
     <type_out>   Output data type ("int" or "float")  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
 - `grid_topo_index` - xxx  
  
   Usage: grid_topo_index <demfile> <units> <vertres> <afile> <bfile> <wifile> <fracfile>  
     <demfile>  DEM with arcinfo header  
     <units>    Length units for DEM lateral dimensions ("deg" or "m")  
     <vertres>  Vertical resolution of DEM  
     <afile>    (output) arcinfo map of pixel upslope contributing areas  
     <bfile>    (output) arcinfo map of pixel slopes  
     <wifile>   (output) arcinfo map of topographic wetness index values (= ln(a/b))  
     <fracfile> (output) file listing the wetness index values and their fractional areas  
  
 - `grid_topo_index_dinf` - xxx  
  
   Usage: grid_topo_index_dinf <demfile> <units> <vertres> <wifile> <afile> <bfile> <fracfile>  
     <demfile>  DEM with arcinfo header  
     <units>    Length units for DEM lateral dimensions ("deg" or "m")  
     <vertres>  Vertical resolution of DEM  
     <wifile>   (output) arcinfo map of topographic wetness index values  
     <afile>    (output) arcinfo map of contributing area values  
     <bfile>    (output) arcinfo map of tanbeta values  
     <fracfile> (output) file listing the wetness index values and their fractional areas  
  
 - `grid_utm2latlon` - xxx  
  
   Usage: grid_utm2latlon <in_grid> <type> <zone> <resolution> <minlon> <maxlon> <minlat> <maxlat> <coord_prec> <data_prec> <out_grid>  
     <in_grid>    Input grid file name  
     <type>       Data type ("int" or "float")  
     <zone>       UTM longitudinal zone  
     <resolution> Output resolution (cellsize, in degrees)  
     <minlon>     Western boundary of output grid (degrees)  
     <maxlon>     Eastern boundary of output grid (degrees)  
     <minlat>     Southern boundary of output grid (degrees)  
     <maxlat>     Northern boundary of output grid (degrees)  
     <coord_prec> Precision of coordinates, i.e. number of decimal places  
     <data_prec>  Precision of data, i.e. number of decimal places (ignored for "int" data)  
     <out_grid>   Output grid file name  
  
