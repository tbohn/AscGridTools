# Descriptions and Usage of Grid Tools

 - gridclip

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

grid_agg.c
grid_agg_class.c
grid_apply_mask.c
gridclip.c
grid_latlon2utm.c
grid_make_mask_thresh.c
grid_math.c
grid_overlay.c
grid_smooth.c
grid_stats.c
grid_subsample.c
grid_topo_index.c
grid_topo_index_dinf.c
grid_utm2latlon.c
