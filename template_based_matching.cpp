//ToDo:
//  XXX Check for cancel more frequently?
//  Somehow let the user load multiple templates at once.
//  Let the user save the splatted-plane data set so that they can compare
//     with the original data set.
//  XXX Problem: If the volume scale is greater than one, we need to rasterize the
//    values into the plane because each volume element will fill more than one, so
//    we end up skipping points when stepping in the template.  We could also do some
//    sort of post-processing blur to fix this, maybe.

//------------------------------------------------------------------------
// This program will read in AFM files (or image files) and perform
// template-based matching on them.  See the file template_afm_estimation.doc
// in the source directory for a full description of the algorithm.
//
// The development of this program was supported by the NIH National Research
// Resource on "Computer-Integrated Systems for Microscopy and Manipulation"
// at UNC; NIH grant number P41 EB002025, and by the University of North
// Carolina at Chapel Hill.
//
// Written by: Russell M. Taylor II
// Created: October 23, 2005
//------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <tcl.h>
#include <tk.h>
#include "Tcl_Linkvar.h"
#ifdef	_WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glut.h>

#include <Topo.h>
#include <BCGrid.h>
#include <BCPlane.h>

// This pragma tells the compiler not to tell us about truncated debugging info
// due to name expansion within the string, list, and vector classes.
#pragma warning( disable : 4786 )
#include <string>
#include <algorithm>
#include <vector>
#include <list>
using namespace std;

//------------------------------------------------------------------------
// Horrible Windows compatibility stuff
#ifdef _WIN32
const char *READFLAGS = "rb";
const char *WRITEFLAGS = "cwb";
#else
const char *READFLAGS = "r";
const char *WRITEFLAGS = "cw";
#endif

#ifndef	_WIN32
#include <sys/time.h>
static void Sleep(const long milliSeconds)
{
  struct timeval delay;
  delay.tv_sec = milliSeconds / 1000;
  delay.tv_usec = (milliSeconds - 1000*delay.tv_sec) * 1000;
  select(0, NULL, NULL, NULL, &delay);
}
#endif

#ifndef	M_PI
const double M_PI = 3.1415926535898;
#endif

//------------------------------------------------------------------------
// Constant and type definitions
const char *Version_string = "01.05";
typedef double	Color[4]; // RGBA between zero and one

//----------------------------------------------------------------------------
// Class to deal with reading from and writing to image buffers.

class image_buffer {
public:
  image_buffer(const unsigned numX, const unsigned numY) :
    d_numX(numX), d_numY(numY),
    d_buffer(NULL)
  {
    d_buffer = new unsigned char[numX * numY * 4];
    if (d_buffer == NULL) {
      fprintf(stderr,"image_buffer::image_buffer(): Out of memory\n");
    }
  };
  
  ~image_buffer() { if (d_buffer) { delete d_buffer; d_buffer = NULL; } };

  // See if the image_buffer is set up properly
  bool doing_okay(void) const
  {
    return (d_buffer != NULL);
  }

  // Tell what the range is for the image.
  inline unsigned numX(void) const { return d_numX; }
  inline unsigned numY(void) const { return d_numY; }
  inline void read_range(int &minx, int &maxx, int &miny, int &maxy) const
  {
    minx = 0; maxx = numX()-1;
    miny = 0; maxy = numY()-1;
  }

  /// Read a pixel from the image; return true if the pixel
  // was in the image, false if it was not.
  inline bool read_pixel(unsigned x, unsigned y, unsigned rgba, unsigned char &result) const
  {
    if ( (x >= numX()) || (y >= numY()) || (rgba > 3) ) {
      return false;
    } else {
      result = d_buffer[rgba + 4 * ( x + numX() * y )];
      return true;
    }
  }

  /// Read a pixel from the image; Don't check boundaries.
  inline unsigned char read_pixel_nocheck(unsigned x, unsigned y, unsigned rgba) const
  {
    return d_buffer[rgba + 4 * ( x + numX() * y )];
  }

  // Do bilinear interpolation to read from the image, in order to
  // smoothly interpolate between pixel values.
  // All sorts of speed tweaks in here because it is in the inner loop for
  // the spot tracker and other codes.
  inline bool	read_pixel_bilerp(double x, double y, unsigned rgba, double &result) const
  {
    // The order of the following statements is optimized for speed.
    // The double version is used below for xlowfrac comp, ixlow also used later.
    // Slightly faster to explicitly compute both here to keep the answer around.
    double xlow = floor(x); int ixlow = (int)xlow;
    // The double version is used below for ylowfrac comp, ixlow also used later
    // Slightly faster to explicitly compute both here to keep the answer around.
    double ylow = floor(y); int iylow = (int)ylow;
    int ixhigh = ixlow+1;
    int iyhigh = iylow+1;
    double xhighfrac = x - xlow;
    double yhighfrac = y - ylow;
    double xlowfrac = 1.0 - xhighfrac;
    double ylowfrac = 1.0 - yhighfrac;
    unsigned char ll, lh, hl, hh;

    // Combining the if statements into one using || makes it slightly slower.
    // Interleaving the result calculation with the returns makes it slower.
    if (!read_pixel(ixlow, iylow, rgba, ll)) { return false; }
    if (!read_pixel(ixlow, iyhigh, rgba, lh)) { return false; }
    if (!read_pixel(ixhigh, iylow, rgba, hl)) { return false; }
    if (!read_pixel(ixhigh, iyhigh, rgba, hh)) { return false; }
    result = ll * xlowfrac * ylowfrac + 
	     lh * xlowfrac * yhighfrac +
	     hl * xhighfrac * ylowfrac +
	     hh * xhighfrac * yhighfrac;
    return true;
  };

  // Do bilinear interpolation to read from the image, in order to
  // smoothly interpolate between pixel values.
  // All sorts of speed tweaks in here because it is in the inner loop for
  // the spot tracker and other codes.
  // Does not check boundaries to make sure they are inside the image.
  inline double read_pixel_bilerp_nocheck(double x, double y, unsigned rgba) const
  {
    // The order of the following statements is optimized for speed.
    // The double version is used below for xlowfrac comp, ixlow also used later.
    // Slightly faster to explicitly compute both here to keep the answer around.
    double xlow = floor(x); int ixlow = (int)xlow;
    // The double version is used below for ylowfrac comp, ixlow also used later
    // Slightly faster to explicitly compute both here to keep the answer around.
    double ylow = floor(y); int iylow = (int)ylow;
    int ixhigh = ixlow+1;
    int iyhigh = iylow+1;
    double xhighfrac = x - xlow;
    double yhighfrac = y - ylow;
    double xlowfrac = 1.0 - xhighfrac;
    double ylowfrac = 1.0 - yhighfrac;

    // Combining the if statements into one using || makes it slightly slower.
    // Interleaving the result calculation with the returns makes it slower.
    return read_pixel_nocheck(ixlow, iylow, rgba) * xlowfrac * ylowfrac +
	   read_pixel_nocheck(ixlow, iyhigh, rgba) * xlowfrac * yhighfrac +
	   read_pixel_nocheck(ixhigh, iylow, rgba) * xhighfrac * ylowfrac +
	   read_pixel_nocheck(ixhigh, iyhigh, rgba) * xhighfrac * yhighfrac;
  };

  // Provide a pointer to the buffer directly.  This is a const pointer, so that
  // folks can't use it to modify the image data.
  inline const unsigned char *buffer(void) { return d_buffer; };

  // Write a value into one of the colors of an image.  Return true if the pixel
  // was in range, false if not.
  inline bool write_pixel(unsigned x, unsigned y, unsigned rgba, unsigned char val)
  {
    if ( (x >= numX()) || (y >= numY()) || (rgba > 3) ) {
      return false;
    } else {
      d_buffer[rgba + 4 * ( x + numX() * y )] = val;
      return true;
    }
  };

  // Write a value into one of the colors of an image.
  inline void write_pixel_nocheck(int x, int y, unsigned rgba, unsigned char val)
  {
    d_buffer[rgba + 4 * ( x + numX() * y )] = val;
  };

  // Recolor myself by multiplying by the passed-in unit colors.
  inline void recolor(const Color &color) {
    unsigned x, y;
    for (x = 0; x < d_numX; x++) {
      for (y = 0; y < d_numY; y++) {
	write_pixel_nocheck(x,y,0, read_pixel_nocheck(x,y,0) * color[0]);
	write_pixel_nocheck(x,y,1, read_pixel_nocheck(x,y,1) * color[1]);
	write_pixel_nocheck(x,y,2, read_pixel_nocheck(x,y,2) * color[2]);
	write_pixel_nocheck(x,y,3, read_pixel_nocheck(x,y,3) * color[3]);
      }
    }
  };

protected:
  unsigned d_numX;	      //< Number of pixels in X
  unsigned d_numY;	      //< Number of pixels in Y
  unsigned char *d_buffer;    //< Pointer to RGBA buffer for the pixel
};

//----------------------------------------------------------------------------
// Class to deal with storing a list of locations on the image

class Position {
public:
  Position(int x, int y) { d_x = x; d_y = y; };
  int x(void) const { return d_x; };
  int y(void) const { return d_y; };
  void set_x(int x) { d_x = x; };
  void set_y(int y) { d_y = y; };

protected:
  int   d_x;
  int   d_y;
};

//----------------------------------------------------------------------------
// Class to deal with storing a list of component IDs and volumes

class Component_volume {
public:
  Component_volume(double id, double volume) { d_id = id; d_volume = volume; };
  int id(void) const { return d_id; };
  double volume(void) const { return d_volume; };
  void set_id(int id) { d_id = id; };
  void set_volume(double volume) { d_volume = volume; };

protected:
  int   d_id;
  double   d_volume;
};

//------------------------------------------------------------------------
// Static global variables
static Tcl_Interp	        *g_tk_control_interp;
static bool			verbose = false;

const  unsigned			g_num_colors = 5;	  //< Number of colors defined for use
static Color			g_colors[g_num_colors] =  //< Colors to use for the layers
			   {  { 0.57, 0.57, 0.57, 1.0} ,  //< White scaled down to isoluminant
			      { 1.00, 0.00, 0.00, 1.0 } , //< Red
			      { 0.00, 1.00, 0.00, 1.0 } , //< Green
			      { 0.71, 0.71, 0.00, 1.0 } , //< Yellow scaled down to isoluminant
			      { 0.00, 0.00, 1.00, 1.0 } };//< Blue

static int			g_data_window = -1;		  //< Glut window displaying one data set
static unsigned			g_nColsData = 640;	  //< Number of columns in the display window
static unsigned			g_nRowsData = 480;	  //< Number of rows in the display window
static image_buffer		*g_data_image = NULL;	  //< Pointer to the storage for the data set
static BCGrid		        *g_current_data = NULL;	  //< Data set currently displayed in browse window
static const char               *g_plane_to_compute_from = "data";  //< Name of the plane to do our math based on.

static int			g_second_window = -1;     //< Glut window displaying the visualization
static unsigned			g_nColsSecond = 512;	  //< Number of columns in the display window
static unsigned			g_nRowsSecond = 512;	  //< Number of rows in the display window
static image_buffer		*g_second_image = NULL;	  //< Pointer to the storage for the visualization window

static int			g_preview_window = -1;    //< Glut window displaying a preview of the template
static unsigned			g_nColsPreview = 0;	  //< Number of columns in the template window
static unsigned			g_nRowsPreview = 0;	  //< Number of rows in the template window
static image_buffer		*g_preview_image = NULL;  //< Pointer to the storage for the template window
static BCGrid                   *g_preview_grid = NULL;   //< Grid to use for previewing the data

static int			g_splat_window = -1;      //< Glut window displaying a the best-fit template splatted for each
static unsigned			g_nColsSplat = 0;	  //< Number of columns in the splat window
static unsigned			g_nRowsSplat = 0;	  //< Number of rows in the splat window
static image_buffer		*g_splat_image = NULL;    //< Pointer to the storage for the splat window

static vector<BCGrid *>         g_template_images;        //< List of templates to test against.  One plane per template grid.

static bool			g_ready_to_display = false;       //< Don't unless we get an image
static bool			g_data_already_posted = false;    //< Posted redisplay since the last display?
static bool			g_second_already_posted = false;  //< Posted redisplay since the last display?
static bool			g_splat_already_posted = false;   //< Posted redisplay since the last display?

static int		        g_mousePressX, g_mousePressY;     //< Where the mouse was when the button was pressed
static int		        g_whichDragAction = 0;		  //< What action to take for mouse drag

static vector<Position>         g_flattenList;            //< Keep track of the points used to flatten the image
static double                   g_radius = 4;             //< Radius of the flatting detector

static vector<Component_volume> g_componentVolumes;       //< List of ids and volumes of existing components
static vector<Component_volume> g_thresholdedVolumes;     //< List of ids and volumes of existing components

//------------------------------------------------------------------------
// Forward function definitions
static void cleanup();
static void dirtyexit();
static void redisplay_second(void);
static void redisplay_data(void);
static void redisplay_splat(void);
static void myDataDisplayFunc(void);
static void myPreviewDisplayFunc(void);
static void mySplatDisplayFunc(void);
static bool perform_deline_filter(const BCPlane *from, BCPlane *to);
static bool calculate_median_filter(const string from_name, const string to_name);
static void calculate_flatten(const string from_name, const string to_name);
static void calculate_height_threshold(const string volume_plane_name, const string threshold_plane_name);
static void calculate_connected_components(const string threshold_plane_name, const string component_string_name);
static void calculate_component_volumes(const string volume_plane_name, const string component_plane_name, vector<Component_volume> &volumes);
static void calculate_thresholded_volumes(const string in_plane_name, const vector<Component_volume> &in_volumes,
                                          const double norm_min, double norm_max,
                                          const string out_plane_name, vector<Component_volume> &out_volumes);

//------------------------------------------------------------------------
// User interface controls
void  data_filename_changed(char *newvalue, void *);
void  template_filename_changed(char *newvalue, void *);
void  save_filename_changed(char *newvalue, void *);
void  handle_template_name_change(char *newvalue, void *userdata);
void  float_causes_redisplay(float newvalue, void *);
void  int_causes_redisplay(int newvalue, void *);

Tclvar_selector         g_template_filename("template_filename", NULL, NULL, "", template_filename_changed, NULL);
Tclvar_list_of_strings	g_template_names;
Tclvar_selector		g_template_name("template_preview","",&g_template_names,"none");
Tclvar_float_with_scale g_height_threshold("height_threshold","", 0, 1, 0.4, float_causes_redisplay);
Tclvar_float_with_scale g_dilation("dilation","", 0, 5, 0, float_causes_redisplay);
Tclvar_selector         g_save_template_filename("save_template_filename", NULL, NULL, "");

Tclvar_float g_min_norm_volume("min_norm_volume", 0, float_causes_redisplay);
Tclvar_float g_max_norm_volume("max_norm_volume", 1, float_causes_redisplay);

Tclvar_int_with_button	g_quit("quit",NULL);
Tclvar_int_with_button	g_save_cancel("save_cancel",NULL);

Tclvar_int_with_button	g_substrate_fit("substrate_fit","", 1, int_causes_redisplay);
Tclvar_int_with_button	g_median_filter("median_filter","", 1, int_causes_redisplay);
Tclvar_int_with_button	g_show_labels("show_labels","", 0, int_causes_redisplay);

Tclvar_float_with_scale g_angle_step("angle_step","", 1, 90, 10);
Tclvar_float_with_scale g_pixel_step("pixel_step","", 0.25, 5, 5.0);
Tclvar_int_with_button	g_save("save", NULL);
Tclvar_selector         g_save_file_name("savefilename", NULL, NULL, "", save_filename_changed, NULL);

// Used to display progress of the calculation so the user doesn't lose
// hope.
Tclvar_int_with_button  g_comp_id("comp_id", NULL);
Tclvar_int_with_button  g_comp_count("comp_count", NULL);
Tclvar_int_with_button  g_temp_id("temp_id", NULL);
Tclvar_int_with_button  g_temp_count("temp_count", NULL);


void drawStringAtXY(double x, double y, char *string)
{
  void *font = GLUT_BITMAP_TIMES_ROMAN_10;
  int len, i;

  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++) {
    glutBitmapCharacter(font, string[i]);
  }
}


// If the data filename becomes non-empty, then load a new data set and
// add it to the list of available data sets.

void  data_filename_changed(char *newvalue, void *)
{
  if (strlen(newvalue) <= 0) {
    return;
  }

  // We should only be opening one file now, not multiple data files
  if (g_current_data != NULL) {
    char command[1024];
    sprintf(command, "tk_messageBox -type ok -title Error -message {Can only open data file at startup.}");
    if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
        fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                g_tk_control_interp->result);
        return;
    }
    return;
  }

  // Open new data file.
  // Add each to the list of available data entries.
  if (verbose) {
    fprintf(stderr,"Reading data file %s\n", newvalue);
  }

  // XXX For some reason, this fails when opening topo files but works for
  // DI files and image files.  The nano app rebuilt from the same source
  // can open these files.

  TopoFile  fakeTopo;
  BCGrid *new_grid = new BCGrid();
  new_grid = new_grid->loadFile(newvalue, fakeTopo);
  if (new_grid == NULL) {
    fprintf(stderr,"Can't read data image %s\n", newvalue);
    return;
  }
  new_grid->head()->rename("data");
  g_current_data = new_grid;
  g_nColsData = g_current_data->numX();
  g_nRowsData = g_current_data->numY();

  //--------------------------------------------------------------------
  // Calculate the substrate-fit using Gary Bishop's magical method that
  // looks for the most-consistent substrate in Y.  This will actually
  // fix the effects of line-by-line flattening artifacts (where the
  // substrate is lower on lines that have more stuff on them) as well
  // as removing slopes in the image.

  // Make a new plane that will be used for each component in turn.  It will
  // get filled with the proper subset of the flattened plane for each component.
  BCPlane *substrate = g_current_data->addPlaneCopy(g_current_data->getPlaneByName("data"));
  substrate->rename("substrate");
  perform_deline_filter(new_grid->head(), substrate);

  //--------------------------------------------------------------------
  // Information relevant to the data display window has changed,
  // so redraw it.
  redisplay_data();

  //--------------------------------------------------------------------
  // Data relevant to the second display has changed.  Force redisplay.
  redisplay_second();
}

// If the template filename becomes non-empty, load a new template and
// add it to the list of available templates.

void  template_filename_changed(char *newvalue, void *)
{
  if (strlen(newvalue) <= 0) {
    return;
  }

  // Open new template file.
  // Add each to the list of available templates.
  if (verbose) {
    fprintf(stderr,"Reading data file %s\n", newvalue);
  }

  // XXX For some reason, this fails when opening topo files but works for
  // DI files and image files.  The nano app rebuilt from the same source
  // can open these files.

  TopoFile  fakeTopo;
  BCGrid *new_grid = new BCGrid();
  new_grid = new_grid->loadFile(newvalue, fakeTopo);
  if (new_grid == NULL) {
    fprintf(stderr,"Can't read template image %s\n", newvalue);
    return;
  }

  // Make the name out of everything past the last slash or backslash
  char *name_to_use = new char[strlen(newvalue)+1];
  const char *name_so_far = newvalue;
  const char *last_slash = strrchr(newvalue, '/') + 1;
  const char *last_back = strrchr(newvalue, '\\') + 1;
  if (last_slash > name_so_far) { name_so_far = last_slash; }
  if (last_back > name_so_far) { name_so_far = last_back; }
  strcpy(name_to_use, name_so_far);

  // Add the name to the list of available ones.
  g_template_names.Add_entry(name_to_use);

  // Change the name of the first (only) plane in the grid to
  // match the file name.
  new_grid->head()->rename(name_to_use);
  g_template_images.push_back(new_grid);

  // Cause this to be the template that is selected for preview.
  // Call the callback for when it changes
  g_template_name.set_tcl_change_callback(handle_template_name_change, NULL);
  g_template_name = name_to_use;

  // Cause the whole update pipeline to rerun because one of its inputs
  // has changed.
  float_causes_redisplay(0.0, NULL);
}

// If the template name being previewed becomes non-empty and is not
// "none", create a new OpenGL window to preview it in and set it to
// draw.

void  handle_template_name_change(char *newvalue, void *)
{
  // Destroy any existing window and buffer.
  if (g_preview_window != -1) {
    glutDestroyWindow(g_preview_window);
    g_preview_window = -1;
  }
  if (g_preview_image) { delete g_preview_image; g_preview_image = NULL; };

  // If we have an empty name or "none", then we're done.
  if ( (strlen(newvalue) <= 0) || (strcmp(newvalue, "none") == 0) ) {
    return;
  }

  // Look up the information about the template we're going to use.
  // We do this by going through all of the template grids looking for
  // one that has a plane with the name we are looking for.
  g_preview_grid = NULL;
  vector<BCGrid *>::iterator  g;
  for (g = g_template_images.begin(); g != g_template_images.end(); ++g) {
    if ( *( (*g)->head()->name() ) == string(newvalue) ) {
      g_preview_grid = *g;
      break;
    }
  }
  if (g_preview_grid == NULL) {
    return;
  }

  g_nColsPreview = g_preview_grid->numX();
  g_nRowsPreview = g_preview_grid->numY();

  // Create the buffers that Glut will use to display in the window.  This is allocating an
  // RGBA buffer.  It needs to be 4-byte aligned, so we allocated it as a group of
  // words and then cast it to the right type.  We're using RGBA rather than just RGB
  // because it also solves the 4-byte alignment problem caused by funky sizes of image
  // that are RGB images.
  g_preview_image = new image_buffer(g_nColsPreview, g_nRowsPreview);
  if ( (g_preview_image == NULL) || !g_preview_image->doing_okay() ) {
    fprintf(stderr,"Out of memory when allocating image!\n");
    fprintf(stderr,"  (Image is %u by %u)\n", g_nColsPreview, g_nRowsPreview);
    cleanup();
  }

  // Create a new window, positioned under the data window.
  glutInitWindowPosition(175, 10 + 60 + g_nRowsData);
  glutInitWindowSize(g_nColsPreview, g_nRowsPreview);
  g_preview_window = glutCreateWindow(newvalue);

  // Set the rendering callback for the window.
  glutDisplayFunc(myPreviewDisplayFunc);
}

void  int_causes_redisplay(int newvalue, void *)
{
  float_causes_redisplay(newvalue, NULL);
}

void  float_causes_redisplay(float newvalue, void *)
{
  // Figure out which plane to calculate based on
  if (g_substrate_fit) {
    g_plane_to_compute_from = "substrate";
  } else {
    g_plane_to_compute_from = "data";
  }

  // Calculation of various statistics may need to be redone, so do
  // that here.
  calculate_flatten(g_plane_to_compute_from, "flat");
  // Only the connected-components should be affected by the median
  // filter calculation, and it should happen before the threshold.
  // All other calculations should come straight from the flat plane.
  if (g_median_filter) {
    calculate_median_filter("flat", "median_filter");
    calculate_height_threshold("median_filter", "height_threshold");
  } else {
    calculate_height_threshold("flat", "height_threshold");
  }
  calculate_connected_components("height_threshold", "connected_components");
  calculate_component_volumes("flat","connected_components", g_componentVolumes);
  calculate_thresholded_volumes("connected_components", g_componentVolumes,
    g_min_norm_volume, g_max_norm_volume,
    "thresholded_components", g_thresholdedVolumes);

  // Information relevant to the data display window has changed,
  // so redraw it.
  redisplay_data();

  // Data relevant to the second display has changed.  Force redisplay.
  redisplay_second();

  // Data relevant to the splay display has changed.  Force redisplay.
  redisplay_splat();
}

// Helper function to copy from one plane into another, filling in only
// those entries not masked by a connected component, and setting all
// other values to zero.
// Assumptions:

static bool copy_with_threshold(const BCPlane *from, const BCPlane *comp, int comp_id, BCPlane *to)
{
  if ( (from->numX() != comp->numX()) || (from->numX() != to->numX()) ||
       (from->numY() != comp->numY()) || (from->numY() != to->numY()) ) {
    return false;
  }

  int x, y;
  for (x = 0; x < from->numX(); x++) {
    for (y = 0; y < from->numY(); y++) {
      if (comp->scaledValue(x,y) == comp_id) {
        to->setValue(x,y, from->scaledValue(x,y));
      } else {
        to->setValue(x,y, 0.0);
      }
    }
  }

  return true;
}

// Helper function to calculate the center of a component in a plane.
// Finds and returns the center of area in cx,cy.  If there is a problem,
// it returns false.  Each spot with the correct id value is counted
// with the same weight.

static bool find_center_of_component(double &cx, double &cy, const BCPlane *height, int id)
{
  cx = cy = 0.0;

  // Loop through the plane and count up the pixels with
  // the correct id.  When we're done, we'll divide by this.
  int x, y;
  unsigned sum = 0;
  for (x = 0; x < height->numX(); x++) {
    for (y = 0; y < height->numY(); y++) {
      double val = height->scaledValue(x,y);
      if (val == id) {
        sum++;
        cx += x;
        cy += y;
      }
    }
  }

  if (sum == 0) {
    return false;
  }

  cx /= sum;
  cy /= sum;
  return true;
}

// Helper function to calculate the center of mass of a plane.
// Finds and returns the center of mass in cx,cy.  If there is a problem,
// it returns false.

static bool find_center_of_mass(double &cx, double &cy, const BCPlane *height)
{
  cx = cy = 0.0;

  // Loop through the plane and count up the height at each
  // existing point times the X and Y, as well as the total
  // height.  When we're done, we divide each by the total.
  int x, y;
  double sum = 0.0;
  for (x = 0; x < height->numX(); x++) {
    for (y = 0; y < height->numY(); y++) {
      double val = height->scaledValue(x,y);
      sum += val;
      cx += val * x;
      cy += val * y;
    }
  }

  if (sum == 0) {
    return false;
  }

  cx /= sum;
  cy /= sum;
  return true;
}

static inline double square_diff(double a, double b)
{
  return (a-b)*(a-b);
}

// Helper function to compute the offset of the lower-left template
// corner in the target height plane's (h) space, the differential step in X
// and Y within the target space for steps in X and Y in the template
// space, and the relative volume size of template space vs. target
// space.  The template is to be translated in its own space so that
// (tx,ty) lies at the origin, then rotated by angle, then scaled
// so that nanometers in its space matches nanometers in the target
// space, then translated so that its origin lies at (cx,cy) in the
// target space.

static inline bool find_offset_scale_zoom(const BCPlane *h, double cx, double cy,
                                   const BCPlane *t, double tx, double ty, double angle,
                                   double &offset_x, double &offset_y,
                                   double &dx_target_for_x, double &dy_target_for_x,
                                   double &dx_target_for_y, double &dy_target_for_y,
                                   double &volume_scale)
{
  if ( (h == NULL) || (t == NULL) ) {
    fprintf(stderr,"find_offset_scale_zoom(): called with NULL plane\n");
    return false;
  }

  // How many nanometers per single-pixel step in each image
  double nm_per_x_in_t = (t->maxX() - t->minX()) / (t->numX()-1);
  double nm_per_y_in_t = (t->maxY() - t->minY()) / (t->numY()-1);
  double nm_per_x_in_h = (h->maxX() - h->minX()) / (h->numX()-1);
  double nm_per_y_in_h = (h->maxY() - h->minY()) / (h->numY()-1);

  // The normalized step sizes to take into account the rotation
  double nx_from_x =  cos( angle * M_PI/180 );
  double ny_from_x =  sin( angle * M_PI/180 );
  double nx_from_y = -sin( angle * M_PI/180 );
  double ny_from_y =  cos( angle * M_PI/180 );

  // Compute the step size (in x and y) in h's pixel units of a single-pixel
  // step in x and y in t's pixel units.  We get this by scaling the pixel
  // steps to nanometers in X and Y, then rotating by the angle, then scaling
  // the result from nanometers to pixels in h's image.
  double nm_per_x_from_x = nm_per_x_in_t * nx_from_x;
  double nm_per_x_from_y = nm_per_y_in_t * nx_from_y;
  double nm_per_y_from_x = nm_per_x_in_t * ny_from_x;
  double nm_per_y_from_y = nm_per_y_in_t * ny_from_y;

  dx_target_for_x = nm_per_x_from_x / nm_per_x_in_h;
  dx_target_for_y = nm_per_y_from_x / nm_per_x_in_h;
  dy_target_for_x = nm_per_x_from_y / nm_per_y_in_h;
  dy_target_for_y = nm_per_y_from_y / nm_per_y_in_h;

  // Compute the location for the (0,0) pixel from t's image in h's image,
  // in pixel units.  Keep in mind that the origin needs to be shifted in
  // t's space to (tx,ty).
  offset_x = dx_target_for_x * (-tx) + dx_target_for_y * (-ty);
  offset_y = dy_target_for_x * (-tx) + dy_target_for_y * (-ty);

  // Translate the whole deal so that it is centered around cx,cy
  offset_x += cx;
  offset_y += cy;

  // Compute the total volume scale
  volume_scale = (nm_per_x_in_t * nm_per_y_in_t) / (nm_per_x_in_h * nm_per_y_in_h);
  return true;
}


// Adds the non-zero template elements into the plane whose pointer is passed.
// Weights each point in the template by its area, so that the total volume is
// preserved.  Does bilerp splat into for nearest neighbors for each point.
// Returns false on failure, true on success.

static bool splat_template_into_plane(const BCPlane *from, BCPlane *to, double best_x, double best_y, double best_angle)
{
  if ( (from == NULL) || (to == NULL) ) {
    fprintf(stderr,"splat_template_into_plane(): called with NULL plane\n");
    return false;
  }

  double tx, ty;
  find_center_of_mass(tx, ty, from);

  double x0,y0, volume_scale, dx_from_x, dx_from_y, dy_from_x, dy_from_y;
  find_offset_scale_zoom(to, best_x,best_y, from, tx,ty,best_angle, x0,y0, dx_from_x,dy_from_x, dx_from_y,dy_from_y, volume_scale);

  int x,y;
  for (x = 0; x < from->numX(); x++) {
    for (y = 0; y < from->numY(); y++) {
      // Find the coordinate in to's space that tells where this splat goes.
      double sx = x0 + x * dx_from_x + y * dx_from_y;
      double sy = y0 + x * dy_from_x + y * dy_from_y;

      // Find the weightings for the four corners, based on the distance from
      // each and also based on the volume scale between the template and the
      // plane we're splatting into.
      int xlower = floor(sx);
      int ylower = floor(sy);
      int xhigher = xlower + 1;
      int yhigher = ylower + 1;
      double xfrac = sx - xlower;
      double yfrac = sy - ylower;
      double wll = (1 - xfrac) * (1 - yfrac) * volume_scale;
      double whl = (    xfrac) * (1 - yfrac) * volume_scale;
      double wlh = (1 - xfrac) * (    yfrac) * volume_scale;
      double whh = (    xfrac) * (    yfrac) * volume_scale;

      // Make sure we're not writing past the edge of the plane.
      if ( (xlower < 0) || (ylower < 0) || (xhigher >= to->numX()) || (yhigher >= to->numY()) ) {
        continue;
      }

      // Add the values into the four corner pixels.
      to->setValue(xlower, ylower, to->scaledValue(xlower,ylower) + from->scaledValue(x,y)*wll);
      to->setValue(xhigher, ylower, to->scaledValue(xhigher,ylower) + from->scaledValue(x,y)*whl);
      to->setValue(xlower, yhigher, to->scaledValue(xlower,yhigher) + from->scaledValue(x,y)*wlh);
      to->setValue(xhigher, yhigher, to->scaledValue(xhigher,yhigher) + from->scaledValue(x,y)*whh);
    }
  }

  // XXX Problem: If the volume scale is greater than one, we need to rasterize the
  // values into the plane because each volume element will fill more than one, so
  // we end up skipping points when stepping in the template.  We could also do some
  // sort of post-processing blur to fix this, maybe.

  return true;
}

// Computes the sum-of-squares difference between each pixel in the template (t)
// and an interpolated position in the height image (h).  First translates so that
// the center of the template (located at tx,ty) is located at the desired
// center of the height image (cx,cy).  Steps in the image are taken in the
// direction angle (in degrees) away from the basis.
//
// The relative sizes of the scan areas of the two images is taken into account,
// so that the same-sized steps are taken in the template image (single-pixel
// steps) and in the height image (computed-distance steps).

static double compute_template_error(const BCPlane *h, double cx, double cy, double angle, const BCPlane *t, double tx, double ty)
{
  double err = 0.0;

  double x0,y0, volume_scale, dx_from_x, dx_from_y, dy_from_x, dy_from_y;
  find_offset_scale_zoom(h, cx,cy, t, tx,ty,angle, x0,y0, dx_from_x,dy_from_x, dx_from_y,dy_from_y, volume_scale);

  // Loop through the pixels in t's image and sum up the squared errors.
  // Pixels that fall outside h's image should have their values counted
  // as zero (don't make the difference zero or the template just scoots off
  // the screen).  The errors should be VOLUME errors, so need to be scaled
  // by t's pixel sizes in X and Y.
  int x,y;
  double xmax = h->numX()-1;
  double ymax = h->numY()-1;
  for (x = 0; x < t->numX(); x++) {
    for (y = 0; y < t->numY(); y++) {
      double hx = x0 + dx_from_x*x + dx_from_y*y;
      double hy = y0 + dy_from_x*x + dy_from_y*y;
      if ( (hx >= 0) && (hx <= xmax) && (hy >= 0) && (hy <= ymax) ) {
        err += ( square_diff( t->scaledValue(x,y), h->scaledValue(hx,hy) ) ) * volume_scale;
      } else {
        err += t->scaledValue(x,y) * volume_scale;
      }
    }
  }

  return err;
}

// Attempts various position and orientations for the template image against the
// height image and finds the minimum sum-of-square error terms between them.
// This is pixel-wise comparison for every point in the template, with the location
// of each pixel interpolated into the height image.
// Parameters:
//    dxy is the step size to take in pixels in the height image
//    danlge is the rotational step size in degrees
// Return:
//    Return value is the error at the best fit
//    best_x, best_y, best_angle report where this fit occured.

static double compute_min_template_error(const BCPlane *height, double cx, double cy, const BCPlane *t, double dxy, double dangle,
                                         double &best_x, double &best_y, double &best_angle)
{
  // Compute the center of mass of the template, which we'll line up with the
  // center of mass of the image to start with.
  double tx, ty;
  find_center_of_mass(tx, ty, t);

  // Figure out how much to shift at most.  This should be the number of pixels
  // in the height image that corresponds to half the width of the template.
  double nm_per_x_in_t = (t->maxX() - t->minX()) / (t->numX()-1);
  double nm_per_y_in_t = (t->maxY() - t->minY()) / (t->numY()-1);
  double nm_per_x_in_h = (height->maxX() - height->minX()) / (height->numX()-1);
  double nm_per_y_in_h = (height->maxY() - height->minY()) / (height->numY()-1);

  double max_shift_in_x = (t->numX()/2) * nm_per_x_in_t / nm_per_x_in_h;
  double max_shift_in_y = (t->numY()/2) * nm_per_y_in_t / nm_per_y_in_h;

  // Okay, now we start to go winding through strange spaces.  Remember:
  //  (cx,cy) is the center in pixel coordinates within height's pixel space
  //  (tx,ty) is the center in pixel coordinates within t's pixel space.
  //  dxy is the step size in height's pixel units to shift the template (t)
  //  dangle is the step size in degrees to try reorienting the template (t)

  double angle;
  double dx, dy;
  double min_err = compute_template_error(height, cx,cy, 0.0, t, tx,ty);
  best_x = cx;
  best_y = cy;
  best_angle = 0;
  for (angle = 0; angle < 360; angle += dangle) {
    // Check the non-shifted version.
    double err = compute_template_error(height, cx,cy, angle, t, tx,ty);
    if (err < min_err) {
      min_err = err;
      best_x = cx;
      best_y = cy;
      best_angle = angle;
    }

    // Check all of the shifted versions.  We shift away from
    // the cx,cy position by the dx,dy amounts (all of these
    // are in height's pixel units).  We tell it how far to
    // shift t in its own units to center it on the center of
    // mass.
    for (dx = dxy; dx <= max_shift_in_x; dx += dxy) {
      for (dy = dxy; dy <= max_shift_in_y; dy += dxy) {
        // Be sure to check +/- shift in X and Y.  We do it
        // this way so that we test symmetrically around the
        // non-shifted version.
        err = compute_template_error(height, cx+dx,cy+dy, angle, t, tx,ty);
        if (err < min_err) {
          min_err = err;
          best_x = cx+dx;
          best_y = cy+dy;
          best_angle = angle;
        }
        err = compute_template_error(height, cx-dx,cy+dy, angle, t, tx,ty);
        if (err < min_err) {
          min_err = err;
          best_x = cx-dx;
          best_y = cy+dy;
          best_angle = angle;
        }
        err = compute_template_error(height, cx+dx,cy-dy, angle, t, tx,ty);
        if (err < min_err) {
          min_err = err;
          best_x = cx+dx;
          best_y = cy-dy;
          best_angle = angle;
        }
        err = compute_template_error(height, cx-dx,cy-dy, angle, t, tx,ty);
        if (err < min_err) {
          min_err = err;
          best_x = cx-dx;
          best_y = cy-dy;
          best_angle = angle;
        }
      }
    }
  }

  return min_err;
}

// If the save filename becomes non-empty, calculate the optimal fits and
// save the results in the file passed in.  If there are no templates,
// then just save the volume info and move the volume data into the image.

void  save_filename_changed(char *newvalue, void *)
{
  if (strlen(newvalue) <= 0) {
    return;
  }

  // Open new template file.
  // Add each to the list of available templates.
  if (verbose) {
    fprintf(stderr,"Save CSV file %s\n", newvalue);
  }

  // Make sure that the file does not exist by deleting it if it does.
  // The Tcl code had a dialog box that asked the user if they wanted
  // to overwrite, so this is "safe."
  FILE *in_the_way;
  if ( (in_the_way = fopen(newvalue, "r")) != NULL) {
    fclose(in_the_way);
    int err;
    if ( (err=remove(newvalue)) != 0) {
      fprintf(stderr,"Error: could not delete existing logfile %s\n", newvalue);
      perror("   Reason");
      cleanup();
      exit(-1);
    }
  }

  // Open the file if we can.
  FILE *csv_file = NULL;
  if ( NULL == (csv_file = fopen(newvalue, "w")) ) {
    fprintf(stderr,"Cannot open CSV file for writing: %s\n", newvalue);
    g_save = 0;
    return;
  }
  
  // Put the header into the file.
  fprintf(csv_file, "Object_ID,X,Y,Volume,Best_template,Best_error,Best_angle,Second_best_template,Second_best_error,Second_best_angle\n");

  // Create a progress display panel.
  if (Tcl_Eval(g_tk_control_interp, "track_progress") != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", "track_progress",
                  g_tk_control_interp->result);
          return;
  }

  // Make a new plane that will be used for each component in turn.  It will
  // get filled with the proper subset of the flattened plane for each component.
  BCPlane *subset = g_current_data->addPlaneCopy(g_current_data->getPlaneByName("flat"));
  subset->rename("subset");

  // Create a splat plane into which all of the best-fit templates will be splatted
  // for comparison with the data.  If there is already a splat plane, re-use it.  In
  // either case, we set it to zero so that we can splat the stuff into it.
  BCPlane *splat = g_current_data->getPlaneByName("splat");
  if (splat == NULL) {
    splat = g_current_data->addPlaneCopy(g_current_data->getPlaneByName("flat"));
    if (splat == NULL) {
      fprintf(stderr,"Cannot create splat plane\n");
      cleanup();
    }
    splat->rename("splat");
  }
  int x,y;
  for (x = 0; x < splat->numX(); x++) {
    for (y = 0; y < splat->numY(); y++) {
      splat->setValue(x,y, 0.0);
    }
  }

  // Look through each template for each region and find the best match.
  // Keep track of the top two fits so that we can report them.
  unsigned comp;
  for (comp = 0; comp < g_thresholdedVolumes.size(); comp++) {

    // Let the progress bar update itself and the user cancel if they want.
    g_comp_id = comp + 1;
    g_comp_count = g_thresholdedVolumes.size();
    Tclvar_mainloop();
    while (Tk_DoOneEvent(TK_DONT_WAIT)) {};

    // Copy the proper subset of the flattened plane into the subset plane.
    copy_with_threshold(g_current_data->getPlaneByName("flat"),
      g_current_data->getPlaneByName("thresholded_components"), g_thresholdedVolumes[comp].id(),
      subset);

    // Calculate the center of mass for the component, which will
    // be stored in the statistics and used to determine the initial
    // offset for the template.
    double cx, cy;
    find_center_of_mass(cx, cy, subset);

    // Run through each template and find the best and second-best matches
    unsigned t;
    unsigned best_id = 0, second_id = 0;
    double   best_error = 1e10, second_error = 1e10;
    double   best_x = 0, second_x;
    double   best_y = 0, second_y;
    double   best_angle = 0, second_angle;
    for (t = 0; t < g_template_images.size(); t++) {
      // Update the progress panel.  Push changes to
      // Tcl and update them in the display.
      g_temp_id = t + 1;
      g_temp_count = g_template_images.size();
      Tclvar_mainloop();
      while (Tk_DoOneEvent(TK_DONT_WAIT)) {};

      // If we've been asked to cancel, then we break out in the middle here
      // after cleaning things up.
      if (g_save_cancel) {
        // We want to drop completely out of the nested loops here and
        // do the teardown parts.  Proper use of "goto"!
        goto done_checking_components;
        break;
      }

      // Compute the fit.
      double fit_x, fit_y, fit_angle;
      double fit = compute_min_template_error(subset, cx, cy, g_template_images[t]->head(), g_pixel_step, g_angle_step,
        fit_x, fit_y, fit_angle);
      if (fit < best_error) {
        second_error = best_error;
        second_id = best_id;
        second_x = best_x;
        second_y = best_y;
        second_angle = best_angle;

        best_error = fit;
        best_id = t;
        best_x = fit_x;
        best_y = fit_y;
        best_angle = fit_angle;
      } else if (fit < second_error) {
        second_error = fit;
        second_id = t;
        second_x = fit_x;
        second_y = fit_y;
        second_angle = fit_angle;
      }
    }

    // Fill in the parts of the output that don't depend on the
    // templates.
    fprintf(csv_file, "%d,%lg,%lg,%lg",
      static_cast<int>(g_thresholdedVolumes[comp].id()),
      cx,cy,
      g_thresholdedVolumes[comp].volume() );

    // If we have templates, fill in the values for the best and
    // second best, and splat the best into the data set.  Otherwise,
    // splat the subset image itself into the data set.
    if (g_template_images.size() > 0) {
      // Fill in the values for this component.
      fprintf(csv_file, ",%s,%lg,%lg,%s,%lg,%lg",
        g_template_images[best_id]->head()->name()->c_str(), best_error, best_angle,
        g_template_images[second_id]->head()->name()->c_str(), second_error, second_angle);

      // Splat this template at this position and angle into the splat image.
      splat_template_into_plane(g_template_images[best_id]->head(), splat, best_x, best_y, best_angle);
    } else {
      // Splat the subset image right where it came from.
      splat_template_into_plane(subset, splat, cx, cy,0);
    }

    // End the line for this component.
    fprintf(csv_file, "\n");
  }

  // Target for the goto above that breaks out of the loop if the user cancels.
  done_checking_components:

  g_current_data->removePlane("subset");

  // Close the progress display panel.
  if (Tcl_Eval(g_tk_control_interp, "done_tracking_progress") != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", "track_progress",
                  g_tk_control_interp->result);
          return;
  }

  // We're done with the file.
  fclose(csv_file);

  // Turn off saving button.
  g_save = 0;

  // Cause the whole update pipeline to rerun because one of its inputs
  // has changed.
  float_causes_redisplay(0.0, NULL);
}

// More robust fgets() function that strips off the newlines.
// Useful for reading configuration files.
char *read_line_from_file(char *buffer, int maxlen, FILE *f)
{
  char *ret = fgets(buffer, maxlen, f);
  if (ret == NULL) {
    return NULL;
  }

  // Remove any stupid extra characters at the end of the line put in by
  // Windows or some other OS.  There may be up to two of them, so we do
  // this twice to get rid of both if they are there.
  char *last = &buffer[strlen(buffer)-1];
  if ( (*last == '\012') || (*last == '\015') ) { *last = '\0'; }
  last = &buffer[strlen(buffer)-1];
  if ( (*last == '\012') || (*last == '\015') ) { *last = '\0'; }
  return ret;
}

// This is called when someone kills the task by closing Glut or some
// other means we don't have control over.
static void  dirtyexit(void)
{
  static bool did_already = false;

  if (did_already) {
    return;
  } else {
    did_already = true;
  }

  // Done with the camera and other objects.
  printf("Exiting\n");

  //------------------------------------------------------------------------
  // Free up memory buffers
  if (g_data_image) { delete g_data_image; g_data_image = NULL; };
  if (g_second_image) { delete g_second_image; g_second_image = NULL; };
  if (g_preview_image) { delete g_preview_image; g_preview_image = NULL; };
}

static void  cleanup(void)
{
  static bool cleaned_up_already = false;

  if (cleaned_up_already) {
    return;
  } else {
    cleaned_up_already = true;
  }

  // Done with the camera and other objects.
  printf("Cleanly ");

  // Do the dirty-exit stuff, then do any VRPN cleanup needed (for
  // some reason, it cannot be done in dirty in the program we copied
  // from).
  dirtyexit();
}

static  void put_plane_into_buffer(BCPlane *plane, image_buffer *b)
{
    double  offset = plane->minValue();
    double  scale = 255 / (plane->maxValue() - plane->minValue());

    // Copy pixels into the image buffer.
    int r,c;
    unsigned char uns_pix;
    for (r = 0; r < plane->numY(); r++) {
      for (c = 0; c < plane->numX(); c++) {
	uns_pix = scale * (plane->scaledValue(c,r) - offset);
        double value = uns_pix;
        if (value > 255) { value = 255; }
        if (value < 0) { value = 0; }

	// This writes pixels
	// from the first channel into all colors of the image.  It uses
	// RGBA so that we don't have to worry about byte-alignment problems
	// that plagued us when using RGB pixels.
	b->write_pixel(c,r,0, value);
	b->write_pixel(c,r,1, value);
	b->write_pixel(c,r,2, value);
	b->write_pixel(c,r,3, 255);
      }
    }

/* Because we are using the BCGrid class to read things now, we can go ahead
   and not flip -- it will have already flipped.
    // Store the pixels from the image into the frame buffer
    // so that they cover the entire image (starting from upper-left
    // corner, which is at (-1,1) and using a pixel zoom factor that
    // is negative in Y so that we flip the image because everyone but
    // OpenGl indexes things from the upper-left corner).
//    glRasterPos2f(-1, 1);
//    glPixelZoom(1.0, -1.0);
*/
    glRasterPos2f(-1,-1);
    glDrawPixels(plane->numX(), plane->numY(), GL_RGBA, GL_UNSIGNED_BYTE, b->buffer());
}

void myDataDisplayFunc(void)
{
  // Clear the window and prepare to draw in the back buffer
  glDrawBuffer(GL_BACK);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  // If we have a current data set specified, then fill it into the buffer
  // and draw it.
  if (g_current_data) {
    // Read from the first (only) plane in the grid.
    // Figure out the scale and offset to apply to get this into the range 0..255
    BCPlane *plane = g_current_data->getPlaneByName(g_plane_to_compute_from);
    put_plane_into_buffer(plane, g_data_image);
  }

  // If we have any flattening points specified, draw markers for them
  unsigned i;
  for (i = 0; i < g_flattenList.size(); i++) {
    // Normalize center and radius so that they match the coordinates
    // (-1..1) in X and Y.
    double  x = -1.0 + g_flattenList[i].x() * (2.0/g_nColsData);
    double  y = -1.0 + g_flattenList[i].y() * (2.0/g_nRowsData);
    double  dx = g_radius * (2.0/g_nColsData);
    double  dy = g_radius * (2.0/g_nRowsData);
    
    glColor3f(1,0,0);
    glBegin(GL_LINES);
      glVertex2f(x-dx,y);
      glVertex2f(x+dx,y);
      glVertex2f(x,y-dy);
      glVertex2f(x,y+dy);
    glEnd();
  }

  // See if OpenGL had any problems
  GLenum ret;
  if ( (ret = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "myDataDisplayFunc(): OpenGL error detected\n");
  }

  // Swap buffers so we can see it.
  glutSwapBuffers();

  // Have no longer posted redisplay since the last display.
  g_data_already_posted = false;
}

void mySecondDisplayFunc(void)
{

  // Clear the window and prepare to draw in the back buffer
  glDrawBuffer(GL_BACK);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  if (!g_ready_to_display) {
    glutSwapBuffers();
    g_second_already_posted = false;
    return;
  }

  //----------------------------------------------------------------
  // If we don't have three flattening points, we don't display
  // anything here.
  if (g_flattenList.size() != 3) {
    glutSwapBuffers();
    g_second_already_posted = false;
    return;
  }

  //----------------------------------------------------------------
  // Read from the flattened plane in the grid (if it exists).
  // Figure out the scale and offset to apply to get this into the range 0..255
  BCPlane *plane = g_current_data->getPlaneByName("flat");
  BCPlane *height_threshold = g_current_data->getPlaneByName("height_threshold");
  BCPlane *thresholded_component = g_current_data->getPlaneByName("thresholded_components");
  if ( (plane == NULL) || (height_threshold == NULL) || (thresholded_component == NULL) ) {
    glutSwapBuffers();
    g_second_already_posted = false;
    return;
  }
  double  offset = plane->minValue();
  double  scale = 255 / (plane->maxValue() - plane->minValue());

  // Copy pixels into the image buffer.
  unsigned r,c;
  unsigned char uns_pix;
  for (r = 0; r < g_nRowsData; r++) {
    for (c = 0; c < g_nColsData; c++) {
      uns_pix = scale * (plane->scaledValue(c,r) - offset);
      double value = uns_pix;
      if (value > 255) { value = 255; }
      if (value < 0) { value = 0; }

      // This writes pixels
      // from the first channel into all colors of the image.  It uses
      // RGBA so that we don't have to worry about byte-alignment problems
      // that plagued us when using RGB pixels.
      // If this pixel is below the threshold, draw it in blue.
      if (height_threshold->scaledValue(c,r)) {
        // If the pixel is connected to the boundary or is part of
        // a component that was too close to another, draw it in
        // a brighter blue.  If the component is outside the
        // volume range, draw it in medium red.
        if (thresholded_component->scaledValue(c,r) == -4) {
          g_second_image->write_pixel(c,r,0, 128.0);
          g_second_image->write_pixel(c,r,1, 0.0);
          g_second_image->write_pixel(c,r,2, 0.0);
        } else if (thresholded_component->scaledValue(c,r) < 0) {
          g_second_image->write_pixel(c,r,0, 0.0);
          g_second_image->write_pixel(c,r,1, 0.0);
          g_second_image->write_pixel(c,r,2, 128.0);
        } else {
          g_second_image->write_pixel(c,r,0, value);
          g_second_image->write_pixel(c,r,1, value);
          g_second_image->write_pixel(c,r,2, value);
        }
      } else {
        g_second_image->write_pixel(c,r,0, 0.0);
        g_second_image->write_pixel(c,r,1, 0.0);
        g_second_image->write_pixel(c,r,2, 64.0);
      }
      g_second_image->write_pixel(c,r,3, 255);
    }
  }

/* Because we are using the BCGrid class to read things now, we can go ahead
   and not flip -- it will have already flipped.
  // Store the pixels from the image into the frame buffer
  // so that they cover the entire image (starting from upper-left
  // corner, which is at (-1,1) and using a pixel zoom factor that
  // is negative in Y so that we flip the image because everyone but
  // OpenGl indexes things from the upper-left corner).
//  glRasterPos2f(-1, 1);
//  glPixelZoom(1.0, -1.0);
*/
  glDrawPixels(g_nColsSecond, g_nRowsSecond, GL_RGBA, GL_UNSIGNED_BYTE, g_second_image->buffer());

  // See if OpenGL had any problems
  GLenum ret;
  if ( (ret = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "myDataDisplayFunc(): OpenGL error detected\n");
  }

  // Swap buffers so we can see it.
  glutSwapBuffers();
  // Have no longer posted redisplay since the last display.
  g_second_already_posted = false;
}

static void myPreviewDisplayFunc(void)
{
  // Clear the window and prepare to draw in the back buffer
  glDrawBuffer(GL_BACK);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  // If we have a preview data set specified, then fill it into the buffer
  // and draw it.
  if (g_preview_grid && g_preview_image) {
    // Read from the first (only) plane in the grid.
    // Figure out the scale and offset to apply to get this into the range 0..255
    BCPlane *plane = g_preview_grid->head();
    put_plane_into_buffer(plane, g_preview_image);
  }

  // See if OpenGL had any problems
  GLenum ret;
  if ( (ret = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "myPreviewDisplayFunc(): OpenGL error detected\n");
  }

  // Swap buffers so we can see it.
  glutSwapBuffers();
}

void mySplatDisplayFunc(void)
{
  // Clear the window and prepare to draw in the back buffer
  glDrawBuffer(GL_BACK);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  // If we have a splat data set specified, then fill it into the buffer
  // and draw it.
  if (g_current_data) {
    // Read from the splat plane in the grid.
    // Figure out the scale and offset to apply to get this into the range 0..255
    BCPlane *plane = g_current_data->getPlaneByName("splat");
    if (plane) {
      put_plane_into_buffer(plane, g_splat_image);
    }
  }

  // If we have objects, label them with their ID.
  if (g_show_labels) {
    char  label[50];
    unsigned i;
    BCPlane *components = g_current_data->getPlaneByName("thresholded_components");
    if (components != NULL) {
      for (i = 0; i < g_thresholdedVolumes.size(); i++) {
        int id = g_thresholdedVolumes[i].id();
        double cx,cy;

        // Find the center of mass, so we know where to put the label
        find_center_of_component(cx, cy, components, id);

        // Draw the label slightly offset from the center of mass
        sprintf(label, "%d", id);
        glColor3f(1,1,1);
        double offsetx = -7;
        double offsety = offsetx;
        double wx = -1 + 2*( (cx+offsetx) / (g_nColsSplat-1));
        double wy = -1 + 2*( (cy+offsety) / (g_nRowsSplat-1));
        drawStringAtXY(wx, wy, label);
      }
    }
  }

  // See if OpenGL had any problems
  GLenum ret;
  if ( (ret = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "myDataDisplayFunc(): OpenGL error detected\n");
  }

  // Swap buffers so we can see it.
  glutSwapBuffers();

  // Have no longer posted redisplay since the last display.
  g_splat_already_posted = false;
}

static void redisplay_second(void)
{
  // Post a redisplay so that Glut will draw the new image
  if (!g_second_already_posted && (g_second_window != -1)) {
    glutSetWindow(g_second_window);
    g_second_already_posted = true;
    glutPostRedisplay();
  }
}

static void redisplay_data(void)
{
  if (!g_data_already_posted && (g_data_window != -1)) {
    glutSetWindow(g_data_window);
    g_data_already_posted = true;
    glutPostRedisplay();
  }
}

static void redisplay_splat(void)
{
  if (!g_splat_already_posted && (g_splat_window != -1)) {
    glutSetWindow(g_splat_window);
    g_splat_already_posted = true;
    glutPostRedisplay();
  }
}

void myDataReshapeFunc(int width, int height)
{
  //------------------------------------------------------------------------
  // We don't allow reshaping of the
  // window.  We avoid it by calling a reshape function that puts it
  // back every time the user tries to change it.  Otherwise, reset the
  // viewport so that it covers the new window size.
  glutReshapeWindow(g_nColsData, g_nRowsData);
  redisplay_data();
}

void mySecondReshapeFunc(int width, int height)
{
  //------------------------------------------------------------------------
  // We don't allow reshaping of the
  // window.  We avoid it by calling a reshape function that puts it
  // back every time the user tries to change it.  Otherwise, reset the
  // viewport so that it covers the new window size.
  glutReshapeWindow(g_nColsSecond, g_nRowsSecond);
  redisplay_data();
}

void mySplatReshapeFunc(int width, int height)
{
  //------------------------------------------------------------------------
  // We don't allow reshaping of the
  // window.  We avoid it by calling a reshape function that puts it
  // back every time the user tries to change it.  Otherwise, reset the
  // viewport so that it covers the new window size.
  glutReshapeWindow(g_nColsSplat, g_nRowsSplat);
  redisplay_data();
}

// The status functions are called when a window is hidden or brought into
// view. Whenever the window is exposed at all, redraw.
void myDataStatusFunc(int status)
{
  if ( (status == GLUT_FULLY_RETAINED) || (status == GLUT_PARTIALLY_RETAINED) ) {
    myDataDisplayFunc();
  }
}

// The status functions are called when a window is hidden or brought into
// view. Whenever the window is exposed at all, redraw.
void mySecondStatusFunc(int status)
{
  if ( (status == GLUT_FULLY_RETAINED) || (status == GLUT_PARTIALLY_RETAINED) ) {
    // XXX For some reason, this causes the window to redisplay twice...
    //mySecondDisplayFunc();
  }
}

// The status functions are called when a window is hidden or brought into
// view. Whenever the window is exposed at all, redraw.
void mySplatStatusFunc(int status)
{
  if ( (status == GLUT_FULLY_RETAINED) || (status == GLUT_PARTIALLY_RETAINED) ) {
    // XXX For some reason, this causes the window to redisplay twice...
    //mySecondDisplayFunc();
  }
}

void myIdleFunc(void)
{
  //------------------------------------------------------------
  // This must be done in any Tcl app, to allow Tcl/Tk to handle
  // events.  This must happen at the beginning of the idle function
  // so that the camera-capture and video-display routines are
  // using the same values for the global parameters.

  while (Tk_DoOneEvent(TK_DONT_WAIT)) {};

  //------------------------------------------------------------
  // This is called once every time through the main loop.  It
  // pushes changes in the C variables over to Tcl.

  if (Tclvar_mainloop()) {
    fprintf(stderr,"Tclvar Mainloop failed\n");
  }

  //------------------------------------------------------------
  // Sleep for 10ms to avoid eating the whole processor just doing
  // busy-waiting.  This function is defined at the top of this file
  // for non-Windows compilation.
  Sleep(10);

  //------------------------------------------------------------
  // Time to quit?
  if (g_quit) {
    cleanup();
    exit(0);
  }
}

//--------------------------------------------------------------------------
// Find the average value within radius of the x,y coordinate within
// the specified plane.
static double find_average_value_around(const BCPlane *from, int xcenter, int ycenter, double g_radius)
{
  unsigned count = 0;               //< How many pixels have been included
  int dx, dy;                       //< Loop around the center point
  double r2 = g_radius * g_radius;  //< Compare to see if point is within radius
  double sum = 0.0;                 //< Used to run up the sum during averaging

  for (dx = -g_radius; dx <= g_radius; dx++) {
    for (dy = -g_radius; dy <= g_radius; dy++) {
      if ( dx*dx + dy*dy <= r2 ) {
        int x = xcenter + dx;
        int y = ycenter + dy;
        if ( (x >= 0) && (y >= 0) && (x < from->numX()) && (y < from->numY()) ) {
          count++;
          sum += from->scaledValue(x,y);
        }
      }
    }
  }

  if (count == 0) {
    return 0.0;
  } else {
    return sum / count;
  }
}

//--------------------------------------------------------------------------
// Calculate the "flat" plane, which is the input plane sheared so that
// the three points selected all lie at zero.  Average regions around the
// three points to get a more robust estimate.
// The routine must be robust, so that it can be called when any or all
// of its required inputs don't exist.
static void calculate_flatten(const string from_name, const string to_name)
{
  //------------------------------------------------------------
  // If we don't have an input plane, bail.
  if (g_current_data == NULL) { return; }
  BCPlane *from = g_current_data->getPlaneByName(from_name);
  if (from == NULL) { return; }

  //------------------------------------------------------------
  // If we don't have three points to flatten from, bail.
  if (g_flattenList.size() != 3) { return; }

  //------------------------------------------------------------
  // See if the flattened plane exists already. If not, then
  // create it.
  BCPlane *to = g_current_data->getPlaneByName(to_name);
  if (to == NULL) {
    to = g_current_data->addPlaneCopy(from);
    if (to == NULL) {
      fprintf(stderr,"calculate_flatten(): Could not create new plane\n");
      exit(-1);
    }
    to->rename(to_name);
  }

  //------------------------------------------------------------
  // Calculate the three averaged points, with X and Y specified
  // by the pointers and Z specified by the average of surrounding
  // pixel values.
  double x1 = g_flattenList[0].x();
  double y1 = g_flattenList[0].y();
  double z1 = find_average_value_around(from, x1, y1, g_radius);

  double x2 = g_flattenList[1].x();
  double y2 = g_flattenList[1].y();
  double z2 = find_average_value_around(from, x2, y2, g_radius);

  double x3 = g_flattenList[2].x();
  double y3 = g_flattenList[2].y();
  double z3 = find_average_value_around(from, x3, y3, g_radius);

  //------------------------------------------------------------
  // Calculate the plane equation passing through the three points.
  // solve for dx, dy in:
  //  z3 - z1 = dx * ( x3 - x1 ) + dy * ( y3 - y1 )
  //  z2 - z1 = dx * ( x2 - x1 ) + dy * ( y2 - y1 )
  // as so:
  //  dy = (z2 - z1 + (z1 - z3) * k) / (y2 - y1 + (y1 - y3) * k)
  //  dx = (z3 - z1 - dy * (y3 - y1)) / (x3 - x1)
  //  where k = ( x2 - x1 ) / ( x3 - x1 )

  // Are the three points co-linear?  If so, give up.
  if( (y2-y1) * (x3-x1) + (y1-y3) * (x2-x1) == 0 ) {
    return;
  }
 
  // Are pt1 and pt3 co-linear in x?
  // If so, swap pt 2 and pt 3, or division will go bad.
  if( x3 == x1 ) {
      double temp = x3;  x3 = x2;  x2 = temp;
      temp = y3;  y3 = y2;  y2 = temp;
  }

  double k = ( x2 - x1 ) / ( x3 - x1 );
  double dy = (z2 - z1 + (z1 - z3) * k) / (y2 - y1 + (y1 - y3) * k);
  double dx = (z3 - z1 - dy * (y3 - y1)) / (x3 - x1);
  double offset = dx * x1 + dy * y1;

  //------------------------------------------------------------
  // Subtract the plane equation from the original data values,
  // resulting in a new data set that has those three points at
  // zero.  This will fill in the valus for the new plane.
  int x,y;
  for (x = 0; x < from->numX(); x++) {
    for (y = 0; y < from->numY(); y++) {
      double value = from->scaledValue(x,y) - (offset + dx*x + dy*y );
      to->setValue(x,y, value);
    }
  }
}

//--------------------------------------------------------------------------
// Calculate a plane from the flat plane that stores a mask telling
// whether the result is above threshold or not.  The mask holds 0 for
// pixels below threshold and 1 for pixels above threshold.
// The routine needs to be robust to being called before its inputs
// are defined.
static void calculate_height_threshold(const string volume_plane_name, const string threshold_plane_name)
{
  //------------------------------------------------------------
  // If we don't have an input plane, bail.
  if (g_current_data == NULL) { return; }
  BCPlane *from = g_current_data->getPlaneByName(volume_plane_name);
  if (from == NULL) { return; }

  //------------------------------------------------------------
  // See if the height_threshold plane exists already. If not, then
  // create it.
  BCPlane *to = g_current_data->getPlaneByName(threshold_plane_name);
  if (to == NULL) {
    to = g_current_data->addPlaneCopy(from);
    if (to == NULL) {
      fprintf(stderr,"calculate_height_threshold(): Could not create new plane\n");
      exit(-1);
    }
    to->rename(threshold_plane_name);
  }

  double  offset = from->minValue();
  double  scale = from->maxValue() - from->minValue();
  double  scaled_threshold = offset + (g_height_threshold * scale);

  //------------------------------------------------------------
  // Test each value against the scaled threshold.
  int x,y;
  for (x = 0; x < from->numX(); x++) {
    for (y = 0; y < from->numY(); y++) {
      if (from->scaledValue(x,y) < scaled_threshold) {
        to->setValue(x,y, 0);
      } else {
        to->setValue(x,y, 1);
      }
    }
  }
}

//--------------------------------------------------------------------------
// Helper routine that fills all pixels that are connected to
// the indicated pixel according to the "mask" plane (which holds 0 for
// pixels that are off and 1 for pixels that are on).  The pixels are
// set in the "fill" plane.  Recursion failed due to stack overflow when
// the whole image was one region, so the new algorithm keeps a list of
// "border positions" and repeatedly goes through that list until it is
// empty.
static void flood_connected_component(const BCPlane *mask, BCPlane *fill, int x, int y, double value)
{
  list<Position>  boundary;

  //------------------------------------------------------------
  // Fill in the first pixel as requested, and add it to the
  // boundary list.
  fill->setValue(x,y, value);
  boundary.push_front(Position(x,y));

  //------------------------------------------------------------
  // Go looking for ways to expand the boundary so long as the
  // boundary is not empty.  Once each boundary pixel has been
  // completely handled, delete it from the boundary list.  
  while (!boundary.empty()) {
    list<Position>::iterator  i;

    // We add new boundary pixels to the beginning of the list so that
    // we won't see them while traversing the list: each time through,
    // we only consider the positions that were on the boundary when
    // we started through that time.
    i = boundary.begin();
    while (i != boundary.end()) {

      // Move to the location of the pixel and check around it.
      x = i->x();
      y = i->y();

      // If any of the four neighbors are valid new boundary elements,
      // mark them and insert them into the boundary list.
      if ( x-1 >= 0 ) {
        if (mask->scaledValue(x-1,y) && (fill->scaledValue(x-1,y) != value) ) {
          fill->setValue(x-1,y, value);
          boundary.push_front(Position(x-1,y));
        }
      }
      if ( x+1 < mask->numX() ) {
        if (mask->scaledValue(x+1,y) && (fill->scaledValue(x+1,y) != value) ) {
          fill->setValue(x+1,y, value);
          boundary.push_front(Position(x+1,y));
        }
      }
      if ( y-1 >= 0 ) {
        if (mask->scaledValue(x,y-1) && (fill->scaledValue(x,y-1) != value) ) {
          fill->setValue(x,y-1, value);
          boundary.push_front(Position(x,y-1));
        }
      }
      if ( y+1 < mask->numY() ) {
        if (mask->scaledValue(x,y+1) && (fill->scaledValue(x,y+1) != value) ) {
          fill->setValue(x,y+1, value);
          boundary.push_front(Position(x,y+1));
        }
      }
      
      // Get rid of this entry, stepping to the next one before doing so.
      list<Position>::iterator j = i;
      ++i;
      boundary.erase(j);
    }
  }
}


//--------------------------------------------------------------------------
// Go through and fill in all of the elements whose value is val1 with
// val2.
static void replace_component(BCPlane *fill, double val1, double val2)
{
  int x,y;
  for (x = 0; x < fill->numX(); x++) {
    for (y = 0; y < fill->numY(); y++) {
      if (fill->scaledValue(x,y) == val1) {
        fill->setValue(x,y, val2);
      }
    }
  }
}

//--------------------------------------------------------------------------
// Dilate the region with the specified value by one pixel.  Increase both
// the mask and the index.  If this makes two regions overlap, set the
// region values to -2 (unless one is the boundary, in which case it
// stays -1).
// For speed reasons, the calling routine must pass in a boundary list
// that contains at least the points surrounding the specified region.

static void dilate(BCPlane *mask, BCPlane *fill, const vector<Position> superboundary, double value)
{
  unsigned i;
  int x,y;
  vector<Position>  boundary;

  // Allocate enough space to fill in a bunch of boundary entries
  // to avoid taking a long time in memory allocation.
  boundary.reserve(superboundary.size() / 10 + 1);

  //------------------------------------------------------------
  // Before we do anything, make a list of points that are right
  // next to the specified region.  We use this to speed up all
  // of the following steps and also to ensure that we only do
  // operations on points that were on the original boundary.
  for (i = 0; i < superboundary.size(); i++) {
    x = superboundary[i].x();
    y = superboundary[i].y();

    // Skip points already in the region.
    if ( fill->scaledValue(x,y) == value) {
      continue;
    }

    // Locate points with neighbors in the region.
    if ( x-1 >= 0 ) {
      if (fill->scaledValue(x-1,y) == value) {
        boundary.push_back(Position(x,y));
        continue;
      }
    }
    if ( x+1 < mask->numX() ) {
      if (fill->scaledValue(x+1,y) == value) {
        boundary.push_back(Position(x,y));
        continue;
      }
    }
    if ( y-1 >= 0 ) {
      if (fill->scaledValue(x,y-1) == value) {
        boundary.push_back(Position(x,y));
        continue;
      }
    }
    if ( y+1 < mask->numY() ) {
      if (fill->scaledValue(x,y+1) == value) {
        boundary.push_back(Position(x,y));
        continue;
      }
    }
  }

  //------------------------------------------------------------
  // First, check points right next to the specified region to
  // see if they are marked with a different value; if so, we need
  // to clobber all of them to -2 unless they are the boundary.
  // Go through the whole list and clobber all of the other ones,
  // marking the passed-in one for deletion after the pass is complete.
  // This order is to ensure that we get all potentially-overlapping
  // regions before marking our target region invalid.
  // XXX If we later allow regions to be merged, make sure to both
  // fill with our value and insert all of the points into the active
  // boundary list, or do whatever is correct to keep the boundary
  // list up to date.
  bool overlapped = false;
  for (i = 0; i < boundary.size(); i++) {
    x = boundary[i].x();
    y = boundary[i].y();
    double val;
    if ( (val = fill->scaledValue(x,y)) != -3) {
      overlapped = true;
      if (val >= 0) {
        replace_component(fill, val, -2);
      }
    }
  }
  if (overlapped) {
    replace_component(fill, value, -2);
    return;
  }

  //------------------------------------------------------------
  // If no overlap, go through and make all empty (value -3) neighbor cells
  // to the specified region as belonging to that region and expand the
  // mask for those regions.
  for (i = 0; i < boundary.size(); i++) {
    x = boundary[i].x();
    y = boundary[i].y();
    mask->setValue(x,y, 1);
    fill->setValue(x,y, value);
  }
}


//--------------------------------------------------------------------------
// Calculate a plane from the height_threshold plane that stores a mask telling
// which connected component each pixel belongs to.  The mask holds -3 for
// pixels that are below the threshold, -2 for pixels in components that are too
// close to each other, -1 for pixes that are on the border or connected to
// the border, and 0-N for actual components.  The mask holds -4 for regions
// whose volume is outside the range of desired volumes.
// The routine needs to be robust to being called before its inputs
// are defined.
static void calculate_connected_components(const string threshold_plane_name, const string component_string_name)
{
  //------------------------------------------------------------
  // If we don't have an input plane, bail.
  if (g_current_data == NULL) {
    return;
  }
  BCPlane *from = g_current_data->getPlaneByName(threshold_plane_name);
  if (from == NULL) {
    return;
  }

  //------------------------------------------------------------
  // See if the connected_components plane exists already. If not, then
  // create it.
  BCPlane *to = g_current_data->getPlaneByName(component_string_name);
  if (to == NULL) {
    to = g_current_data->addPlaneCopy(from);
    if (to == NULL) {
      fprintf(stderr,"calculate_connected_components(): Could not create new plane\n");
      exit(-1);
    }
    to->rename(component_string_name);
  }

  double  offset = from->minValue();
  double  scale = from->maxValue() - from->minValue();
  double  scaled_threshold = offset + (g_height_threshold * scale);

  //------------------------------------------------------------
  // Fill the whole array with -3, which marks pixels as being masked
  // by the height_threshold.
  int x,y;
  for (x = 0; x < from->numX(); x++) {
    for (y = 0; y < from->numY(); y++) {
      to->setValue(x,y, -3);
    }
  }

  //------------------------------------------------------------
  // Mark the border pixels as belonging to the border component,
  // and also mark all pixels which are above threshold and connected
  // to a border pixel.

  for (x = 0; x < from->numX(); x++) {
    flood_connected_component(from, to, x,0, -1);
    flood_connected_component(from, to, x,from->numY()-1, -1);
  }
  for (y = 0; y < from->numY(); y++) {
    flood_connected_component(from, to, 0,y, -1);
    flood_connected_component(from, to, from->numX()-1,y, -1);
  }

  //------------------------------------------------------------
  // Pass through and mark unlabeled pixels into successive
  // components starting with index 0.
  int index = -1;
  for (x = 0; x < from->numX(); x++) {
    for (y = 0; y < from->numY(); y++) {
      if (from->scaledValue(x,y) && (to->scaledValue(x,y) == -3)) {
        index++;
        flood_connected_component(from, to, x,y, index);
      }
    }
  }

  //------------------------------------------------------------
  // Dilate the regions other than the boundary region by the amount
  // of the global dilation setting.  During each dilation step for
  // each region, expand the mask around the existing regions but
  // check to see if the dilation causes an overlap of two regions.
  // If it does, set the regions each to -2 (unless one is the boundary,
  // which is left alone.
  int i, j;
  for (j = 0; j < g_dilation; j++) {

    // First we get a list of all candidate boundary points:
    // those whose neighbors differ in value from itself.  This
    // is a superset of all points, and is passed into the dilation
    // code to speed things up.
    vector<Position>  boundary;
    for (x = 0; x < to->numX(); x++) {
      for (y = 0; y < to->numY(); y++) {

        // Find out what our value is.
        double value = to->scaledValue(x,y);

        // Locate points with neighbors in the region.
        if ( x-1 >= 0 ) {
          if (to->scaledValue(x-1,y) != value) {
            boundary.push_back(Position(x,y));
            continue;
          }
        }
        if ( x+1 < to->numX() ) {
          if (to->scaledValue(x+1,y) != value) {
            boundary.push_back(Position(x,y));
            continue;
          }
        }
        if ( y-1 >= 0 ) {
          if (to->scaledValue(x,y-1) != value) {
            boundary.push_back(Position(x,y));
            continue;
          }
        }
        if ( y+1 < to->numY() ) {
          if (to->scaledValue(x,y+1) != value) {
            boundary.push_back(Position(x,y));
            continue;
          }
        }
      }
    }

    // Now dilate each region.
    for (i = 0; i <= index; i++) {
      dilate(from, to, boundary, i);
    }

    // Clear the boundary for next time.
    boundary.clear();
  }
}

//--------------------------------------------------------------------------
// Calculate a list of volumes for each component.  Note that the component
// indices may have gaps (due to regions being destroyed during dilation).
// The list of components is not sorted.
static void calculate_component_volumes(const string volume_plane_name, const string component_plane_name, vector<Component_volume> &volumes)
{
  //------------------------------------------------------------
  // Nothing there yet!
  volumes.clear();

  //------------------------------------------------------------
  // See if the connected_components plane exists already. Also make sure
  // the flattened plane exists, which we will use to determine the volume.
  if (g_current_data == NULL) {
    return;
  }
  BCPlane *flat = g_current_data->getPlaneByName(volume_plane_name);
  BCPlane *comps = g_current_data->getPlaneByName(component_plane_name);
  if ( (flat == NULL) || (comps == NULL) ) {
    return;
  }

  //------------------------------------------------------------
  // Add up the volumes of each pixel for each component found in
  // the image.  Keep track of this in the global list of volumes.
  // Insert entries as needed.

  int x,y;
  for (x = 0; x < comps->numX(); x++) {
    for (y = 0; y < comps->numY(); y++) {

      // See if the id for this point is for a surviving cluster.
      double id = comps->scaledValue(x,y);
      if (id >= 0) {
        // Find the relevant ID in the list or add it to the end.
        unsigned i;
        int which = -1;
        for (i = 0; i < volumes.size(); i++) {
          if (volumes[i].id() == id) {
            which = i;
            break;
          }
        }
        if (which == -1) {
          volumes.push_back(Component_volume(id,0));
          which = volumes.size() - 1;
        }

        // Compute the volume of the pixel and add it to the
        // component volume.
        double volume = flat->scaledValue(x,y) / (flat->derangeX() * flat->derangeY());
        volumes[which].set_volume(g_componentVolumes[which].volume() + volume);
      }
    }
  }
}

static void calculate_thresholded_volumes(const string in_plane_name, const vector<Component_volume> &in_volumes,
                                          const double norm_min, double norm_max,
                                          const string out_plane_name, vector<Component_volume> &out_volumes)
{
  //------------------------------------------------------------
  // Nothing there yet!
  out_volumes.clear();

  //------------------------------------------------------------
  // Make sure we have a list to go through
  if (in_volumes.size() == 0) {
    return;
  }

  //------------------------------------------------------------
  // See if the connected_components plane exists.  If not, we're
  // done.
  if (g_current_data == NULL) {
    return;
  }
  BCPlane *comps = g_current_data->getPlaneByName(in_plane_name);
  if (comps == NULL) {
    return;
  }

  //------------------------------------------------------------
  // See if the output plane exists already. If not, then
  // create it.  If so, delete it and reconstruct it so that
  // we end up with a copy of the original plane in there.
  BCPlane *to = g_current_data->getPlaneByName(out_plane_name);
  if (to == NULL) {
    to = g_current_data->addPlaneCopy(comps);
    if (to == NULL) {
      fprintf(stderr,"calculate_flatten(): Could not create new plane\n");
      exit(-1);
    }
  } else {
    g_current_data->removePlane(out_plane_name);
    to = g_current_data->addPlaneCopy(comps);
    if (to == NULL) {
      fprintf(stderr,"calculate_flatten(): Could not create new plane\n");
      exit(-1);
    }
  }
  to->rename(out_plane_name);

  //------------------------------------------------------------
  // Go through the list of components and find the min and max
  // volumes.  Use this to figure a scale and offset for the
  // thresholds.  Use that to figure out a scaled min and max
  // value that corresponds to the sliders.
  double min = in_volumes[0].volume();
  double max = min;
  unsigned i;
  for (i = 0; i < in_volumes.size(); i++) {
    double volume = in_volumes[i].volume();
    if (volume < min) { min = volume; }
    if (volume > max) { max = volume; }
  }
  double offset = min;
  double range = (max - min);
  double scaled_min = offset + norm_min * range;
  double scaled_max = offset + norm_max * range;

  //------------------------------------------------------------
  // Go through the list of components again.  For each one whose
  // volume is within the range, insert it into the list of output
  // volumes.  For ones outside the range, change all entries in
  // the output plane to -4 to indicate they were out of the
  // volume range.
  for (i = 0; i < in_volumes.size(); i++) {
    double volume = in_volumes[i].volume();
    if ( (volume < scaled_min) || (volume > scaled_max) ) {
      replace_component(to, in_volumes[i].id(), -4);
    } else {
      out_volumes.push_back(in_volumes[i]);
    }
  }
}

//--------------------------------------------------------------------------
// Helper routine to get the Y coordinate right when going between camera
// space and openGL space.
static double	flip_y(double y)
{
  return g_nRowsData - 1 - y;
}


//-------------------------------------------------------------------
// See if there is a valid component at the stated location. 
// If so, the user is asked for the name
// of a template file (UNCA format) to save the data in, then
// the appropriate component is splatted into a new plane and
// a subset of the plane extracted that includes all nonzero
// elements.  This is saved into the template file, which is
// then re-loaded immediately into the template list and made
// the current template.

bool  make_template_from_component_at(unsigned x, unsigned y)
{
  //---------------------------------------------------------
  // Find the plane that holds the connected components by
  // number and look up the number at the selected pixel.
  int comp_id;
  BCPlane *flat = g_current_data->getPlaneByName("flat");
  BCPlane *thresholded_component = g_current_data->getPlaneByName("thresholded_components");
  if ( (flat == NULL) || (thresholded_component == NULL) ) {
    fprintf(stderr, "make_template_from_component_at(): Can't find component plane\n");
    return false;
  }
  comp_id = thresholded_component->scaledValue(x,y);

  // If we don't have a valid (non-negative) component, bail.
  if (comp_id < 0) {
    return false;
  }

  //---------------------------------------------------------
  // Ask the user for the name of a file into which to save
  // the component and wait until they respond or cancel.
  // If they cancel, bail (but reset quit so we don't exit.
  if (Tcl_Eval(g_tk_control_interp, "ask_user_for_save_template_filename") != TCL_OK) {
    fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", "ask_user_for_save_template__filename", g_tk_control_interp->result);
    cleanup();
    exit(-1);
  }
  // Wait until we get a filename or are asked to quit
  do {
    // This pushes changes in the C variables over to Tcl.
    while (Tk_DoOneEvent(TK_DONT_WAIT)) {};
    if (Tclvar_mainloop()) {
      fprintf(stderr,"Tclvar Mainloop failed\n");
    }
  } while ( (g_save_template_filename == NULL) && !g_quit);
  if (g_quit) {
    g_quit = 0;
    return false;
  }

  //---------------------------------------------------------
  // Make a new plane to splat the component into; clear the
  // elements of the plane and then
  // copy the component's nonzero elements in.  Find the axis-aligned
  // bounding box of its nonzero elements.
  BCPlane *splat = g_current_data->getPlaneByName("temp_splat");
  if (splat == NULL) {
    splat = g_current_data->addPlaneCopy(g_current_data->getPlaneByName("flat"));
    if (splat == NULL) {
      fprintf(stderr, "make_template_from_component_at(): Can't create splat plane\n");
      cleanup();
    }
    splat->rename("temp_splat");
  }
  int tx,ty;
  int maxx = 0, maxy = 0;
  int minx = splat->numX()-1, miny = splat->numY() - 1;
  for (tx = 0; tx < splat->numX(); tx++) {
    for (ty = 0; ty < splat->numY(); ty++) {
      if (comp_id == thresholded_component->scaledValue(tx,ty)) {
        splat->setValue(tx,ty, flat->scaledValue(tx,ty));
        if (tx < minx) { minx = tx; }
        if (tx > maxx) { maxx = tx; }
        if (ty < miny) { miny = ty; }
        if (ty > maxy) { maxy = ty; }
      } else {
        splat->setValue(tx,ty, 0.0);
      }
    }
  }

  //---------------------------------------------------------
  // Save the bounding-box subset of the plane into a UNCA
  // file with the name the user gave us.  First create a new
  // grid with a smaller number of elements and a different
  // size (proper subset) of the existing grid, then copy
  // the bounding box into it.  Finally, save the plane from
  // the new grid into a UNCA file with the specified name.
  int num_x = maxx-minx+1;
  int num_y = maxy-miny+1;
  double min_x, max_x, min_y, max_y;
  g_current_data->gridToWorld(minx, maxx, min_x, max_x);
  g_current_data->gridToWorld(miny, maxy, min_y, max_y);
  BCGrid  save_grid( num_x,num_y, min_x, max_x, min_y, max_y, READ_FILE );
  BCPlane *save_plane = save_grid.addNewPlane("data", flat->units()->c_str(), 0);
  if (save_plane == NULL) {
    fprintf(stderr, "make_template_from_component_at(): Can't create save plane\n");
    return false;
  }
  for (tx = 0; tx < num_x; tx++) {
    for (ty = 0; ty < num_y; ty++) {
      save_plane->setValue(tx,ty, splat->scaledValue(minx + tx, miny + ty));
    }
  }

  // Save the plane into a UNCA file.
  FILE *outfile = fopen((char*)(g_save_template_filename), "w");
  if (outfile == NULL) {
    fprintf(stderr, "make_template_from_component_at(): Can't open output file\n");
    save_grid.removePlane(save_plane->name()->c_str());
    return false;
  }
  if (save_grid.writeUNCAFile(outfile, save_plane) != 0) {
    fprintf(stderr, "make_template_from_component_at(): Can't save plane\n");
    save_grid.removePlane(save_plane->name()->c_str());
    return false;
  }
  fclose(outfile);
  save_grid.removePlane(save_plane->name()->c_str());

  //---------------------------------------------------------
  // Load the named file as a new template and make it the
  // current template.
  template_filename_changed((char*)(g_save_template_filename), NULL);
  return true;
}

//-------------------------------------------------------------------
// Clicking the mouse in the second window selects one of the
// regions to use as a template.

void secondMouseCallbackForGLUT(int button, int state, int x, int y)
{
  // Flip the y value to make it match the image coordinates.
  y = flip_y(y);

  switch(button) {
    // The right and middle mouse buttons do nothing.
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
        // Nothing on press
      } else {
        // Nothing on release
      }
      break;

    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) {
        // Nothing on press
      } else {
        // Nothing on release
      }
      break;

    // The left mouse button selects a component if the mouse is
    // over one and then starts the action.
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {
        make_template_from_component_at(x,y);
      } else {
        // Nothing on release
      }
      break;
  }
}

//----------------------------------------------------------------------------
// I got this function from Gary Bishop.  He got it from a NASA tech
// report.  It does a non-causal low-pass filter of a vector in several
// passes, returning the filtered vector.
// The vector (of length 'len') is modified in-place.  The filter cutoff
// is specified in samples/cycle as 'freq'.  The number of passes to make
// with the algorithm are specified in 'passes'.
// Note: Original function specified parameters as dt (distance between
//       sample points [dist/sample]) and fc (frequency of cutoff in Hz
//       [cycles/dist]).  This version calculates those parameters from
//       'freq [sample/cycle]' and runs the original code.

static void low_pass_filter_line(double *vec, int len, double freq, int passes)
{
    double dt = 1.0;			// Assume unit sample spacing
    double fc = dt/freq;

    double wc = 2 * M_PI * fc;
    double tau = sqrt(pow(2.0, 1.0/(2*passes)) - 1) / wc;

    double k1 = exp(-dt/tau);
    double k2 = -(k1 - tau/dt * (1 - k1));
    double k3 = 1 - tau/dt * (1 - k1);

    for(int p=0; p<passes; p++) {
        double g0 = (vec[0]+vec[1])/2;
	int i;
        for(i=0; i<len-1; i++) {
            double g1 = k1*g0 + k2*vec[i] + k3*vec[i+1];
            vec[i] = g0;
            g0 = g1;
        }
        vec[len-1] = g0;
        g0 = (vec[len-1]+vec[len-2])/2;
        for(i=len-1; i>0; i--) {
            double g1 = k1*g0 + k2*vec[i] + k3*vec[i-1];
            vec[i] = g0;
            g0 = g1;
        }
        vec[0] = g0;
    }
}

/****************************************************************************
        Gary Bishop's "secret sauce" flattening algorithm.
	Reads a block of X by Y floats into an array.  This is an image that
 has been scanned and has line artifacts in it.  The code attempts to remove
 the high-frequency components of this scanline noise from the image without
 removing slope or other information from the actual data.
	The algorithm is to assume that the image M(x,y) is composed of two
 signals: S(x,y) is the surface data and N(y) is the noise (assumed to depend
 only on y, not on x).  We apply a low-pass filter to M(x,y) in only the y
 direction to produce L(x,y) [low-pass version].  We form the high-pass
 array as: H(x,y) = M(x,y) - L(x,y).  This gives us a high-pass (in the y
 direction) version of the original function.  The assumption is that the
 high-frequency components of S(x,y) are zero-mean.  Thus, we find the
 function N(y) as the average over x of H(x,y): N(y) = SUMx(H(x,y))/nx.
 We produce the output O(x,y) = M(x,y) - N(y).
	The parameters are nx, ny (number of samples in x and y), f (cutoff
 frequency for low-pass in pixels) and p (number of passes of the low-pass
 filter).
	The code exits with no output and status -1 if it can't read all
 of its input.  It exits with status -1 if it can't write all of its output.
 If things work, it exits with code 0.
        The from and to plane can be the same.
 ****************************************************************************/

static bool perform_deline_filter_old(const BCPlane *from, BCPlane *to)
{
  // Check our parameters.
  if ( (from == NULL) || (to == NULL) ) {
    fprintf(stderr,"perform_deline_filter(): Called with NULL pointer!\n");
    return false;
  }
  if ( (from->numX() != to->numX()) || (from->numY() != to->numY()) ) {
    fprintf(stderr,"perform_deline_filter(): Called with planes of different sizes!\n");
    return false;
  }

  int nx = from->numX();
  int ny = from->numY();
  float freq = 64;   //< Samples/cycle cutoff frequency
  int passes = 10;   //< How many passes of the above function to iterate
  int x, y;

  // Allocate the needed arrays
  float	*M, *N, *L, *H;
  double *T;
  if ( ((M = new float[nx*ny]) == NULL) ||
       ((L = new float[nx*ny]) == NULL) ||
       ((H = new float[nx*ny]) == NULL) ||
       ((N = new float[ny]) == NULL) ||
       ((T = new double[ny]) == NULL) ) {
	  fprintf(stderr,"perform_deline_filter(): Out of memory\n");
	  return false;
  }

  // Read the data into M(x,y)
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      M[x + y*nx] = from->scaledValue(x,y);
    }
  }

  // Apply the low-pass filter to produce L(x,y) one line at a time
  // Uses T(y) as the temporary vector to hold each for filtering
  for (x = 0; x < nx; x++) {
	  // Get a column of M into T
	  for (y = 0; y < ny; y++) { T[y] = M[x+y*nx]; }

	  // Filter the column
	  low_pass_filter_line(T, ny, freq, passes);

	  // Move the column into L
	  for (y = 0; y < ny; y++) { L[x+y*nx] = T[y]; }
  }

  // Subtract to get H(x,y)
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
	  H[x+y*nx] = M[x+y*nx] - L[x+y*nx];
    }
  }

  // Average across rows to find N(y)
  for (y = 0; y < ny; y++) {
	  N[y] = 0.0;
	  for (x = 0; x < nx; x++) { N[y] += H[x+y*nx]; };
	  N[y] /= nx;
  }

  // Subtract N(y) from M(x,y) to produce the output image
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
	  M[x+y*nx] -= N[y];
    }
  }

  // Copy the result into the output plane
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
       to->setValue(x,y,M[x + y*nx]);
    }
  }

  // Free the memory we allocated
  delete [] T;
  delete [] N;
  delete [] H;
  delete [] L;
  delete [] M;
  return true;
}

/****************************************************************************
        Russ Taylor's "secret sauce" flattening algorithm.
	Reads a block of X by Y floats into an array.  This is an image that
 has been scanned and has line artifacts in it.  The code attempts to remove
 the high-frequency components of this scanline noise from the image without
 removing slope or other information from the actual data.
        It operates by first low-pass filtering each line of the image
 (blurring in X).  For each line after the first, it finds the median step
 between in and the previous line.  This step is used as an offset for all
 values in that line (and all later lines).  The original values from the line
 (not the blurred values) are all shifted by the opposite of this amount.  This
 removes some noise and then offsets so that the median step from one line to
 the next is zero.  Some steps will be negative, and some positive; this median
 value assumes that most of the line is substrate, so that we're still in the
 substrate-step range when we get to the median.
	The parameters are nx, ny (number of samples in x and y), f (cutoff
 frequency for low-pass in pixels) and p (number of passes of the low-pass
 filter).
	The code exits with no output and status -1 if it can't read all
 of its input.  It exits with status -1 if it can't write all of its output.
 If things work, it exits with code 0.
        The from and to plane can be the same.
 ****************************************************************************/

static bool perform_deline_filter(const BCPlane *from, BCPlane *to)
{
  // Check our parameters.
  if ( (from == NULL) || (to == NULL) ) {
    fprintf(stderr,"perform_deline_filter(): Called with NULL pointer!\n");
    return false;
  }
  if ( (from->numX() != to->numX()) || (from->numY() != to->numY()) ) {
    fprintf(stderr,"perform_deline_filter(): Called with planes of different sizes!\n");
    return false;
  }

  int nx = from->numX();
  int ny = from->numY();
  float freq = 10;   //< Samples/cycle cutoff frequency
  int passes = 5;   //< How many passes of the above function to iterate
  int x, y;

  // Allocate two lines, one for the "previous" and one for the "current"
  // line.  These are blurred and then used to find the median step.
  double *previous = new double[nx];
  double *current = new double[nx];
  if ( !current || !previous ) {
    fprintf(stderr,"perform_deline_filter(): Out of memory\n");
    return false;
  }

  // First, just copy the zeroeth line into the output, and put its
  // blurred contents into previous, so that the following loop code
  // will work correctly in the first iteration.
  for (x = 0; x < nx; x++) {
    double val = from->scaledValue(x,0);
    previous[x] = val;
    to->setValue(x,0, val);
  }
  low_pass_filter_line(previous, nx, freq, passes);

  // Then do the processing on each subsequent pass, doing the
  // median offset after blurring each line.  Remember to sum the
  // total median offset since the beginning, because all following
  // lines also need to be offset
  vector<double>  steps;
  double          median;
  double          sum_of_medians = 0;
  for (y = 1; y < ny; y++) {

    // Copy and blur the values for this line
    for (x = 0; x < nx; x++) {
      current[x] = from->scaledValue(x,y);
    }
    low_pass_filter_line(current, nx, freq, passes);

    // Find the median step size.
    steps.clear();
    for (x = 0; x < nx; x++) {
      steps.push_back(current[x] - previous[x]);
    }
    sort(steps.begin(), steps.end());
    median = steps[nx/2];
    sum_of_medians += median;

    // Copy the values, offsetting by the negative of the median
    // step size.
    for (x = 0; x < nx; x++) {
      to->setValue(x,y, from->scaledValue(x,y) - sum_of_medians);
    }

    // Copy the current into the previous for next time
    for (x = 0; x < nx; x++) {
      previous[x] = current[x];
    }
  }

  // Free the memory we allocated
  delete [] previous;
  delete [] current;
  return true;
}

//--------------------------------------------------------------------------
// Performs a 3x3 median filter on the from plane and stores the result
// into the to plane.  These planes can be the same.

static bool calculate_median_filter(const string from_name, const string to_name)
{
  //------------------------------------------------------------
  // If we don't have an input plane, bail.
  if (g_current_data == NULL) { return false; }
  BCPlane *from = g_current_data->getPlaneByName(from_name);
  if (from == NULL) { return false; }

  //------------------------------------------------------------
  // See if the flattened plane exists already. If not, then
  // create it.
  BCPlane *to = g_current_data->getPlaneByName(to_name);
  if (to == NULL) {
    to = g_current_data->addPlaneCopy(from);
    if (to == NULL) {
      fprintf(stderr,"calculate_median_filter(): Could not create new plane\n");
      exit(-1);
    }
    to->rename(to_name);
  }

  if ( (from->numX() != to->numX()) || (from->numY() != to->numY()) ) {
    fprintf(stderr,"calculate_median_filter(): Called with planes of different sizes!\n");
    return false;
  }

  int nx = from->numX();
  int ny = from->numY();
  int x, y;

  // Allocate the needed arrays
  double *M;
  if ( (M = new double[nx*ny]) == NULL) {
	  fprintf(stderr,"calculate_median_filter(): Out of memory\n");
	  return false;
  }

  // Read the data into M(x,y)
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      M[x + y*nx] = from->scaledValue(x,y);
    }
  }

  //---------------------------------------
  // Median-filter M into the output plane.

  // Copy pixels on the boundary of the image.
  for (x = 0; x < nx; x++) {
    to->setValue(x,0,M[x + 0*nx]);
    to->setValue(x,ny-1,M[x + (ny-1)*nx]);
  }
  for (y = 0; y < ny; y++) {
    to->setValue(0,y,M[0 + y*nx]);
    to->setValue(nx-1,y,M[(nx-1) + y*nx]);
  }

  // Filter the pixels not on the boundary of the image.
  vector<double>  values;
  for (x = 1; x < nx-1; x++) {
    for (y = 1; y < ny-1; y++) {

      // Get the nine surrounding values into a linear array.
      int i,j;
      values.clear();
      for (i = -1; i <= 1; i++) {
        for (j = -1; j <= 1; j++) {
          values.push_back(M[(x+i) + (y+j)*nx]);
        }
      }

      // Sort the nine values.
      sort(values.begin(), values.end());

      // Write the middle (median) value into the result.
      to->setValue(x,y, values[4]);
    }
  }

  delete [] M;
  return true;
}

//--------------------------------------------------------------------------
// Helper routine to drag the nearest tracker (if any exist) to the
// specified coordinate.  Returns true on success, false on failure
// (no tracker).
static  bool  drag_nearest_tracker_to(int x, int y)
{
  // Find the nearest point, if there are any points
  int which = -1;
  double dist2 = 1e50;
  unsigned i;
  for (i = 0; i < g_flattenList.size(); i++) {
    Position p = g_flattenList[i];
    double tempd2 = pow(static_cast<double>(x - p.x()), 2) + pow(static_cast<double>(y - p.y()), 2);
    if (tempd2 < dist2) {
      dist2 = tempd2;
      which = i;
    }
  }
  if (which == -1) { return false; }

  // Drag the nearest point to the specified location
  g_flattenList[which].set_x(x);
  g_flattenList[which].set_y(y);
  return true;
}

void dataMouseCallbackForGLUT(int button, int state, int x, int y)
{
  // Record where the button was pressed for use in the motion
  // callback, flipping the Y axis to make the coordinates match
  // image coordinates.
  g_mousePressX = x;
  g_mousePressY = y = flip_y(y);

  switch(button) {
    // The right mouse button will commence dragging the nearest
    // point (if there are any active).
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
	g_whichDragAction = 1;
	drag_nearest_tracker_to(x,y);
      } else {
	// Nothing to do at release.
        g_whichDragAction = 0;
      }
      break;

    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) {
	g_whichDragAction = 0;
      } else {
        g_whichDragAction = 0;
      }
      break;

    // The left mouse button will place a new point (if there
    // are less than three) and commence dragging that point.
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {
        if (g_flattenList.size() < 3) {
          g_flattenList.push_back(Position(x,y));
        }
        drag_nearest_tracker_to(x,y);
	g_whichDragAction = 1;
      } else {
        g_whichDragAction = 0;
      }
      break;
  }
  float_causes_redisplay(0.0, NULL);
}

void dataMotionCallbackForGLUT(int x, int y) {

  // Make mouse coordinates match image coordinates.
  y = flip_y(y);

  // Make sure we stay within the image; don't let them
  // drag it off into space.
  if (x < 0) { x = 0; }
  if (x >= static_cast<int>(g_nColsData)) { x = g_nColsData - 1; }
  if (y < 0) { y = 0; }
  if (y >= static_cast<int>(g_nRowsData)) { y = g_nRowsData - 1; }

  switch (g_whichDragAction) {

  case 0: //< Do nothing on drag.
    break;

  // Pull the closest existing tracker
  // to the location where the mouse button was pressed, and
  // keep pulling it around if the mouse is moved while this
  // button is held down.
  case 1:
    drag_nearest_tracker_to(x,y);
    break;

  default:
    fprintf(stderr,"Internal Error: Unknown drag action (%d)\n", g_whichDragAction);
  }
  float_causes_redisplay(0.0, NULL);
}

void Usage (const char * s)
{
  fprintf(stderr,"Usage: %s [-v] [-maxwindow size] [datafile]\n",s);
  fprintf(stderr,"     -v: Verbose mode, print out lots of info along the way\n");
  fprintf(stderr,"     datafile: Names of zero or one data files to load\n");
  exit(-1);
}

int main (int argc, char * argv[])
{
  //------------------------------------------------------------------------
  // Parse the command line
  int	realparams = 0;
  int i;
  i = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-v") == 0) {
	verbose = true;
    } else if (argv[i][0] == '-') {	// Unknown flag
	Usage(argv[0]);
    } else switch (realparams) {		// Non-flag parameters
      case 0: // Don't forget to break if we have real arguments added later
      case 1:
      default:
	// Load as many data files as are listed, filling them into a list and
	// checking if they are all the same resolution.  Set the rows and
	// columns to that resolution.
	data_filename_changed(argv[i], NULL);
	realparams++;
    }
    i++;
  }
  if (realparams > 1) {
    Usage(argv[0]);
  }

  //------------------------------------------------------------------
  // Set up exit handler to make sure we clean things up no matter
  // how we are quit.  We hope that we exit in a good way and so
  // cleanup() gets called, but if not then we do a dirty exit.
  atexit(cleanup);

  //------------------------------------------------------------------
  // Generic Tcl startup.  Getting and interpreter and mainwindow.

  char		command[256];
  Tk_Window       tk_control_window;
  g_tk_control_interp = Tcl_CreateInterp();

  /* Start a Tcl interpreter */
  if (Tcl_Init(g_tk_control_interp) == TCL_ERROR) {
          fprintf(stderr,
                  "Tcl_Init failed: %s\n",g_tk_control_interp->result);
          return(-1);
  }

  /* Start a Tk mainwindow to hold the widgets */
  if (Tk_Init(g_tk_control_interp) == TCL_ERROR) {
	  fprintf(stderr,
	  "Tk_Init failed: %s\n",g_tk_control_interp->result);
	  return(-1);
  }
  tk_control_window = Tk_MainWindow(g_tk_control_interp);
  if (tk_control_window == NULL) {
          fprintf(stderr,"%s\n", g_tk_control_interp->result);
          return(-1);
  }

  //------------------------------------------------------------------
  // Loading the particular definition files we need.  russ_widgets is
  // required by the Tclvar_float_with_scale class.  simple_magnet_drive
  // is application-specific and sets up the controls for the integer
  // and float variables.

  /* Load the Tcl scripts that handle widget definition and
   * variable controls */
  sprintf(command, "source russ_widgets.tcl");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }

  //------------------------------------------------------------------
  // Put the version number into the main window.
  sprintf(command, "label .versionlabel -text Template_matching_v:%s", Version_string);
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }
  sprintf(command, "pack .versionlabel");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }

  //------------------------------------------------------------------
  // Put the "Thank-you Ware" button into the main window.
  sprintf(command, "package require http");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }
  sprintf(command, "set thanks_text \"Say Thank You!\"");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }
  sprintf(command, "button .thankyouware -textvariable thanks_text -command { ::http::geturl \"http://www.cs.unc.edu/Research/nano/cismm/thankyou/yourewelcome.htm?program=cismm_template_based_matching&Version=%s\" ; set thanks_text \"Paid for by NIH/NIBIB\" }", Version_string);
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }
  sprintf(command, "pack .thankyouware");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }

  //------------------------------------------------------------------
  // Load the specialized Tcl code needed by this program.  This must
  // be loaded before the Tclvar_init() routine is called because it
  // puts together some of the windows needed by the variables.
  sprintf(command, "source template_based_matching.tcl");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }

  //------------------------------------------------------------------
  // This routine must be called in order to initialize all of the
  // variables that came into scope before the interpreter was set
  // up, and to tell the variables which interpreter to use.  It is
  // called once, after the interpreter exists.

  // Initialize the variables using the interpreter
  if (Tclvar_init(g_tk_control_interp)) {
	  fprintf(stderr,"Can't do init!\n");
	  return -1;
  }
  sprintf(command, "wm geometry . +10+10");
  if (Tcl_Eval(g_tk_control_interp, command) != TCL_OK) {
          fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", command,
                  g_tk_control_interp->result);
          return(-1);
  }

  //------------------------------------------------------------
  // Create a callback for a variable that will hold the data
  // name and then create a dialog box that will ask the user
  // to either fill it in or quit.
  Tclvar_selector filename("data_filename", NULL, NULL, "", data_filename_changed, NULL);

  //------------------------------------------------------------------
  // If we don't have any input images, then throw a Tcl dialog asking
  // the user for the name of a file to use and wait until they respond.
  if (g_current_data == NULL) {

    if (Tcl_Eval(g_tk_control_interp, "ask_user_for_data_filename") != TCL_OK) {
      fprintf(stderr, "Tcl_Eval(%s) failed: %s\n", "ask_user_for_data_filename", g_tk_control_interp->result);
      cleanup();
      exit(-1);
    }
    // Wait until we get a filename or are asked to quit
    do {
      // This pushes changes in the C variables over to Tcl.
      while (Tk_DoOneEvent(TK_DONT_WAIT)) {};
      if (Tclvar_mainloop()) {
	fprintf(stderr,"Tclvar Mainloop failed\n");
      }
    } while ( (g_current_data == NULL) && !g_quit);

    // If we were asked to quit, do so
    if (g_quit) {
      cleanup();
      exit(0);
    }
  }

  g_nColsSplat = g_nColsSecond = g_nColsData;
  g_nRowsSplat = g_nRowsSecond = g_nRowsData;

  //------------------------------------------------------------------
  // Set initial values for the flatten-list targets.
  g_flattenList.push_back(Position(g_nColsData/8, g_nRowsData/8));
  g_flattenList.push_back(Position(g_nColsData/2, g_nRowsData*7/8));
  g_flattenList.push_back(Position(g_nColsData*7/8, g_nRowsData/8));

  //------------------------------------------------------------------
  // Set up callback handler to change the data set that is shown in
  // the dataset display window.
  g_template_name.set_tcl_change_callback(handle_template_name_change, NULL);

  //------------------------------------------------------------------
  // We now have at least one image, so we can go ahead and display.
  g_ready_to_display = true;

  //------------------------------------------------------------------
  // Initialize GLUT and create the windows that will display the
  // data set and the second visualization.
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowPosition(175, 10);
  glutInitWindowSize(g_nColsData, g_nRowsData);
  g_data_window = glutCreateWindow("Original or Substrate-Fit Image");
  glutInitWindowPosition(190 + g_nColsData, 10);
  glutInitWindowSize(g_nColsSecond, g_nRowsSecond);
  g_second_window = glutCreateWindow("Flattened");
  glutInitWindowPosition(190 + g_nColsData, 10 + 60 + g_nRowsData);
  glutInitWindowSize(g_nColsSplat, g_nRowsSplat);
  g_splat_window = glutCreateWindow("Best Fit");

  // Create the buffers that Glut will use to display in the windows.  This is allocating an
  // RGBA buffer.  It needs to be 4-byte aligned, so we allocated it as a group of
  // words and then cast it to the right type.  We're using RGBA rather than just RGB
  // because it also solves the 4-byte alignment problem caused by funky sizes of image
  // that are RGB images.
  g_splat_image = new image_buffer(g_nColsSplat, g_nRowsSplat);
  if ( (g_splat_image == NULL) || !g_splat_image->doing_okay() ) {
    fprintf(stderr,"Out of memory when allocating image!\n");
    fprintf(stderr,"  (Image is %u by %u)\n", g_nColsSplat, g_nRowsSplat);
    cleanup();
  }
  g_second_image = new image_buffer(g_nColsSecond, g_nRowsSecond);
  if ( (g_second_image == NULL) || !g_second_image->doing_okay() ) {
    fprintf(stderr,"Out of memory when allocating image!\n");
    fprintf(stderr,"  (Image is %u by %u)\n", g_nColsSecond, g_nRowsSecond);
    cleanup();
  }
  g_data_image = new image_buffer(g_nColsData, g_nRowsData);
  if ( (g_data_image == NULL) || !g_data_image->doing_okay() ) {
    fprintf(stderr,"Out of memory when allocating image!\n");
    fprintf(stderr,"  (Image is %u by %u)\n", g_nColsData, g_nRowsData);
    cleanup();
  }

  //------------------------------------------------------------------------
  // Cause a redisplay/recalculation to get things started.
  float_causes_redisplay(0.0, NULL);

  //------------------------------------------------------------------------
  // Set up callback handlers and put Glut in control.
  glutSetWindow(g_splat_window);
  glutDisplayFunc(mySplatDisplayFunc);
  glutReshapeFunc(mySplatReshapeFunc);
  glutWindowStatusFunc(mySplatStatusFunc);

  glutSetWindow(g_second_window);
  glutDisplayFunc(mySecondDisplayFunc);
  glutReshapeFunc(mySecondReshapeFunc);
  glutWindowStatusFunc(mySecondStatusFunc);
  glutMouseFunc(secondMouseCallbackForGLUT);

  glutSetWindow(g_data_window);
  glutDisplayFunc(myDataDisplayFunc);
  glutReshapeFunc(myDataReshapeFunc);
  glutWindowStatusFunc(myDataStatusFunc);
  glutMotionFunc(dataMotionCallbackForGLUT);
  glutMouseFunc(dataMouseCallbackForGLUT);

  glutIdleFunc(myIdleFunc);

  //------------------------------------------------------------------------
  // Hand off control to Glut and do everything else in callback handlers
  // (not the idle handler).
  glutMainLoop();

  // We never get here because glutMainloop never returns.
  cleanup();
  return 0;
}
