// libplot.d
// This D module provides a binding to GNU libplot,
// a shared library for 2-dimensional vector graphics.
//
// It is mostly a copy and paste from the public header file plot.h
// so it retains its GPL3 licence.
//
// Peter J. 2018-03-11
//

module libplot;
import core.stdc.stdio;

extern (C):
enum char[8] PL_LIBPLOT_VER_STRING = "4.4";
enum uint PL_LIBPLOT_VER = 404;

extern const char[8] pl_libplot_ver;   /* need room for 99.99aa */

/* The functions in the C binding deal with `plPlotter' and
   `plPlotterParams' objects.  They are the same as the `Plotter' and
   `PlotterParams' objects of the C++ binding.  Internally, they are called
   `plPlotterStruct' and `plPlotterParamsStruct'.  In the context of this
   header file, they are opaque. */
struct plPlotterStruct {}
alias plPlotterStruct plPlotter;
struct plPlotterParamsStruct {}
alias plPlotterParamsStruct plPlotterParams;

/* THE C API */

/* Constructor/destructor for the plPlotter type.  Parameter values are
   specified at creation time via a plPlotterParams instance.  There is no
   copy constructor. */
plPlotter * pl_newpl_r (const(char) *type, FILE *infile, FILE *outfile,
                            FILE *errfile, const plPlotterParams *plotter_params);
int pl_deletepl_r (plPlotter *plotter);

/* Constructor/destructor/copy constructor for the plPlotterParams type,
   any instance of which stores parameters that are used when creating a
   plPlotter. */
plPlotterParams * pl_newplparams ();
int pl_deleteplparams (plPlotterParams *plotter_params);
plPlotterParams * pl_copyplparams (const plPlotterParams *plotter_params);

/* A function for setting a single Plotter parameter in a plPlotterParams
   instance.  */
int pl_setplparam (plPlotterParams *plotter_params, const(char) *parameter, void *value);

/* THE PLOTTER METHODS */

/* 13 functions in traditional (pre-GNU) libplot */
int pl_arc_r (plPlotter *plotter, int xc, int yc, int x0, int y0, int x1, int y1);
int pl_box_r (plPlotter *plotter, int x0, int y0, int x1, int y1);
int pl_circle_r (plPlotter *plotter, int x, int y, int r);
int pl_closepl_r (plPlotter *plotter);
int pl_cont_r (plPlotter *plotter, int x, int y);
int pl_erase_r (plPlotter *plotter);
int pl_label_r (plPlotter *plotter, const(char) *s);
int pl_line_r (plPlotter *plotter, int x0, int y0, int x1, int y1);
int pl_linemod_r (plPlotter *plotter, const(char) *s);
int pl_move_r (plPlotter *plotter, int x, int y);
int pl_openpl_r (plPlotter *plotter);
int pl_point_r (plPlotter *plotter, int x, int y);
int pl_space_r (plPlotter *plotter, int x0, int y0, int x1, int y1);

/* 46 additional functions in GNU libplot, plus 1 obsolete function
   [pl_outfile_r]. */
FILE* pl_outfile_r (plPlotter *plotter, FILE* outfile);/* OBSOLETE */
int pl_alabel_r (plPlotter *plotter, int x_justify, int y_justify, const(char) *s);
int pl_arcrel_r (plPlotter *plotter, int dxc, int dyc, int dx0, int dy0, int dx1, int dy1);
int pl_bezier2_r (plPlotter *plotter, int x0, int y0, int x1, int y1, int x2, int y2);
int pl_bezier2rel_r (plPlotter *plotter, int dx0, int dy0, int dx1, int dy1, int dx2, int dy2);
int pl_bezier3_r (plPlotter *plotter, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3);
int pl_bezier3rel_r (plPlotter *plotter, int dx0, int dy0, int dx1, int dy1, int dx2, int dy2, int dx3, int dy3);
int pl_bgcolor_r (plPlotter *plotter, int red, int green, int blue);
int pl_bgcolorname_r (plPlotter *plotter, const(char) *name);
int pl_boxrel_r (plPlotter *plotter, int dx0, int dy0, int dx1, int dy1);
int pl_capmod_r (plPlotter *plotter, const(char) *s);
int pl_circlerel_r (plPlotter *plotter, int dx, int dy, int r);
int pl_closepath_r (plPlotter *plotter);
int pl_color_r (plPlotter *plotter, int red, int green, int blue);
int pl_colorname_r (plPlotter *plotter, const(char) *name);
int pl_contrel_r (plPlotter *plotter, int x, int y);
int pl_ellarc_r (plPlotter *plotter, int xc, int yc, int x0, int y0, int x1, int y1);
int pl_ellarcrel_r (plPlotter *plotter, int dxc, int dyc, int dx0, int dy0, int dx1, int dy1);
int pl_ellipse_r (plPlotter *plotter, int x, int y, int rx, int ry, int angle);
int pl_ellipserel_r (plPlotter *plotter, int dx, int dy, int rx, int ry, int angle);
int pl_endpath_r (plPlotter *plotter);
int pl_endsubpath_r (plPlotter *plotter);
int pl_fillcolor_r (plPlotter *plotter, int red, int green, int blue);
int pl_fillcolorname_r (plPlotter *plotter, const(char) *name);
int pl_fillmod_r (plPlotter *plotter, const(char) *s);
int pl_filltype_r (plPlotter *plotter, int level);
int pl_flushpl_r (plPlotter *plotter);
int pl_fontname_r (plPlotter *plotter, const(char) *s);
int pl_fontsize_r (plPlotter *plotter, int size);
int pl_havecap_r (plPlotter *plotter, const(char) *s);
int pl_joinmod_r (plPlotter *plotter, const(char) *s);
int pl_labelwidth_r (plPlotter *plotter, const(char) *s);
int pl_linedash_r (plPlotter *plotter, int n, const int *dashes, int offset);
int pl_linerel_r (plPlotter *plotter, int dx0, int dy0, int dx1, int dy1);
int pl_linewidth_r (plPlotter *plotter, int size);
int pl_marker_r (plPlotter *plotter, int x, int y, int type, int size);
int pl_markerrel_r (plPlotter *plotter, int dx, int dy, int type, int size);
int pl_moverel_r (plPlotter *plotter, int x, int y);
int pl_orientation_r (plPlotter *plotter, int direction);
int pl_pencolor_r (plPlotter *plotter, int red, int green, int blue);
int pl_pencolorname_r (plPlotter *plotter, const(char) *name);
int pl_pentype_r (plPlotter *plotter, int level);
int pl_pointrel_r (plPlotter *plotter, int dx, int dy);
int pl_restorestate_r (plPlotter *plotter);
int pl_savestate_r (plPlotter *plotter);
int pl_space2_r (plPlotter *plotter, int x0, int y0, int x1, int y1, int x2, int y2);
int pl_textangle_r (plPlotter *plotter, int angle);

/* 32 floating point counterparts to some of the above (all GNU additions) */
double pl_ffontname_r (plPlotter *plotter, const(char) *s);
double pl_ffontsize_r (plPlotter *plotter, double size);
double pl_flabelwidth_r (plPlotter *plotter, const(char) *s);
double pl_ftextangle_r (plPlotter *plotter, double angle);
int pl_farc_r (plPlotter *plotter, double xc, double yc, double x0, double y0, double x1, double y1);
int pl_farcrel_r (plPlotter *plotter, double dxc, double dyc, double dx0, double dy0, double dx1, double dy1);
int pl_fbezier2_r (plPlotter *plotter, double x0, double y0, double x1, double y1, double x2, double y2);
int pl_fbezier2rel_r (plPlotter *plotter, double dx0, double dy0, double dx1, double dy1, double dx2, double dy2);
int pl_fbezier3_r (plPlotter *plotter, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);
int pl_fbezier3rel_r (plPlotter *plotter, double dx0, double dy0, double dx1, double dy1, double dx2, double dy2, double dx3, double dy3);
int pl_fbox_r (plPlotter *plotter, double x0, double y0, double x1, double y1);
int pl_fboxrel_r (plPlotter *plotter, double dx0, double dy0, double dx1, double dy1);
int pl_fcircle_r (plPlotter *plotter, double x, double y, double r);
int pl_fcirclerel_r (plPlotter *plotter, double dx, double dy, double r);
int pl_fcont_r (plPlotter *plotter, double x, double y);
int pl_fcontrel_r (plPlotter *plotter, double dx, double dy);
int pl_fellarc_r (plPlotter *plotter, double xc, double yc, double x0, double y0, double x1, double y1);
int pl_fellarcrel_r (plPlotter *plotter, double dxc, double dyc, double dx0, double dy0, double dx1, double dy1);
int pl_fellipse_r (plPlotter *plotter, double x, double y, double rx, double ry, double angle);
int pl_fellipserel_r (plPlotter *plotter, double dx, double dy, double rx, double ry, double angle);
int pl_flinedash_r (plPlotter *plotter, int n, const double *dashes, double offset);
int pl_fline_r (plPlotter *plotter, double x0, double y0, double x1, double y1);
int pl_flinerel_r (plPlotter *plotter, double dx0, double dy0, double dx1, double dy1);
int pl_flinewidth_r (plPlotter *plotter, double size);
int pl_fmarker_r (plPlotter *plotter, double x, double y, int type, double size);
int pl_fmarkerrel_r (plPlotter *plotter, double dx, double dy, int type, double size);
int pl_fmove_r (plPlotter *plotter, double x, double y);
int pl_fmoverel_r (plPlotter *plotter, double dx, double dy);
int pl_fpoint_r (plPlotter *plotter, double x, double y);
int pl_fpointrel_r (plPlotter *plotter, double dx, double dy);
int pl_fspace_r (plPlotter *plotter, double x0, double y0, double x1, double y1);
int pl_fspace2_r (plPlotter *plotter, double x0, double y0, double x1, double y1, double x2, double y2);

/* 6 floating point operations with no integer counterpart (GNU additions) */
int pl_fconcat_r (plPlotter *plotter, double m0, double m1, double m2, double m3, double m4, double m5);
int pl_fmiterlimit_r (plPlotter *plotter, double limit);
int pl_frotate_r (plPlotter *plotter, double theta);
int pl_fscale_r (plPlotter *plotter, double x, double y);
int pl_fsetmatrix_r (plPlotter *plotter, double m0, double m1, double m2, double m3, double m4, double m5);
int pl_ftranslate_r (plPlotter *plotter, double x, double y);

/* THE OLD (non-thread-safe) C API */

/* 3 functions specific to the old C API.  (For construction/destruction
   and selection of Plotters, and setting of Plotter parameters.  The fact
   that a single Plotter is globally `selected' makes the old API
   non-thread-safe.) */
int pl_newpl (const(char) *type, FILE *infile, FILE *outfile, FILE *errfile);
int pl_selectpl (int handle);
int pl_deletepl (int handle);

/* A function for setting parameters of Plotters that will subsequently be
   created.  This also makes the old API non-thread-safe. */
int pl_parampl (const(char) *parameter, void *value);

/* THE PLOTTER METHODS */
/* In the old API, the Plotter to be acted on is specified by first calling 
   selectpl(). */

/* 13 functions in traditional (pre-GNU) libplot */
int pl_arc (int xc, int yc, int x0, int y0, int x1, int y1);
int pl_box (int x0, int y0, int x1, int y1);
int pl_circle (int x, int y, int r);
int pl_closepl ();
int pl_cont (int x, int y);
int pl_erase ();
int pl_label (const(char) *s);
int pl_line (int x0, int y0, int x1, int y1);
int pl_linemod (const(char) *s);
int pl_move (int x, int y);
int pl_openpl ();
int pl_point (int x, int y);
int pl_space (int x0, int y0, int x1, int y1);

/* 46 additional functions in GNU libplot, plus 1 obsolete function
   [pl_outfile]. */
FILE* pl_outfile (FILE* outfile);/* OBSOLETE */
int pl_alabel (int x_justify, int y_justify, const(char) *s);
int pl_arcrel (int dxc, int dyc, int dx0, int dy0, int dx1, int dy1);
int pl_bezier2 (int x0, int y0, int x1, int y1, int x2, int y2);
int pl_bezier2rel (int dx0, int dy0, int dx1, int dy1, int dx2, int dy2);
int pl_bezier3 (int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3);
int pl_bezier3rel (int dx0, int dy0, int dx1, int dy1, int dx2, int dy2, int dx3, int dy3);
int pl_bgcolor (int red, int green, int blue);
int pl_bgcolorname (const(char) *name);
int pl_boxrel (int dx0, int dy0, int dx1, int dy1);
int pl_capmod (const(char) *s);
int pl_circlerel (int dx, int dy, int r);
int pl_closepath ();
int pl_color (int red, int green, int blue);
int pl_colorname (const(char) *name);
int pl_contrel (int x, int y);
int pl_ellarc (int xc, int yc, int x0, int y0, int x1, int y1);
int pl_ellarcrel (int dxc, int dyc, int dx0, int dy0, int dx1, int dy1);
int pl_ellipse (int x, int y, int rx, int ry, int angle);
int pl_ellipserel (int dx, int dy, int rx, int ry, int angle);
int pl_endpath ();
int pl_endsubpath ();
int pl_fillcolor (int red, int green, int blue);
int pl_fillcolorname (const(char) *name);
int pl_fillmod (const(char) *s);
int pl_filltype (int level);
int pl_flushpl ();
int pl_fontname (const(char) *s);
int pl_fontsize (int size);
int pl_havecap (const(char) *s);
int pl_joinmod (const(char) *s);
int pl_labelwidth (const(char) *s);
int pl_linedash (int n, const int *dashes, int offset);
int pl_linerel (int dx0, int dy0, int dx1, int dy1);
int pl_linewidth (int size);
int pl_marker (int x, int y, int type, int size);
int pl_markerrel (int dx, int dy, int type, int size);
int pl_moverel (int x, int y);
int pl_orientation (int direction);
int pl_pencolor (int red, int green, int blue);
int pl_pencolorname (const(char) *name);
int pl_pentype (int level);
int pl_pointrel (int dx, int dy);
int pl_restorestate ();
int pl_savestate ();
int pl_space2 (int x0, int y0, int x1, int y1, int x2, int y2);
int pl_textangle (int angle);

/* 32 floating point counterparts to some of the above (all GNU additions) */
double pl_ffontname (const(char) *s);
double pl_ffontsize (double size);
double pl_flabelwidth (const(char) *s);
double pl_ftextangle (double angle);
int pl_farc (double xc, double yc, double x0, double y0, double x1, double y1);
int pl_farcrel (double dxc, double dyc, double dx0, double dy0, double dx1, double dy1);
int pl_fbezier2 (double x0, double y0, double x1, double y1, double x2, double y2);
int pl_fbezier2rel (double dx0, double dy0, double dx1, double dy1, double dx2, double dy2);
int pl_fbezier3 (double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);
int pl_fbezier3rel (double dx0, double dy0, double dx1, double dy1, double dx2, double dy2, double dx3, double dy3);
int pl_fbox (double x0, double y0, double x1, double y1);
int pl_fboxrel (double dx0, double dy0, double dx1, double dy1);
int pl_fcircle (double x, double y, double r);
int pl_fcirclerel (double dx, double dy, double r);
int pl_fcont (double x, double y);
int pl_fcontrel (double dx, double dy);
int pl_fellarc (double xc, double yc, double x0, double y0, double x1, double y1);
int pl_fellarcrel (double dxc, double dyc, double dx0, double dy0, double dx1, double dy1);
int pl_fellipse (double x, double y, double rx, double ry, double angle);
int pl_fellipserel (double dx, double dy, double rx, double ry, double angle);
int pl_flinedash (int n, const(double)* dashes, double offset);
int pl_fline (double x0, double y0, double x1, double y1);
int pl_flinerel (double dx0, double dy0, double dx1, double dy1);
int pl_flinewidth (double size);
int pl_fmarker (double x, double y, int type, double size);
int pl_fmarkerrel (double dx, double dy, int type, double size);
int pl_fmove (double x, double y);
int pl_fmoverel (double dx, double dy);
int pl_fpoint (double x, double y);
int pl_fpointrel (double dx, double dy);
int pl_fspace (double x0, double y0, double x1, double y1);
int pl_fspace2 (double x0, double y0, double x1, double y1, double x2, double y2);

/* 6 floating point operations with no integer counterpart (GNU additions) */
int pl_fconcat (double m0, double m1, double m2, double m3, double m4, double m5);
int pl_fmiterlimit (double limit);
int pl_frotate (double theta);
int pl_fscale (double x, double y);
int pl_fsetmatrix (double m0, double m1, double m2, double m3, double m4, double m5);
int pl_ftranslate (double x, double y);


/* UNDOCUMENTED FONT API CALLS */
/* These are used by the graphics programs in the plotutils package (e.g.,
   `graph') to access the font tables within libplot, so that the user can
   be given lists of font names. */

void *_pl_get_hershey_font_info (plPlotter *plotter);
void *_pl_get_ps_font_info (plPlotter *plotter);
void *_pl_get_pcl_font_info (plPlotter *plotter);
void *_pl_get_stick_font_info (plPlotter *plotter);


/* THE GLOBAL VARIABLES IN GNU LIBPLOT */
/* There are two: user-settable error handlers (not yet documented). */
// extern (C) int function(const(char) *msg) pl_libplot_warning_handler;
// extern (C) int function(const(char) *msg) pl_libplot_error_handler;
// If we include the two lines above, we get the linker error:
// /usr/bin/ld: pl_libplot_warning_handler:
//     TLS definition in drifting_eye.o section .tbss mismatches non-TLS definition
//     in /usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib/libplot.so section .bss
// In the original /usr/include/plot.h these lines were:
// extern int (*pl_libplot_warning_handler) (const char *msg);
// extern int (*pl_libplot_error_handler) (const char *msg);

/* Useful definitions, included in both plot.h and plotter.h. */

/* Symbol types for the marker() function, extending over the range 0..31.
   (1 through 5 are the same as in the GKS [Graphical Kernel System].)

   These are now defined as enums rather than ints.  Cast them to ints if
   necessary. */
enum 
{ M_NONE, M_DOT, M_PLUS, M_ASTERISK, M_CIRCLE, M_CROSS, 
  M_SQUARE, M_TRIANGLE, M_DIAMOND, M_STAR, M_INVERTED_TRIANGLE, 
  M_STARBURST, M_FANCY_PLUS, M_FANCY_CROSS, M_FANCY_SQUARE, 
  M_FANCY_DIAMOND, M_FILLED_CIRCLE, M_FILLED_SQUARE, M_FILLED_TRIANGLE, 
  M_FILLED_DIAMOND, M_FILLED_INVERTED_TRIANGLE, M_FILLED_FANCY_SQUARE,
  M_FILLED_FANCY_DIAMOND, M_HALF_FILLED_CIRCLE, M_HALF_FILLED_SQUARE,
  M_HALF_FILLED_TRIANGLE, M_HALF_FILLED_DIAMOND,
  M_HALF_FILLED_INVERTED_TRIANGLE, M_HALF_FILLED_FANCY_SQUARE,
  M_HALF_FILLED_FANCY_DIAMOND, M_OCTAGON, M_FILLED_OCTAGON 
};

/* ONE-BYTE OPERATION CODES FOR GNU METAFILE FORMAT. These are now defined
   as enums rather than ints.  Cast them to ints if necessary.

   There are 85 currently recognized op codes.  The first 10 date back to
   Unix plot(5) format. */

enum
{  
/* 10 op codes for primitive graphics operations, as in Unix plot(5) format. */
  O_ARC		=	'a',  
  O_CIRCLE	=	'c',  
  O_CONT	=	'n',
  O_ERASE	=	'e',
  O_LABEL	=	't',
  O_LINEMOD	=	'f',
  O_LINE	=	'l',
  O_MOVE	=	'm',
  O_POINT	=	'p',
  O_SPACE	=	's',
  
/* 42 op codes that are GNU extensions */
  O_ALABEL	=	'T',
  O_ARCREL	=	'A',
  O_BEZIER2	=       'q',
  O_BEZIER2REL	=       'r',
  O_BEZIER3	=       'y',
  O_BEZIER3REL	=       'z',
  O_BGCOLOR	=	'~',
  O_BOX		=	'B',	/* not an op code in Unix plot(5) */
  O_BOXREL	=	'H',
  O_CAPMOD	=	'K',
  O_CIRCLEREL	=	'G',
  O_CLOSEPATH	=	'k',
  O_CLOSEPL	=	'x',	/* not an op code in Unix plot(5) */
  O_COMMENT	=	'#',
  O_CONTREL	=	'N',
  O_ELLARC	=	'?',
  O_ELLARCREL	=	'/',
  O_ELLIPSE	=	'+',
  O_ELLIPSEREL	=	'=',
  O_ENDPATH	=	'E',
  O_ENDSUBPATH	=	']',
  O_FILLTYPE	=	'L',
  O_FILLCOLOR	=	'D',
  O_FILLMOD	=	'g',
  O_FONTNAME	=	'F',
  O_FONTSIZE	=	'S',
  O_JOINMOD	=	'J',
  O_LINEDASH	= 	'd',
  O_LINEREL	=	'I',
  O_LINEWIDTH	=	'W',
  O_MARKER	=	'Y',
  O_MARKERREL	=	'Z',
  O_MOVEREL	=	'M',
  O_OPENPL	=	'o',	/* not an op code in Unix plot(5) */
  O_ORIENTATION	=	'b',
  O_PENCOLOR	=	'-',
  O_PENTYPE	=	'h',
  O_POINTREL	=	'P',
  O_RESTORESTATE=	'O',
  O_SAVESTATE	=	'U',
  O_SPACE2	=	':',
  O_TEXTANGLE	=	'R',

/* 30 floating point counterparts to many of the above.  They are not even
   slightly mnemonic. */
  O_FARC	=	'1',
  O_FARCREL	=	'2',
  O_FBEZIER2	=       '`',
  O_FBEZIER2REL	=       '\'',
  O_FBEZIER3	=       ',',
  O_FBEZIER3REL	=       '.',
  O_FBOX	=	'3',
  O_FBOXREL	=	'4',
  O_FCIRCLE	=	'5',
  O_FCIRCLEREL	=	'6',
  O_FCONT	=	')',
  O_FCONTREL	=	'_',
  O_FELLARC	=	'}',
  O_FELLARCREL	=	'|',
  O_FELLIPSE	=	'{',
  O_FELLIPSEREL	=	'[',
  O_FFONTSIZE	=	'7',
  O_FLINE	=	'8',
  O_FLINEDASH	= 	'w',
  O_FLINEREL	=	'9',
  O_FLINEWIDTH	=	'0',
  O_FMARKER	=	'!',
  O_FMARKERREL	=	'@',
  O_FMOVE	=	'$',
  O_FMOVEREL	=	'%',
  O_FPOINT	=	'^',
  O_FPOINTREL	=	'&',
  O_FSPACE	=	'*',
  O_FSPACE2	=	';',
  O_FTEXTANGLE	=	'(',

/* 3 op codes for floating point operations with no integer counterpart */
  O_FCONCAT		=	'\\',
  O_FMITERLIMIT		=	'i',
  O_FSETMATRIX		=	'j'
};
