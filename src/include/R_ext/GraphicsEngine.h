/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-8 The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* Used by graphics.c, grid and by third-party graphics devices */

#ifndef R_GRAPHICSENGINE_H_
#define R_GRAPHICSENGINE_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * The current graphics engine (including graphics device) API version
 * MUST be integer
 * 
 * This number should be bumped whenever there are changes to 
 * GraphicsEngine.h or GraphicsDevice.h so that add-on packages
 * that compile against these headers (graphics systems such as
 * graphics and grid;  graphics devices such as gtkDevice, RSvgDevice)
 * can detect any version mismatch.
 *
 * Version 1:  Introduction of the version number.
 * Version 2:  GEDevDesc *dd dropped from GEcontourLines().
 * Version 3:  R_GE_str2col() added to API.
 * Version 4:  UTF-8 text hooks, useRotatedTextInContour,
 *             add newFrameConfirm() to NewDevDesc.
 *             New API: GEaddDevice[2] GEgetDevice, GEkillDevice,
 *             ndevNumber.
 * Version 5:  Clean up 1.4.0/2.0.0 changes!
 *             Remove newDevStruct from GEDevDesc and NewDevDesc.
 *             Remove asp, dot(), hold(), open() from NewDevDesc.
 *             Move displayList, DLlastElt, savedSnapshot from
 *             NewDevDesc to GEDevDesc.
 *             Add 'ask' to GEDevDesc.
 */

#define R_GE_version 5

int R_GE_getVersion(void);

void R_GE_checkVersionOrDie(int version);

/* The graphics engine will only accept locations and dimensions 
 * in native device coordinates, but it provides the following functions
 * for converting between a couple of simple alternative coordinate 
 * systems and device coordinates:
 *    DEVICE = native units of the device
 *    NDC = Normalised device coordinates 
 *    INCHES = inches (!)
 *    CM = centimetres (!!)
 */

typedef enum {
 GE_DEVICE	= 0,	/* native device coordinates (rasters) */
 GE_NDC	= 1,	/* normalised device coordinates x=(0,1), y=(0,1) */
 GE_INCHES = 2,
 GE_CM     = 3
} GEUnit;

#define MAX_GRAPHICS_SYSTEMS 24

typedef enum {
    /* In response to this event, the registered graphics system
     * should allocate and initialise the systemSpecific structure
     */
    GE_InitState = 0,
    /* This event gives the registered system a chance to undo
     * anything done in the initialisation.
     */
    GE_FinaliseState = 1,
    /* This is sent by the graphics engine prior to initialising
     * the display list.  It give the graphics system the chance
     * to squirrel away information it will need for redrawing the
     * the display list
     */
    GE_SaveState = 2,
    /* This is sent by the graphics engine prior to replaying the
     * display list.  It gives the graphics system the chance to
     * restore any information it saved on the GE_SaveState event
     */
    GE_RestoreState = 6,
    /* Copy system state information to the current device.
     * This is used when copying graphics from one device to another
     * so all the graphics system needs to do is to copy across
     * the bits required for the display list to draw faithfully
     * on the new device.
     */
    GE_CopyState = 3,
    /* Create a snapshot of the system state that is sufficient
     * for the current "image" to be reproduced
     */
    GE_SaveSnapshotState = 4,
    /* Restore the system state that is saved by GE_SaveSnapshotState
     */
    GE_RestoreSnapshotState = 5,
    /* When replaying the display list, the graphics engine 
     * checks, after each replayed action, that the action 
     * produced valid output.  This is the graphics system's
     * chance to say that the output is crap (in which case the
     * graphics engine will abort the display list replay).
     */
    GE_CheckPlot = 7,
    /* The device wants to scale the current pointsize
     * (for scaling an image)
     * This is not a nice general solution, but a quick fix for 
     * the Windows device.
     */
    GE_ScalePS = 8
} GEevent;

/*
 *  Some line end/join constants
 */
typedef enum {
  GE_ROUND_CAP  = 1,
  GE_BUTT_CAP   = 2,
  GE_SQUARE_CAP = 3
} R_GE_lineend;

typedef enum {
  GE_ROUND_JOIN = 1,
  GE_MITRE_JOIN = 2,
  GE_BEVEL_JOIN = 3
} R_GE_linejoin;

/* 
 * A structure containing graphical parameters 
 *
 * This is how graphical parameters are passed from graphics systems
 * to the graphics engine AND from the graphics engine to graphics
 * devices.
 *
 * Devices are not *required* to honour graphical parameters
 * (e.g., alpha transparency is going to be tough for some)
 */
typedef struct {
    /*
     * Colours
     *
     * NOTE:  Alpha transparency included in col & fill
     */
    int col;             /* pen colour (lines, text, borders, ...) */
    int fill;            /* fill colour (for polygons, circles, rects, ...) */
    double gamma;        /* Gamma correction */
    /* 
     * Line characteristics
     */
    double lwd;          /* Line width (roughly number of pixels) */
    int lty;             /* Line type (solid, dashed, dotted, ...) */
    R_GE_lineend lend;   /* Line end */
    R_GE_linejoin ljoin; /* line join */
    double lmitre;       /* line mitre */
    /*
     * Text characteristics
     */
    double cex;          /* Character expansion (font size = fontsize*cex) */
    double ps;           /* Font size in points */
    double lineheight;   /* Line height (multiply by font size) */
    int fontface;        /* Font face (plain, italic, bold, ...) */
    char fontfamily[201]; /* Font family */
} R_GE_gcontext;

typedef R_GE_gcontext* pGEcontext;
    

#include <R_ext/GraphicsDevice.h> /* needed for NewDevDesc */

typedef struct _GEDevDesc GEDevDesc;

typedef SEXP (* GEcallback)(GEevent, GEDevDesc *, SEXP);

typedef struct { 
    /* An array of information about each graphics system that
     * has registered with the graphics engine.
     * This is used to store graphics state for each graphics
     * system on each device.
     */
    void *systemSpecific;
    /* 
     * An array of function pointers, one per graphics system that
     * has registered with the graphics engine.
     *
     * system_Callback is called when the graphics engine wants
     * to give a graphics system the chance to play with its
     * device-specific information (stored in systemSpecific)
     * There are two parameters:  an "event" to tell the graphics
     * system why the graphics engine has called this function,
     * and the systemSpecific pointer.  The graphics engine
     * has to pass the systemSpecific pointer because only
     * the graphics engine will know what array index to use.
     */
    GEcallback callback;
} GESystemDesc;

struct _GEDevDesc {
    /* 
     * Stuff that the devices can see (and modify).
     * All detailed in GraphicsDevice.h
     */
    pDevDesc dev;
    /*
     * Stuff about the device that only the graphics engine sees
     * (the devices don't see it).
     */
    Rboolean displayListOn;  /* toggle for display list status */
    SEXP displayList;        /* display list */
    SEXP DLlastElt;          /* A pointer to the end of the display list
				to avoid tranversing pairlists */
    SEXP savedSnapshot;      /* The last element of the display list
			      * just prior to when the display list
			      * was last initialised
			      */
    Rboolean dirty;          /* Has the device received any output? */
    Rboolean recordGraphics; /* Should a graphics call be stored
			      * on the display list?
			      * Set to FALSE by do_recordGraphics,
			      * do_dotcallgr, and do_Externalgr 
			      * so that nested calls are not
			      * recorded on the display list
			      */
    /* 
     * Stuff about the device that only graphics systems see.
     * The graphics engine has no idea what is in here.
     * Used by graphics systems to store system state per device.
     */
    GESystemDesc *gesd[MAX_GRAPHICS_SYSTEMS];

    /* per-device setting for 'ask' (use NewFrameConfirm) */
    Rboolean ask;
};

typedef GEDevDesc* pGEDevDesc;

/* functions from devices.c for use by graphics devices */

#define desc2GEDesc		Rf_desc2GEDesc
/* map NewDevDesc to enclosing GEDevDesc */
pGEDevDesc desc2GEDesc(pDevDesc dd);
int GEdeviceNumber(pGEDevDesc);
pGEDevDesc GEgetDevice(int);
void GEaddDevice(pGEDevDesc);
void GEaddDevice2(pGEDevDesc, const char *);
void GEkillDevice(pGEDevDesc);
pGEDevDesc GEcreateDevDesc(pDevDesc dev);

void GEdestroyDevDesc(pGEDevDesc dd);
void *GEsystemState(pGEDevDesc dd, int index);
void GEregisterWithDevice(pGEDevDesc dd);
void GEregisterSystem(GEcallback callback, int *systemRegisterIndex);
void GEunregisterSystem(int registerIndex);
SEXP GEhandleEvent(GEevent event, pDevDesc dev, SEXP data);

#define fromDeviceX		GEfromDeviceX
#define toDeviceX		GEtoDeviceX
#define fromDeviceY		GEfromDeviceY
#define toDeviceY		GEtoDeviceY
#define fromDeviceWidth		GEfromDeviceWidth
#define toDeviceWidth		GEtoDeviceWidth
#define fromDeviceHeight	GEfromDeviceHeight
#define toDeviceHeight		GEtoDeviceHeight

double fromDeviceX(double value, GEUnit to, pGEDevDesc dd); 
double toDeviceX(double value, GEUnit from, pGEDevDesc dd);
double fromDeviceY(double value, GEUnit to, pGEDevDesc dd); 
double toDeviceY(double value, GEUnit from, pGEDevDesc dd);
double fromDeviceWidth(double value, GEUnit to, pGEDevDesc dd); 
double toDeviceWidth(double value, GEUnit from, pGEDevDesc dd);
double fromDeviceHeight(double value, GEUnit to, pGEDevDesc dd); 
double toDeviceHeight(double value, GEUnit from, pGEDevDesc dd);

/*-------------------------------------------------------------------
 *
 *  COLOUR CODE is concerned with the internals of R colour representation
 *
 *  From colors.c, used in par.c, grid/src/gpar.c
 */

#define RGBpar			Rf_RGBpar
#define RGBpar3			Rf_RGBpar3
#define col2name                Rf_col2name
#define name2col		Rf_name2col

/* Convert an element of a R colour specification (which might be a
   number or a string) into an internal colour specification. */
unsigned int RGBpar(SEXP, int);
unsigned int RGBpar3(SEXP, int, unsigned int);

/* Convert an internal colour specification to/from a colour name */
const char *col2name(unsigned int col); /* used in par.c, grid */
unsigned int name2col(const char *);    /* used by plotmath.c */

/* Convert either a name or a #RRGGBB[AA] string to internal.
   Because people were using it, it also converts "1", "2" ...
   to a colour in the palette, and "0" to transparent white.
*/
unsigned int R_GE_str2col(const char *s);



/*
 *	Some Notes on Line Textures
 *
 *	Line textures are stored as an array of 4-bit integers within
 *	a single 32-bit word.  These integers contain the lengths of
 *	lines to be drawn with the pen alternately down and then up.
 *	The device should try to arrange that these values are measured
 *	in points if possible, although pixels is ok on most displays.
 *
 *	If newlty contains a line texture description it is decoded
 *	as follows:
 *
 *		ndash = 0;
 *		for(i=0 ; i<8 && newlty & 15 ; i++) {
 *			dashlist[ndash++] = newlty & 15;
 *			newlty = newlty>>4;
 *		}
 *		dashlist[0] = length of pen-down segment
 *		dashlist[1] = length of pen-up segment
 *		etc
 *
 *	An integer containing a zero terminates the pattern.  Hence
 *	ndash in this code fragment gives the length of the texture
 *	description.  If a description contains an odd number of
 *	elements it is replicated to create a pattern with an
 *	even number of elements.  (If this is a pain, do something
 *	different its not crucial).
 *
 */

/*--- The basic numbered & names line types; Here device-independent:
  e.g. "dashed" == "44",  "dotdash" == "1343"
*/

/* NB: was also in Rgraphics.h in R < 2.7.0 */
#define LTY_BLANK	-1
#define LTY_SOLID	0
#define LTY_DASHED	4 + (4<<4)
#define LTY_DOTTED	1 + (3<<4)
#define LTY_DOTDASH	1 + (3<<4) + (4<<8) + (3<<12)
#define LTY_LONGDASH	7 + (3<<4)
#define LTY_TWODASH	2 + (2<<4) + (6<<8) + (2<<12)

R_GE_lineend GE_LENDpar(SEXP value, int ind);
SEXP GE_LENDget(R_GE_lineend lend);
R_GE_linejoin GE_LJOINpar(SEXP value, int ind);
SEXP GE_LJOINget(R_GE_linejoin ljoin);

void GESetClip(double x1, double y1, double x2, double y2, pGEDevDesc dd);
void GENewPage(const pGEcontext gc, pGEDevDesc dd);
void GELine(double x1, double y1, double x2, double y2, 
	    const pGEcontext gc, pGEDevDesc dd);
void GEPolyline(int n, double *x, double *y, 
		const pGEcontext gc, pGEDevDesc dd);
void GEPolygon(int n, double *x, double *y, 
	       const pGEcontext gc, pGEDevDesc dd);
SEXP GEXspline(int n, double *x, double *y, double *s, Rboolean open, 
	       Rboolean repEnds, Rboolean draw,
	       const pGEcontext gc, pGEDevDesc dd);
void GECircle(double x, double y, double radius,
	      const pGEcontext gc, pGEDevDesc dd);
void GERect(double x0, double y0, double x1, double y1,
	    const pGEcontext gc, pGEDevDesc dd);
void GEText(double x, double y, const char * const str, cetype_t enc,
	    double xc, double yc, double rot,
	    const pGEcontext gc, pGEDevDesc dd);
void GEMode(int mode, pGEDevDesc dd);
void GESymbol(double x, double y, int pch, double size,
	      const pGEcontext gc, pGEDevDesc dd);
void GEPretty(double *lo, double *up, int *ndiv);
void GEMetricInfo(int c, const pGEcontext gc, 
		  double *ascent, double *descent, double *width,
		  pGEDevDesc dd);
double GEStrWidth(const char *str, cetype_t enc, 
		  const pGEcontext gc, pGEDevDesc dd);
double GEStrHeight(const char *str, cetype_t enc,
		  const pGEcontext gc, pGEDevDesc dd);
int GEstring_to_pch(SEXP pch);

/*-------------------------------------------------------------------
 *
 *  LINE TEXTURE CODE is concerned with the internals of R
 *  line texture representation.
 */
unsigned int GE_LTYpar(SEXP, int);
SEXP GE_LTYget(unsigned int);

/* 
 * From plotmath.c 
 */
double GEExpressionWidth(SEXP expr, 
			 const pGEcontext gc, pGEDevDesc dd);
double GEExpressionHeight(SEXP expr, 
			  const pGEcontext gc, pGEDevDesc dd);
void GEMathText(double x, double y, SEXP expr,
		double xc, double yc, double rot, 
		const pGEcontext gc, pGEDevDesc dd);
/* 
 * (End from plotmath.c)
 */

/* 
 * From plot3d.c : used in package clines
 */
SEXP GEcontourLines(double *x, int nx, double *y, int ny,
		    double *z, double *levels, int nl);
/* 
 * (End from plot3d.c)
 */

/* 
 * From vfonts.c
 */
double R_GE_VStrWidth(const char *s, cetype_t enc, const pGEcontext gc, pGEDevDesc dd);

double R_GE_VStrHeight(const char *s, cetype_t enc, const pGEcontext gc, pGEDevDesc dd);
void R_GE_VText(double x, double y, const char * const s, cetype_t enc, 
		double x_justify, double y_justify, double rotation,
		const pGEcontext gc, pGEDevDesc dd);
/* 
 * (End from vfonts.c)
 */

/* Also in Graphics.h */
#define	DEG2RAD 0.01745329251994329576

pGEDevDesc GEcurrentDevice(void);
Rboolean GEdeviceDirty(pGEDevDesc dd);
void GEdirtyDevice(pGEDevDesc dd);
Rboolean GEcheckState(pGEDevDesc dd);
Rboolean GErecording(SEXP call, pGEDevDesc dd);
void GErecordGraphicOperation(SEXP op, SEXP args, pGEDevDesc dd);
void GEinitDisplayList(pGEDevDesc dd);
void GEplayDisplayList(pGEDevDesc dd);
void GEcopyDisplayList(int fromDevice);
SEXP GEcreateSnapshot(pGEDevDesc dd);
void GEplaySnapshot(SEXP snapshot, pGEDevDesc dd);
void GEonExit(void);
void GEnullDevice(void);


#ifdef __cplusplus
}
#endif

#endif /* R_GRAPHICSENGINE_ */
