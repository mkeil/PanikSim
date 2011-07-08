#pragma once 
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>



#include "types.h"

#define SMALL 1
#define OK 0


void g_setcolor(int color, X11props_t *X11props); 

void h_erase(int xsize, int ysize, X11props_t *X11props) ;
void h_init(int xsize, int ysize, X11props_t *X11props) ;

void g_init(char* window_name, char* icon_name, int x, int y,
            int width, int height, unsigned int border_width, 
			X11props_t *X11props); 

void g_font( char* sw, char* fontname, X11props_t *X11props) ;

void g_close(X11props_t *X11props);

void h_show(int xsize, int ysize, X11props_t *X11props) ;

void g_win( char* sw, char* window_name, char* icon_name, int x, int y,
            int width, int height, unsigned int border_width, X11props_t *X11props );
			