#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

Display *display;
Window win;
GC gc;
Pixmap pix1;              /* hidden windows */
int xx,yy,wwidth,hheight;
unsigned int bborder_width;

Colormap cmap;
/* int nn; */
Visual *vis;

int screen;

#define SMALL 1
#define OK 0


void g_setcolor(int color)
{

    switch( color ) {
    case 1:
        XSetForeground(display,gc,BlackPixel(display,screen));
        break;
    case 0:
        XSetForeground(display,gc,WhitePixel(display,screen));
        break;
    default:
        XSetForeground(display,gc,color);
        break;
    }
}


void h_erase(int xsize, int ysize)
{

    g_setcolor(0);
    XFillRectangle(display,pix1,gc,0,0,xsize,ysize);
    g_setcolor(1);
}


void h_init(int xsize, int ysize)
{

    /* initialize the hidden window */
    pix1=XCreatePixmap(display,win,xsize,ysize,DefaultDepth(display,screen));
    h_erase(xsize,ysize);
}


void g_init(char* window_name, char* icon_name, int x, int y,
            int width, int height, unsigned int border_width)
{

    Pixmap icon_pixmap;
    XSizeHints size_hints;
    char *display_name=NULL;
    unsigned long valuemask =0;
    XGCValues values;

    int aargc;

    xx = x ;
    yy = y;
    wwidth = width;
    hheight = height;
    bborder_width = border_width;

    if((display=XOpenDisplay(display_name))==NULL) {
        printf("basicwin:  cannot connect to X Server %s\n",
               XDisplayName(display_name));
        exit(-1);
    }

    screen=DefaultScreen(display);

    win=XCreateSimpleWindow(display,RootWindow(display,screen),x,y,width,height,
                            border_width,BlackPixel(display,screen),
                            WhitePixel(display,screen));

    size_hints.flags=PPosition|PSize|PMinSize;
    size_hints.x=x;
    size_hints.y=y;
    size_hints.width=width;
    size_hints.height=height;
    size_hints.min_width=50;
    size_hints.min_height=50;

    aargc= 0;

    XSetStandardProperties(display,win,window_name,icon_name,icon_pixmap,NULL,aargc,&size_hints);

    XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask|StructureNotifyMask);

    gc=XCreateGC(display,win,valuemask,&values);
    XSetForeground(display,gc,BlackPixel(display,screen));
    XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinMiter);
    XSetBackground(display,gc,4);

    h_init(width,height);

    XMapWindow(display,win);
    cmap = DefaultColormap(display,screen);
}


void g_font( char* sw, char* fontname )
{

    XFontStruct* font_info;

    if(strcmp(sw,"open")==0) {
        if( !( font_info = XLoadQueryFont(display,fontname) ) ) {
            fprintf( stderr, "Error: cannot open %s font.\n", fontname);
            fprintf( stderr, "Using fixed font.\n" );
            fflush( stderr );
        } else {
            XSetFont( display, gc, font_info->fid );
        }
    } else { /* i.e. if(strcmp(sw,"close")==0) */
        XUnloadFont( display, font_info->fid );
    }

}


void g_close()
{
    XCloseDisplay(display);
    XFreePixmap(display,pix1);
    exit(1);
}



void h_show(int xsize, int ysize)
{

    XCopyArea(display,pix1,win,gc,0,0,xsize,ysize,0,0);
    XSync(display,False);
}


void g_win( char* sw, char* window_name, char* icon_name, int x, int y,
            int width, int height, unsigned int border_width )
{

    if(strcmp(sw,"open")==0) {
        g_init(window_name, icon_name, x, y, width, height, border_width );
    } else { /* i.e. if (strcmp(sw,"close")==0) */
        g_close();
    }
}


