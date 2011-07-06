#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>

#include <glog/logging.h>

#include "Xlib_mod.h"




void g_setcolor(int color, X11props_t *X11props)
{
    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    GC gc = X11props -> gc;
    int screen = X11props -> screen;

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


    // lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> gc = gc;
    X11props -> screen = screen;
}


void h_erase(int xsize, int ysize, X11props_t *X11props)
{
    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    GC gc = X11props -> gc;
    Pixmap pix1 = X11props -> pix1;


    g_setcolor(0, X11props);

    // lokale Werte geupdatet
    display = X11props -> display;
    gc = X11props -> gc;
    pix1 = X11props -> pix1;

    XFillRectangle(display,pix1,gc,0,0,xsize,ysize);

    // Prop Rec updaten = lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> gc = gc;
    X11props -> pix1 = pix1 ;

    g_setcolor(1, X11props);

}


void h_init(int xsize, int ysize, X11props_t *X11props)
{
    LOG (INFO) << "h_init gestartet.";

    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    Window win = X11props -> win;
    Pixmap pix1 = X11props -> pix1;
    int screen = X11props -> screen;
    LOG (INFO) << "Probs lokal zur Verfügung gestellt.";

    // LOG (INFO) << "Display: " << display;

    /* initialize the hidden window */
    // LOG (INFO) << "XCreatePixmap wird aufgerufen."; LOG (INFO) << "display: " << display; LOG (INFO) << "win: " << win; LOG (INFO) << "xsize: " << xsize; LOG (INFO) << "ysize: " << ysize;

    pix1=XCreatePixmap(display,win,xsize,ysize,DefaultDepth(display,screen));



    LOG (INFO) << "h_erase wird aufgerufen.";
    // Update des Prop Records, d.i. in diesem Fall das selbe wie "lokale Veränderungen speichern", weil danach an den lokalen Werten nix mehr gemacht wird
    X11props -> display = display;
    X11props -> win = win;
    X11props -> pix1 = pix1 ;
    X11props -> screen = screen;

    h_erase(xsize,ysize, X11props);
}


void g_init(char* window_name, char* icon_name, int x, int y,
            int width, int height, unsigned int border_width,
            X11props_t *X11props)
{
    LOG (INFO) << "g_init gestartet.";
    // benötigte Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    Window win = X11props -> win;
    GC gc = X11props -> gc;
    int xx,yy,wwidth,hheight;
    xx = X11props -> xx;
    yy = X11props -> yy;
    wwidth = X11props -> wwidth;
    hheight = X11props -> hheight;
    unsigned int bborder_width = X11props -> bborder_width;
    Colormap cmap = X11props -> cmap;
    int screen = X11props -> screen;
    LOG (INFO) << "Probs lokal zuer Verfügung gestellt.";

    // lokale Variablen
    Pixmap icon_pixmap;
    XSizeHints size_hints;
    char *display_name=NULL;
    unsigned long valuemask =0;
    XGCValues values;
    int aargc;
    LOG (INFO) << "lokale Variablen deklariert.";



    xx = x ;
    yy = y;
    wwidth = width;
    hheight = height;
    bborder_width = border_width;
    LOG (INFO) << "lokale Werte gesetzt.";

    if((display=XOpenDisplay(display_name))==NULL) {
        printf("basicwin:  cannot connect to X Server %s\n",
               XDisplayName(display_name));
        exit(-1);
    }
    LOG (INFO) << "X11 Server erreichbar. Display: " << display;

    screen=DefaultScreen(display);
    LOG (INFO) << "screen gesetzt: " << screen;

    win=XCreateSimpleWindow(display,RootWindow(display,screen),x,y,width,height,
                            border_width,BlackPixel(display,screen),
                            WhitePixel(display,screen));

    LOG (INFO) << "win gesetzt: " << win;

    size_hints.flags=PPosition|PSize|PMinSize;
    size_hints.x=x;
    size_hints.y=y;
    size_hints.width=width;
    size_hints.height=height;
    size_hints.min_width=50;
    size_hints.min_height=50;
    LOG (INFO) << "size_hints gesetzt.";

    aargc= 0;

    XSetStandardProperties(display,win,window_name,icon_name,icon_pixmap,NULL,aargc,&size_hints);

    XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask|StructureNotifyMask);

    gc=XCreateGC(display,win,valuemask,&values);
    XSetForeground(display,gc,BlackPixel(display,screen));
    XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinMiter);
    XSetBackground(display,gc,4);

    LOG (INFO) << "h_init wird aufgerufen";
    // Prop Rec aktualisieren;
    X11props -> display = display;
    X11props -> win = win;
    X11props -> gc = gc;
    X11props -> xx = xx;
    X11props -> yy = yy;
    X11props -> wwidth = wwidth;
    X11props -> hheight = hheight;
    X11props -> bborder_width = bborder_width;
    X11props -> cmap = cmap ;
    X11props -> screen = screen;

    h_init(width,height, X11props);

    // lokalen Werte aktualisieren
    display = X11props -> display;
    win = X11props -> win;
    gc = X11props -> gc;
    xx = X11props -> xx;
    yy = X11props -> yy;
    wwidth = X11props -> wwidth;
    hheight = X11props -> hheight;
    bborder_width = X11props -> bborder_width;
    cmap = X11props -> cmap;
    screen = X11props -> screen;

    LOG (INFO) << "h_init beendet.";

    XMapWindow(display,win);
    cmap = DefaultColormap(display,screen);

    LOG (INFO) << "funktionsrumpf beendet.";

    // lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> win = win;
    X11props -> gc = gc;
    X11props -> xx = xx;
    X11props -> yy = yy;
    X11props -> wwidth = wwidth;
    X11props -> hheight = hheight;
    X11props -> bborder_width = bborder_width;
    X11props -> cmap = cmap ;
    X11props -> screen = screen;
}


void g_font( char* sw, char* fontname, X11props_t *X11props)
{
    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    GC gc = X11props -> gc;

    // lokale Variablen
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

    // lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> gc = gc;
}


void g_close(X11props_t *X11props)
{
    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    Pixmap pix1 = X11props -> pix1;
    XCloseDisplay(display);
    XFreePixmap(display,pix1);
    exit(1);

    // lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> pix1 = pix1 ;

}



void h_show(int xsize, int ysize, X11props_t *X11props)
{
    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    Window win = X11props -> win;
    GC gc = X11props -> gc;
    Pixmap pix1 = X11props -> pix1;

    XCopyArea(display,pix1,win,gc,0,0,xsize,ysize,0,0);
    XSync(display,False);

    // lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> win = win;
    X11props -> gc = gc;
    X11props -> pix1 = pix1 ;
}


void g_win( char* sw, char* window_name, char* icon_name, int x, int y,
            int width, int height, unsigned int border_width, X11props_t *X11props )
{
    LOG (INFO) << "g_win gestartet.";

    if(strcmp(sw,"open")==0) {
        LOG (INFO) << "g_init wird aufgerufen.";
        g_init(window_name, icon_name, x, y, width, height, border_width,X11props );

    } else { /* i.e. if (strcmp(sw,"close")==0) */
        g_close(X11props);
    }

}


