Installationshinweis: ****
- 1.) CMake starten und darauf achten, dass alle nötigen Bibliotheken erreichbar sind (X11, CUDA, GLOG). Die CMakeLists.txt befindet sich im Hauptverzeichnis.
- 2.) Anschließend das makefile starten. 

Das Makefile erzeugt zwei ausführbare Programme. Die Version mit fester Stepsize läuft schneller ist entspricht aber nicht dem originalem Programm. 

Einschränkungen: ****
- Das Programm wurde nur für CUDA Karten mit eine CC 2.0 getestet. 
	
Hinweise: ****
Die Datei panic.par enthält alle wichtigen Parameter für die Simulation. Die Datei muß im selben Verzeichnis sein, wie das ausführbare Programm sein. 
Wenn die Partikelanzahl (N0) stark erhöht wird, sollte entsprechend die Raumgröße (RoomXSize und RoomYSize) angepaßt werden, weil sonst die Initialsierung der Partikel nicht fertig wird. Außerdem sollte der Parameter X11_Magn verkleinert werden, damit das Bild ordentlich dargestellt wird. 

Um lange Simulationen zu erleben, muß der Parameter MaxUpdNum erhöht werden. 

Um die Säule ein- bzw. auszuschalten dient der Parameter Colum Switch. 

Ausführliche Beschreibungen zu den Parametern befindet sich im Ordner doku\Parameter.pdf

Im ordner Samples befinden sich verschiedene Beispiele. Diese Datein ins Arbeitsverzeichnis kopieren und damit die orginale panic.par ersetzten. 

theoretische Grundlagen: ****

Die Simluation basiert auf folgender Arbeit

Dirk Helbing, Illes J. Farkas, and Tamas Vicsek
Simulating dynamical features of escape panic
Nature 407, 487-490 (2000)


Inhalt des Verzeichnisses: ****
\.git
\cmake
\CMakeLists.txt	
\doku
\panic.par	- Parameter, die Einfluß auf die Simulation haben
\readme.txt	- diese Datei hier
\src


\doku\PanicPackage.tar	- orginale sequentielle Implementierung des Programms
\doku\Parameter.pdf		- Beschreibung der Parameter in Panic.par
\doku\Simulating dynamical features of escape panic.pdf	- theoretischer Artikel, der die Grundlage des Programms bildet
\doku\Variablen im Hauptprogramm.pdf	- 

\src\base.c			- definiert grundlegende Funktione
\src\deviceFunc.cu	- enthält die Devicefunktionen
\src\deviceFunc.h
\src\hostFunc.cu	- enthält die nötigen Funktionen für den Host
\src\hostFunc.h
\src\kernels.cu		- enthält alle Kernels
\src\kernels_constant_stepsize.cu		- enthält alle Kernels, hier wird die feste StepSize verwendet
\src\kernels.h
\src\PanicSimMain.cu	- das Hauptprogramm
\src\types.h		
\src\Xlib_mod.cu		- Funktionen für die Kommunikation mit der X11 Bibliothek
\src\Xlib_mod.h
