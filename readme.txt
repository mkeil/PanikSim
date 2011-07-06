Bei Aufruf von CMake ist darauf achten, dass externe Bibliotheken erreichabr sind (X11, CUDA, GTest, GLog)

Hinweise:  die Datei panic.par enthält alle wichtigen Parameter für die Simulation. Die Datei muß im selben Verzeichnis sein, wie das ausführbare Programm. 

Wenn die Partikelanzahl (N0) stark erhöht wird, sollte entsprechend die Raumgröße (RoomXSize und RoomYSize) angepaßt werden, weil sonst die Initialsierung der Partikel nicht fertig wird. Außerdem sollte der Parameter X11_Magn verkleinert werden, damit das Bild ordentlich dargestellt wird. 