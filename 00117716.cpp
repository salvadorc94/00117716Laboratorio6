#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "display_tools.h"
#include "tools.h"
#include "sel.h"
#include "assembly.h"


int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;

    cout << "METODOS DE LOS ELEMENTOS FINITOS EN DOS DIMENSIONES CON FUNCIONES DE FORMA\n"
         << "LINEALES Y PESOS DE GALERKIN - 00117716\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    //cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);

    zeroes(K,m.getSize(NODES));
    zeroes(b,m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);

    applyNeumann(m,b);

    applyDirichlet(m,K,b);

    //showMatrix(K);
    zeroes(T,b.size());
    calculate(K,b,T);

    //cout << "La respuesta es: \n";
    //showVector(T);

    writeResults(m,T,filename);

    return 0;
}
