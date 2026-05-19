/* Datei:  tdse.c
   Author: BCM
   Datum:  27.06.2009 
   Datum:  11.06.2012
   Datum:  21.06.2016 weiteren Ausgabemodus hinzugefuegt
   Datum:  2025 Aktualisiert auf C++ und HDF5 (TL)
*/

/* Berechnung der Zeitentwicklung einer Wellenfunktion in einem 1-dimensionalen
   Oszillatorpotential durch Loesung der zeitabhaengigen Schroedingergleichung
   mit der Crank-Nicholson Methode.
*/

/* Zunaechst Definitionen der Standardbibliotheken einbinden mit #include */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>    /* mathematische Funktionen */ 
#include <complex>  /* Komplexe Arithmetik */
#include <vector>

/* Alle Daten werden in einer H5-Datei speichern */
#include <highfive/H5File.hpp>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

/* Vorwaertsdeklaration der Funktionen: */

void GetUserParam( int, char *[] );            /* Eingabe ueber die Kommandozeile: */
void InitPsi( std::vector<std::complex<double>> & );              /* Anfangswellenfunktion */
void InitPot( std::vector<double> &);                      /* Potential */
void InitCN( long, std::vector<double> &, std::vector<std::complex<double>> &, std::vector<std::complex<double>> & ); 
                                               /* Initialisierung der Feldvariablen 
                                                  fuer CN tridiagonalmatrix */
void CN( long, std::vector<std::complex<double>> &,               /* berechne CN Zeitschritt */
         std::vector<std::complex<double>> &, std::vector<std::complex<double>> &, 
         std::vector<std::complex<double>> &, std::vector<std::complex<double>> & ); 
std::complex<double> Energie( std::vector<std::complex<double>> &, std::vector<double> );            /* "Energie"-Erwartungswert */
double Norm( std::vector<std::complex<double>> & );               /* Normierung */

/* Globale Variablen ( koennen mit GetUserParam geaendert werden ): */

double x_mn   =  -10.0;    /* Linker Rand '-L' */
double x_mx   =   10.0;    /* Rechter Rand '-R' */
double h      =   0.05;    /* Ortschrittweite '-h' */
int n_x;                   /* # Ortspunkte */ 
int n_t       =  5000;    /* # Zeitschritte '-T' */ 
double d      =  0.001;    /* Zeitschrittweite '-d' */
int D         =     25;    /* Periode zur Ausgabe '-D' */
double b      =    0.75;    /* Breite der Anfangswellenfunktion '-b' */
double x_0    =   0.0;    /* Mittlere Position der Anfangswellenfunktion '-x' */
double k_0    =    3.0;    /* Mittlere Impuls der Anfangswellenfunktion '-k' */

const std::complex<double> I (0.0,1.0);  

/* ======================= MAIN ===========================*/

int main( int argc, char *argv[] ) {
 
    long dim;                             /* zur Speichplatzreservierung der Felder */
    std::vector<std::complex<double>> Psi;
    std::vector<std::complex<double>> p,q,alpha,beta;
    std::vector<double> V;                            /* Potential */
    std::vector<double> xn;     /* Gitter */
    int n, i;                             /* Indizes Gitterpunkte t=n.d; x=x_mn+i.h */
    double t, x;                          /* Zeit / Ort */
    std::complex<double> Eval;                  /* Energieerwartungswert */
    double Nval;                          /* Normalisierung */
    std::vector<int> shape;

    std::string h5file="tdse.h5"; // Name der H5-Datei

    // Öffnen die H5 file
    //HighFive::File h5(h5file, HighFive::File::ReadWrite | HighFive::File::OpenOrCreate);
    HighFive::File h5(h5file, HighFive::File::ReadWrite | HighFive::File::Overwrite);
    HighFive::DataSet dataset; // Container zum Speichern von Dingen in einer H5-Datei
    
    GetUserParam( argc, argv ); /* Parameter ueber die Kommandozeile */

    /* Zeitschrittweite in Einheiten der Periode: */
    d *= 2 * M_PI; /* M_PI wird von math.h bereitgestellt */
    
    n_x = ( x_mx - x_mn + 0.1 * h ) / h; /* # Gitterpunkte x */
    dim = n_x + 1;                       /* # Feldvariable */

    /* Ausgabe der Parameter: */
    std::printf( "# tdse:\n" );
    std::printf( "# ==== \n" );
    std::printf( "# Gitter: x_min=%f; x_max=%f; %i Punkte; %i Zeitpunkte\n", 
            x_mn, x_mx, n_x, n_t );
    std::printf( "# Schrittweiten: dx=%f; dt=%f\n", h, d );
    std::printf( "# Psi: b= %f; x_0=%f; k_0=%f\n", b, x_0, k_0 );

    Psi.resize(dim);    /* Speicherplatzreservierung Psi */
    beta.resize(dim);   /* Diagonale der tridiagonalen Matrix */
    alpha.resize(dim);  /* Nebendiagonale der ... */
    p.resize(dim);      /* Arbeitsspeicherplatzreservierung */
    q.resize(dim);      /* Arbeitsspeicherplatzreservierung */
    V.resize(dim);      /* Speicherplatzreservierung V */
    xn.resize(dim);     /* Gitterpunkte */
    
    InitPsi( Psi ); /* Initialisiere Psi */
    InitPot( V );   /* Initialisiere V */
    for ( i=0; i<=n_x; i++ )  xn[i] = x_mn + i * h; /* mach Gitter */
    InitCN( n_x, V, alpha, beta ); /* Initialisiere Vektoren fuer CN */

    dataset = h5.createDataSet<double>("Gitter",  HighFive::DataSpace::From(xn));
    dataset.write(xn);
    dataset = h5.createDataSet<double>("Potential",  HighFive::DataSpace::From(V));
    dataset.write(V);

    /* Ausgabe zu t=0 */
    t = 0;
    n=0;
    Eval = Energie( Psi, V); /* Berechne Energieerwartungswert */
    Nval = Norm( Psi );      /* Berechne Normierung */
    dataset = h5.createDataSet<std::complex<double>>(std::to_string(n)+"/Energie",  HighFive::DataSpace::From(Eval));
    dataset.write(Eval);
    dataset = h5.createDataSet<double>(std::to_string(n)+"/Normierung",  HighFive::DataSpace::From(Nval));
    dataset.write(Nval);
    dataset = h5.createDataSet<double>(std::to_string(n)+"/Zeit",  HighFive::DataSpace::From(t));
    dataset.write(t);
    dataset = h5.createDataSet<std::complex<double>>(std::to_string(n)+"/Psi",  HighFive::DataSpace::From(Psi));
    dataset.write(Psi);
    
    for (n=1; n<=n_t; n++) {
      t = n * d / 2 / M_PI; /* Zeit in Einheiten der Periode */
      CN( n_x, Psi, alpha, beta, p, q ); /* Zeitschritt nach Crank-Nicholson */
      Eval = Energie( Psi, V); /* Berechne Energieerwartungswert */
      Nval = Norm( Psi );      /* Berechne Normierung */
      dataset = h5.createDataSet<std::complex<double>>(std::to_string(n)+"/Energie",  HighFive::DataSpace::From(Eval));
      dataset.write(Eval);
      dataset = h5.createDataSet<double>(std::to_string(n)+"/Normierung",  HighFive::DataSpace::From(Nval));
      dataset.write(Nval);
      dataset = h5.createDataSet<double>(std::to_string(n)+"/Zeit",  HighFive::DataSpace::From(t));
      dataset.write(t);
      dataset = h5.createDataSet<std::complex<double>>(std::to_string(n)+"/Psi",  HighFive::DataSpace::From(Psi));
      dataset.write(Psi);
      if ( (n%D)==0 ) { /* Ausgabe am jeden D-ten Zeitschritt */
	std::printf( "%10.4f %20.8e%20.8e%20.8e\n", t, 
		     Nval, std::real( Eval ), std::imag( Eval ) );
      }
    }

    return 0;
}

/* =========================== INITPSI ==================================== */
void InitPsi( std::vector<std::complex<double>> & Psi ) {
    /* Initialisiere die Anfangswellenfuntion:
       Psi(x) = N exp(-(x-x_0)^2/(2 b^2)) exp(i k_0 x)
    */

    int i; /* Gitterindex */
    double x, norm, gauss;

    norm = 1.0 / sqrt( b * sqrt(M_PI) ); /* M_PI wird von math.h bereitgestellt */

    for (i=0; i<=n_x; i++) {
        x = x_mn + i * h;
        gauss = norm * exp( -0.5 * ( x - x_0 ) * ( x - x_0 ) / b / b );
        Psi[i] = std::exp( (k_0 * x) * I ) * gauss; /* I imag. Einheit aus complex.h */
    }

}
/* =========================== InitPot ==================================== */
void InitPot( std::vector<double> & V ) {

    int i; /* Gitterindex */
    double x;

    for (i=0; i<=n_x; i++) {
        x = x_mn + i * h;
        V[i] = 0.5 * x * x; /* Oszillatorpotential */
    }

}
/* =========================== InitCN ==================================== */
void InitCN( long N, std::vector<double> & V, std::vector<std::complex<double>> &alpha, std::vector<std::complex<double>> & beta ) {
    /* Initialisiere Diagonalelemente beta und Nebendiagonalelemente alpha
       fuer die Tridiagonalmatrix im CN-Verfahren */

    int i;
    double h22d, h22;

    h22 = 2 * h * h;
    h22d = h22 / d;

    /* berechne alpha(i)=0, i=1,...,N */
    alpha[0] = (std::complex<double>) 0.0;
    for ( i=1; i<N; i++ ) {
        alpha[i] = (std::complex<double>) 1.0;
    }
    alpha[N] = (std::complex<double>) 0.0;

    /* berechne beta(i) */
    beta[0] = (std::complex<double>) 2.0 * h22d * I;
    for ( i=1; i<N; i++ ) {
        beta[i] = (std::complex<double>) 2.0 * h22d * I - 2.0 - h22 * V[i]; 
    }
    beta[N] = (std::complex<double>) 2.0 * h22d * I;

}

/* =========================== CN ==================================== */
void CN( long N, 
         std::vector<std::complex<double>> &Psi, 
         std::vector<std::complex<double>> &alpha, 
         std::vector<std::complex<double>> &beta, 
         std::vector<std::complex<double>> &p, 
         std::vector<std::complex<double>> &q ) {
    /* loese das tridiagonale Gleichungssystem im CN-Schritt */

    int i;
    std::complex<double> b, chi;

    b = (std::complex<double>) 8.0 * h * h / d * I; /* Faktor fuer rechte Seite */

    /* Ruckwaertssubstitution */

    p[N-1] = - alpha[N] / beta[N];
    q[N-1] = b * Psi[N] / beta[N];

    for ( i=N-1; i>0; i-- ) {
        p[i-1] = - alpha[i] / ( beta[i] + alpha[i] * p[i] );
        q[i-1] = ( b * Psi[i] - alpha[i] * q[i] ) / ( beta[i] + alpha[i] * p[i] );
    }
    /* Forwaertssubstitution */
    
    chi = ( b * Psi[0] - alpha[0] * q[0] ) /
          ( beta[0] + alpha[0] * p[0] );

    Psi[0] = chi - Psi[0];

    for ( i=0; i<N; i++ ) {
        chi = p[i] * chi + q[i];
        Psi[i+1] = chi - Psi[i+1];
    }
    
}

/* =========================== Energie ==================================== */
std::complex<double> Energie( std::vector<std::complex<double>> & Psi, std::vector<double> V ) {
    /* Berechne Erwartungswert < Psi, H Psi > */
    
    std::complex<double> d2psi_dx2;
    std::complex<double> term, sum;
    int i;
/* berechne                              
  R  /     2                      \        
 /   |  1 d                       |------ ~
 | dx|- - --- Psi(x) + V(x) Psi(x)|Psi(x) =
 /   |  2   2                     |       
L    \    dx                      /       
n_x-1
 -- /                                              \
 \  |  1 Psi(i-1)-2 Psi(i)+Psi(i+1)                |------
 /  |- - -------------------------- + h V(i) Psi(i)|Psi(i)
 -- |  2          h                                |
i=1 \                                              /
*/
    sum = (std::complex<double>) 0.0;
    for ( i=1; i<n_x; i++ ) {
        term = h* V[i] * Psi[i];
        d2psi_dx2 = Psi[i-1] - 2.0 * Psi[i] + Psi[i+1];
        d2psi_dx2 *= -0.5 / h;
        term += d2psi_dx2;
        term *= conj( Psi[i] );
        sum += term;
    }

    return sum;
}

/* =========================== Norm ==================================== */
double Norm( std::vector<std::complex<double>> & Psi ) {
    
    double sum;
    int i;
/* berechne     n_x-1
  R               -- 
 /            2 ~ \            2 
 | dx |Psi(x)|  = /  h |Psi(i)| 
 /                --
L                i=1
*/
    sum = 0.0;
    for ( i=1; i<n_x; i++ ) {
        /* berechne ( * h ! )
                2
           |Psi|
        */
      sum += h * std::abs( std::conj(Psi[i]) * Psi[i] );
    }
    return sum;
}

/* =========================== GETUSERPARAM ==================================== */

void GetUserParam( int argc, char *argv[] ){

/* Variablen / Konstanten: */
    int i;
    char* endptr;
    const char* usage = 
        "tdse [-b <Breite> -x <Position> -k <Impuls>\
 -d <delta t> -h <delta x>  -L <x_min> -R <x_max> -T <# t> -D <Ausgabe nach o Schritte>]";
    const char* error_message =
        "# FEHLER(GetuserParam): falsche Option: ";

    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
        /* Es gibt immer mindestens 1 Parameter: argv[0] = Kommandoname */
        for (i=1; i<argc; i++){
            /* Parameter sind 2 Charaktere lang und sollten mit '-' anfaengen ... */
            if ( (strlen(argv[i])==2) && (argv[i][0] == '-') ) { 
                switch (argv[i][1]) { 
                    case 'b':
                        b = strtod( argv[++i], &endptr); /* String --> double */
                        /* ... sollte mit Lehrzeichen beendet werden ... */
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'x':
                        x_0 = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'k':
                        k_0 = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'd':
                        d = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'h':
                        h = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'L':
                        x_mn = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'R':
                        x_mx = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'T':
                        n_t = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    case 'D':
                        D = (int) strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
                            std::printf(error_message);
                            std::printf(usage);
                            std::printf("\n");
                            exit(1);
                        }
                        break;
                    // case 'o':
                    //     o = (int) strtol( argv[++i], &endptr, 10);
                    //     if ( (!isspace(*endptr) && (*endptr) != 0) ) {
                    //         std::printf(error_message);
                    //         std::printf(usage);
                    //         std::printf("\n");
                    //         exit(1);
                    //     }
                    //     break;
                    default:
                        std::printf(error_message);
                        std::printf(usage);
                        std::printf("\n");
                        exit(1);
                }
            } else {
                std::printf(error_message);
                std::printf(usage);
                std::printf("\n");
                exit(1);
            } /* end of: if-then-else */
        } /* end-of: for */ 
    } /* end-of: if */

}
/* Kompilieren mit:
gcc -l m tdse.c -o tdse

# tdse:
# ====
# Gitter: x_min=-10.000000; x_max=10.000000; 200 Punkte; 100 Zeitpunkte
# Schrittweiten: dx=0.100000; dt=0.628319
# Psi: b= 1.000000; x_0=0.000000; k_0=0.000000
    1.0000       1.00000000e+00      4.99687760e-01     -6.60892954e-20
    2.0000       1.00000000e+00      4.99687760e-01     -2.34461932e-18
    3.0000       1.00000000e+00      4.99687760e-01     -1.32199192e-17
    4.0000       1.00000000e+00      4.99687760e-01     -1.23670045e-17
    5.0000       1.00000000e+00      4.99687760e-01     -4.87715459e-18
    6.0000       1.00000000e+00      4.99687760e-01      1.15196983e-17
    7.0000       1.00000000e+00      4.99687760e-01      2.28989358e-17
    8.0000       1.00000000e+00      4.99687760e-01      3.10498119e-17
    9.0000       1.00000000e+00      4.99687760e-01     -4.12563991e-17
   10.0000       1.00000000e+00      4.99687760e-01     -7.52125367e-17
*/
