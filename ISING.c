#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define FicheroInit "Datos/init.txt" //Fichero desde el que se inicializan las variables
#define FicheroOutput "Datos/output.dat" //Ruta de destino en guardaConfiguracion()
#define NormRANu (2.3283063671E-10F)

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

int L, V; //Lado red NOTA: V se inicializa en configuracionInicial()
int semilla_rand; //semilla para inicializar la funci�n Random
int n_pasos_term, n_cambios, n_medidas; //pasos de la termalizacion, cuantos cambios por medida, n� de medidas.
double b_ini, b_f, delta_b;  //valores de beta inicial, final e incremento.
int n_bloques; //n�mero de bloques que se toman para realizar el an�lisis de errores
char rutadest[30]; //Ruta de valores de output

int *xp, *yp, *xm, *ym; //vectores de desplazamiento
int *s; //vector de espines
double prob[5]; //Vector donde se guardan las probabilidades para cada beta

double *medida_magnet; //Vector donde se guardan las medidas de magnetizaci�n.
double *medida_ener; //Vector donde se guardan las medidas de energ�a

double **resultados; //Vector de resultados finales con tama�o [7][N]

double *calor_especifico;
double *susceptibilidad;

int modo_de_ejecucion;

double beta_a_analizar; //Usada en el analisis, no en la simulaci�n

//###########################################\\
//       FUNCIONES DE INICIALIZACI�N         \\
//###########################################\\


float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);

    ind_ran++;
    r=ir1*NormRANu;

    return r;

}

void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;



    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;

    }
    ind_ran=ig1=ig2=ig3=0;
}





double media(double *serie, int Numero) {
    int i,j;
    double x=0,y=0,z=0;
    for(i=0;i<Numero;i++) {
        x+=serie[i]; }

    return x/Numero;
}

 double cuadradosmedia (double *serie, int Numero) {

     double y;
     int j;
     for(j=0;j<Numero;j++)
        y+=(serie[j]*serie[j]);

return y/Numero; }

double varianza(double *serie, int Numero){

    double mediacuadrados, Media, z;
    mediacuadrados = cuadradosmedia(serie,Numero);

    Media=media(serie,Numero);
     z=(mediacuadrados-(Media*Media));
     return z;
}

double error_media_bloques(double *datos,int n_datos,int n_intervalos){

    //declaraci�n de variables para su uso.
    int interv_size=n_datos/n_intervalos,i,j;
    double datosauxiliar[interv_size],mediasint[n_intervalos];
    double varianza_int;

    for(j=0; j<n_intervalos; j++){
        // En cada iteraci�n del primer for se almacena un intervalo en el vector auxiliar, y se almacena su media en el vector de medias.

        for(i=0; i<interv_size; i++)
            datosauxiliar[i] = datos[j*interv_size+i];

        mediasint[j] = media(datosauxiliar, interv_size);
    }

    varianza_int=varianza(mediasint, n_intervalos);
    // A la funci�n que calcule la varianza le paso el vector con las medias de los intervalos y la cantida de medias calculadas(tantas como intervalos haya)

    return sqrt(varianza_int/interv_size);

}



double vmax (double *v, int N) {

    double b=v[0];
    int i;

    for (i=1; i<N; i++) {
        if (v[i]>b) {b=v[i];}
        }

    return b;}

double vmin (double *v, int N) {

    double a=v[0];
    int i;

    for (i=1; i<N; i++) {
        if (v[i]<a) {a=v[i];}
        }

    return a;}

void guardaHist (double *H, int n, double a, double delta, const char* nombre) {

FILE *f;
f=fopen(nombre,"w");
int i;
for (i=0; i<n; i++) {fprintf (f, "%f\t%f\n", a+i*delta+delta/2, H[i]);}
fclose(f);

}

void histograma (double *datos, int N, int n, double a, double b, const char* nombre) {

    int i=0;
    double delta=(b-a)/n;
    int k;
    double* H = (double*) malloc(N*sizeof(double));

    for (i=0; i<N; i++) {
    if (datos[i]==b) {H[n-1]++;}
    else {
    k=floor((datos[i]-a)/delta);
    H[k]++;}
    }

//Normaliza:

    for (i=0; i<n; i++) {
    H[i]=H[i]/(delta*N);}

    guardaHist(H, n, a, delta, nombre);
}

//Inicializa variables seg�n el fichero de inicializaci�n
void leeDatos() {


FILE *file_input;
char fname[10240];
  extern  double b_ini, b_f, delta_b;

    //lee las variables

      sprintf(fname,"caseI/configureI.input");
      file_input = fopen(fname,"r"); // "r"= only read

      fscanf(file_input,"%*s %d" ,&L);
      fscanf(file_input,"%*s %d" ,&semilla_rand);
      fscanf(file_input,"%*s %lf" ,&b_ini);
      fscanf(file_input,"%*s %lf" ,&b_f);
      fscanf(file_input,"%*s %lf" ,&delta_b);
      fscanf(file_input,"%*s %d" ,&n_pasos_term);
      fscanf(file_input,"%*s %d" ,&n_cambios);
      fscanf(file_input,"%*s %d" ,&n_medidas);
      fscanf(file_input,"%*s %d" ,&n_bloques);
      fscanf(file_input,"%*s %d" ,&modo_de_ejecucion);
      fscanf(file_input,"%*s %lf" ,&beta_a_analizar);
      fclose(file_input);



    if(beta_a_analizar < 0) beta_a_analizar = 0.5*log(1+sqrt(2)); //beta cr�tica

    printf("Datos leidos de Fichero:\n");
    printf("L = %d \nsemilla_rand = %d \nb_ini = %lf \nb_f = %lf \nn_delta = %lf \nn_pasos_term = %d \nn_cambios = %d \nn_medidas = %d \nn_bloques = %d ",
           L, semilla_rand, b_ini, b_f, delta_b, n_pasos_term, n_cambios, n_medidas, n_bloques);
}

//Inicializa los vectores de desplazamiento
void inicializaVectoresMovimiento() {

    int i;

    xp = (int*) malloc(L*sizeof(int));
    if(xp==NULL) { printf("Error al reservar memoria)"); exit(2); }

    yp = (int*) malloc(L*sizeof(int));
    if(xp==NULL) { printf("Error al reservar memoria)"); exit(3); }

    xm = (int*) malloc(L*sizeof(int));
    if(xp==NULL) { printf("Error al reservar memoria)"); exit(4); }

    ym = (int*) malloc(L*sizeof(int));
    if(xp==NULL) { printf("Error al reservar memoria)"); exit(5); }

    *(xp + L-1) = -(L-1);
    *(yp + L-1) = -L*(L-1);
    *(xm) = L-1;
    *(ym) = L*(L-1);

    for(i=1; i<L; i++)
    {
        *(xp + i-1) = 1;
        *(yp + i-1) = L;
        *(xm + i) = -1;
        *(ym + i) = -L;
    }
}

//Genera una configuraci�n aleatoria (convendr�a revisarla para que haga otras cosas seg�n un flag)
void configuracionInicial() {

    int i;
    V = L*L;

    s = (int*) malloc(V*sizeof(int)); //Reserva espacio para V y lo guarda en un puntero int.
    if(s==NULL) { printf("Error al reservar memoria)"); exit(6); }

    for(i=0; i<V; i++) {
        if(Random() > 0.5) *(s+i) = 1;
        else *(s+i) = -1;
    }
}

//Reserva memoria para los vectores medida_magnet y medida_ener
void inicializaMedidas()
{
    medida_magnet = (double *) malloc(n_medidas*sizeof(double));
    if(s==NULL) { printf("Error al reservar memoria)"); exit(7); }

    medida_ener = (double *) malloc(n_medidas*sizeof(double));
    if(s==NULL) { printf("Error al reservar memoria)"); exit(8); }
    for(int i=0;i<4;i++)
    {
        //printf( "%lf", medida_ener[99]);
        //printf("\n");
        // printf("%d",n_medidas);
        // printf("\n");

    }

}

//Reserva memoria para el vector resultados. N es el n�mero de par�metros a guardar
void inicializaVectorResultados(int N, int n_resultados)
{
    int i;
    resultados = (double **) malloc(N*sizeof(double*));
    if(resultados==NULL) { printf("Error al reservar memoria)"); exit(9); }

    for(i=0; i<N; i++)    {

        resultados[i] = (double*) malloc(n_resultados*sizeof(double));
        if(resultados[i]==NULL) { printf("Error al reservar memoria)"); exit(10); }
    }

    calor_especifico = (double *) malloc(n_resultados*sizeof(double));
    susceptibilidad = (double *) malloc(n_resultados*sizeof(double));
}



//Calcula la magnetizaci�n de la red
double calcularMagnetizacion()
{
    double magnetizacion = 0;
    int n;

    for(n=0; n<V; n++)
        magnetizacion += s[n];

    magnetizacion /= V;

    return magnetizacion;
}

//Calcula la energ�a de la red
double calcularEnergia()
{
    double energia = 0;
    int i,j;
    int n = 0;

    for(j=0; j<L; j++)
        for(i=0; i<L; i++)
        {
            energia += s[n]*(s[n+xp[i]] + s[n+yp[j]]);
            n++;
        }

    energia /= -2*V; //EL 2 REPRESENTA LA DIMENSION

    return energia;
}

/*Realiza las medidas de la magnetizaci�n y la energ�a de la red y las guarda en la posicion numero_medida de los respectivos vectores */

void mide(int numero_medida)
{
    medida_magnet[numero_medida] = calcularMagnetizacion();
    medida_ener[numero_medida] = calcularEnergia();
}

//Calcula las probabilidades posibles para una beta dada (guardadas en el vector prob)
void calcula_prob(double beta) {

    prob[0]=exp(-beta*(-8));
    prob[1]=exp(-beta*(-4));
    prob[2]=exp(-beta*0);
    prob[3]=exp(-beta*4);
    prob[4]=exp(-beta*8);
}

//Eval�a toda la red y la cambia seg�n su probabilidad (metropolis)
void cambiaRed()
{
    int ind, n=0, i, j;

    for(i=0; i<L; i++)
        for(j=0; j<L; j++)
        {
            ind = s[n]*(s[n+xp[j]]+s[n+yp[i]]+s[n+xm[j]]+s[n+ym[i]])/2+2; //Calcula el �ndice del vector prob que corresponde

            if(Random()<prob[ind]) //Eval�a probabilidad de cambiar el spin
                s[n] = -s[n];

              n++;
        }
}

void Calcula_valores_medios(double *media_E, double *media_m, double *med_cuad_E, double *med_cuad_m,double *error_e,double *error_m){

    int i, n_intervalos=100;
    double energia,magnetizacion = 0,energiacuadrado,magnecuadrado;

    *media_E=media(medida_ener,n_medidas );

     for(i=0; i<n_medidas ; i++)
         magnetizacion += fabs(medida_magnet[i]);



    *media_m= magnetizacion/n_medidas ;
    *med_cuad_E = cuadradosmedia(medida_ener,n_medidas );
    *med_cuad_m = cuadradosmedia(medida_magnet,n_medidas );

    *error_e=error_media_bloques(medida_ener, n_medidas , n_intervalos);
    *error_m=error_media_bloques(medida_magnet,n_medidas ,n_intervalos);
}

void construyeVectorMedidas(int n_beta, double beta)
{
    double media_E, media_m, med_cuad_E, med_cuad_m, error_e, error_m;

 Calcula_valores_medios (&media_E, &media_m, &med_cuad_E, &med_cuad_m, &error_e, &error_m);

    resultados[0][n_beta] = beta;
    resultados[1][n_beta] = media_E;
    resultados[2][n_beta] = error_e;
    resultados[3][n_beta] = media_m;
    resultados[4][n_beta] = error_m;
    resultados[5][n_beta] = med_cuad_E;
    resultados[6][n_beta] = med_cuad_m;

    //resultados[7][n_beta] = 2*V*(resultados[5][n_beta]-(resultados[1][n_beta]*resultados[1][n_beta]));
    //resultados[8][n_beta] = V*(resultados[6][n_beta]-(resultados[3][n_beta]*resultados[3][n_beta]));
}

void escribeMedidas(int n_beta)
{
    int i;
    for(i=0; i<7; i++)
        printf("%g\t", resultados[i][n_beta]);
        printf("\n");
}

/*void escribeMedidasFichero(int n_beta)
{


    int i;
    printf("me vas a escrbir");
        for(i=0; i<7; i++)fprintf(fp," %14.14lf  ",resultados[i][n_beta]);

        fprintf(fp,"%14.14lf %14.14lf",calor_especifico[n_beta], susceptibilidad[n_beta]);

        fprintf(fp,"\n");


    fclose(fp);
}

*/

void iteraBetas(double b_ini, double b_f)
{

    int n_beta, N_pasos; //N�mero de beta (para el �ndice de los vectores de medidas), N�mero de pasos de beta a realizar
    int n_pasos_term_hechos, n_cambios_hechos, n_medidas_hechas; //Controlan cuantas veces se ha hecho cada cosa en los respectivos bucles
    double beta = b_ini, incremento_b;  //beta que se analiza e incremento de b teniendo en cuenta el signo

    if(b_ini < b_f) incremento_b = delta_b; //Sentido creciente
    else incremento_b = -delta_b;          //Sentido decreciente

    int i;
    FILE *fp;
    fp=fopen("outputI-files/listISING.out","w");

    for(n_beta = 0, N_pasos = fabs((b_f - b_ini)/delta_b); n_beta < N_pasos; n_beta++, beta+=incremento_b)
    {
        calcula_prob(beta); //Construye vector de probabilidades


        for(n_pasos_term_hechos = 0; n_pasos_term_hechos<n_pasos_term; n_pasos_term_hechos++) //Termaliza
            cambiaRed();

        for(n_medidas_hechas = 0; n_medidas_hechas < n_medidas; n_medidas_hechas++)
        {
            for(n_cambios_hechos = 0; n_cambios_hechos < n_cambios; n_cambios_hechos++) //Cambia la red sin medir
                cambiaRed();

            mide(n_medidas_hechas); //Toma medidas y las guarda en los vectores
        }
        construyeVectorMedidas(n_beta, beta);
        escribeMedidas(n_beta);

    //printf("me vas a escrbir");
        for(i=0; i<7; i++)fprintf(fp," %14.14lf  ",resultados[i][n_beta]);



        fprintf(fp,"\n");


        (n_beta);
    }
    fclose(fp);
}

//Realiza el histograma y la evolucion de la magnetizacion para una beta dada
void evaluaParametros(double beta, int rep_cada_iteraciones)
{
    calcula_prob(beta);
    int i;
    FILE *f = fopen("Datos/Evolucion Magnetizacion.dat", "wt");
    FILE *g = fopen("Datos/Evolucion Energia.dat", "wt");

    if(modo_de_ejecucion != 3)
        for(i=0; i<n_pasos_term; i++)
            cambiaRed();

    for(i=0; i<n_medidas; i++)
    {
        cambiaRed();
        mide(i);
        if(i%10000 == 0)
            printf("\n%d\n",i);
    }

    histograma(medida_ener, n_medidas, 50, vmin(medida_ener, n_medidas), vmax(medida_ener, n_medidas), "Datos/Histograma energia.dat");
    histograma(medida_magnet, n_medidas, 50, vmin(medida_magnet, n_medidas), vmax(medida_magnet, n_medidas), "Datos/Histograma magnetizacion.dat");
    printf("holaaaaaaaa");
   for(i=0; i<n_medidas; i++)
   {
        if(i%rep_cada_iteraciones == 0)
       {
            fprintf(f, "%lf\n", medida_magnet[i]);
            fprintf(g, "%lf\n", medida_ener[i]);
       }
   }
   fclose(f);
   fclose(g);

}
void init()
{
    leeDatos();
   // printf( "%lf", n_medidas);
     //    printf("\n")
    ini_ran(semilla_rand);
    inicializaVectoresMovimiento();
    configuracionInicial();
    inicializaMedidas();
    inicializaVectorResultados(7, fabs((b_f - b_ini)/delta_b));


}
int main (){
    init();

    if(modo_de_ejecucion <= 1)
    {
        printf("\n\nBeta\t  <e>\t\terrorE\t\t<|m|>\terrorM\t\t<e2>\t\t<m2>\n");
        iteraBetas(b_ini, b_f);
        if (modo_de_ejecucion == 1) iteraBetas(b_f, b_ini); //Sentido inverso
    }

    if(modo_de_ejecucion >= 2)
        evaluaParametros(beta_a_analizar, 100);

    system("pause");
    return 0;
    }

//Guarda la configuracion actual en una columna
//Solo usada en debug (por ahora)
void guardaConfiguracion() {

    int i;
    FILE *f = fopen(FicheroOutput, "wt");

    for(i=0; i<V; i++)
        fprintf(f,"%d\n",*(s+i));

    fclose(f);
}


