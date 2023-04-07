using namespace std;

/*Declaración de constantes*/

const int    N = 1000;                //Number of particles.
const double L = 150;                 //Lado del cuadrado

const string archivo  = "data/evolution.txt";//Cambiar el nombre del archivo con el L que corresponda. 
const string archivo1 = "data/imax.txt";
const string archivo2 = "data/epid.txt";

const int n_inmunes = 0;            // Cantidad de inmunes

const double delta_time = 0.1;        // Evoltion time step of the system.

const double time_i= 0, time_f = 1500;  //Intervalo temporal de evolución del sistema.
const double time_prom = 10;          //Tiempo a partir del cual consideramos que se llegó al equilibrio. 

const double active_velocity = 0.05;   //Velocidad de las partículas cuando no están interactuando.
 
const double infinity = 1000000000;
const double radius   = 1;              //radio de la partícula
const double diameter = 2*radius;     //diámetro de la partícula
const double error    = 0.000001;

/*Parámetros de la simulación*/

const double beta=1;
const double gamma_friction=3.92*active_velocity;

/*Tiempos característicos y probabilidades por unidad de tiempo*/

const double tau_t = 1  , p_transmision = (1/tau_t)*delta_time; //sane---->infected 
const double tau_i = 200, p_infection   = (1/tau_i)*delta_time; //infected--->refractary 
const double tau_r = 500, p_recfractary = (1/tau_r)*delta_time; //refractary--->sane

const double alpha=100, p_rotation    = (1/alpha)*delta_time; 

const double p_init  = 0.0;//población inicial infectada. 
const double p_rinit = 0.0;  //población inicial de refractarios 
const double p_dinit = 0.1;  //población inicial de refractarios dirigidos. 

const double Pi     = 3.14159265358979323846;
const double dos_Pi = 2*Pi;

const double k_powerl= 2.09500; //1.7713537054317126; // v=0.1 //Constante de la power-law distribution
const double v_min   = 0.01, 
             v_max   = 4.;

/***************************************************************************************/
/***************************************************************************************/

//Generador de números aleatorios entre 0 y 1//

mt19937::result_type seed = time(0);
mt19937 gen(seed);                             //Standard mersenne_twister_engine seeded time(0)
uniform_real_distribution<double> dis(0., 1.); // dis(gen), número aleatorio real entre 0 y 1. 

/************************************************************************************************/
   
/* Definimos la clase partículas y sus métodos */

class particle{
    
    int state;
    
    public:
    
    double x,y;
    double velocity;
    double angle;

        
    particle(); //constructor of a default particle. 
    particle(double x1, double y1, double vel, double ang);
    
    bool is_inmune()     { return  (state==3)  ;} 
    bool is_refractary() { return  (state==2)  ;} //Ask for the state 2,1,0.
    bool is_infected()   { return  (state==1)  ;}
    bool is_sane()       { return  (state<1 )  ;}
    
    void set_inmune()     {state=3 ;}
    void set_refractary() {state=2 ;}
    void set_infected()   {state=1 ;}
    void set_sane()       {state=0 ;}
};

/***************************************************************************************/
/***************************************************************************************/

/* Constructor for a default particle */

particle::particle(){
    x=0;
    y=0;
    velocity=0;
    angle=0;
  
}

/***************************************************************************************/

/* Constructor of a partcile in a given phase-state (x,p) of the system */

particle::particle(double x1, double y1, double vel, double ang){
    x=x1;
    y=y1;
    velocity=vel;
    angle=ang;

}

/***************************************************************************************/

/* Create a particle in tha random phase-state  */

particle create_particle(void){
    
    double x,y,velocity,angle;
     
    x=dis(gen)*L;
    y=dis(gen)*L;
    angle=dis(gen)*dos_Pi;

    //velocity=-active_velocity*log(1.-dis(gen)); //distribución exponencial
    velocity= pow( dis(gen)*( pow(v_max,1-k_powerl) - pow(v_min,1-k_powerl))+pow(v_min,1-k_powerl), 1./(1.-k_powerl)); //power_law
    //velocity=active_velocity;
    
    particle A(x,y,velocity,angle);
    
    if( dis(gen) < p_init){ A.set_infected();}
    else A.set_sane();
    if ( dis(gen) < p_rinit and !A.is_infected() ){A.set_refractary();}
   
    return A;
}

/*******************************************/
/*******************************************/

/* Function for the boundary condition  and Integer boundary condition function */

double b_condition(double a){
    return fmod((fmod(a,L)+L),L);   
}

int my_mod(int a, int b){
	return ((a%b)+b)%b;
}

/*********************************************************************************/
/*********************************************************************************/


/*Distance between particles function.*/

double distance(particle A, particle B){
        double x1,x2,y1,y2,res;
        res = infinity;
        x2 = B.x; y2 = B.y;
        for(int i=-1;i<2;i++) for(int j=-1;j<2;j++){
            x1 = A.x + i*L;
            y1 = A.y + j*L;
            res = min(res, pow((x1-x2),2) + pow((y1-y2),2));
        }
        return sqrt(res);
}
/***************************************************************************************/

double distance_x(particle A, particle B){

        double x1,x2,res;
        int j=0;

        vector<double> dx;
        dx.resize(3,0);

        res = infinity;
        x2 = B.x;

        for(int i=-1;i<2;i++){

            x1 = A.x + i*L; 
            dx[i+1]=x1-x2;

            if (abs(dx[i+1]) < res ) {

                res=  abs(dx[i+1]);
                j=i;   
            } //if           
        }//for

        return dx[j+1];
}

double distance_y(particle A, particle B){
    
        double y1,y2,res;
        int j=0;

        vector<double> dy;
        dy.resize(3,0);

        res = infinity;
        y2 = B.y;

        for(int i=-1;i<2;i++){

            y1 = A.y + i*L; 
            dy[i+1]=y1-y2;

            if (abs(dy[i+1]) < res ) {

                res=  abs(dy[i+1]);
                j=i;   
            } //if       
        }//for

        return dy[j+1];
}

double distance1(double dx, double dy){
    
    double d=sqrt( pow(dx,2) + pow(dy,2) );

    return d;
}


/* Interact */

bool interact(particle A, particle B){ return (distance(A,B)<diameter);} //repensar esta función

/***************************************************************************************/
/***************************************************************************************/


/* Evolution time step function of the particle */
/* Campo de interacción */


vector<double> campo( vector<particle> system, vector<int> &index){
 
    vector<double> field; //Campo de salida
    vector<double> potencial;
    
    field.resize(2);
    potencial.resize(2,0); //inicia vector tamaño 2 en 0.

    for(int i=1; i < index.size(); i++){

        double dx_0i=distance_x(system[index[0]],system[index[i]]);
        double dy_0i=distance_y(system[index[0]],system[index[i]]);

        double d_0i=distance1(dx_0i , dy_0i);

        potencial[0]=pow( d_0i,-3 )*dx_0i + potencial[0];
        potencial[1]=pow( d_0i,-3 )*dy_0i + potencial[1];

    }
    
    for(int i=0; i<potencial.size(); i++) {potencial[i]=gamma_friction*potencial[i];}

    field[0] = system[index[0]].velocity*cos(system[index[0]].angle) + potencial[0]; 
    field[1] = system[index[0]].velocity*sin(system[index[0]].angle) + potencial[1]; 
    
    return field;
}

/***************************************************************************************/



particle evolution(vector<particle> &system, vector<int> &index, bool inter){
   
    particle A=system[index[0]];
    bool flag=true; //Flag de infección usada en la dinámica de la epidemia.

    /* Dinámica espacial del agente*/
    
    if (inter){

        vector<double> k; 
        k.resize(2);
        
        k=campo(system,index);  //campo del sistema
    
        A.x = b_condition(A.x + delta_time*k[0]);
        A.y = b_condition(A.y + delta_time*k[1]);        

        }

    else{

        A.x = b_condition(A.x + A.velocity*cos(A.angle)*delta_time);
        A.y = b_condition(A.y + A.velocity*sin(A.angle)*delta_time);
        
    }

    /****************************************************************/
    
    /*Dinámica de la epedemia*/

    for (int i=1; i<index.size(); i++){
        if (A.is_sane() && system[index[i]].is_infected()){
            if (dis(gen) < p_transmision) {

                A.set_infected();
                flag=false;
            }
        }
    }


    //if ( A.is_refractary() ){ if(dis(gen)<p_recfractary){A.set_sane();} }
    if (A.is_infected() && flag){ if(dis(gen)<p_infection){A.set_refractary();} }

    /******************************************************************************/
    
    /*El agente cambiar de dirección en A.angle +/- pi/2*/ 

    if (dis(gen)<p_rotation) {

        double ruido= 2 * dis(gen) * Pi;

        A.angle+=ruido; 
    }

    /****************************************************/    

    return A;       
}



/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/





