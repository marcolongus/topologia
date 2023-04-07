#include <bits/stdc++.h>
#include "classparticle.h"

#define forn(i,a,b) for(int i=a; i<b; i++)

using namespace std;

int main(void){

    ofstream states(archivo); //En este archivo guardamos #infect, refrect, sanos y tiempo.	
	ofstream pico(archivo1);
	ofstream epid(archivo2);
	ofstream anim("data/animacion.txt");
	ofstream topology("data/topologia.txt");
	ofstream mipsframe("data/mips.txt");

	//#pragma omp parallel for
    for (int n_simul=0; n_simul<1 ; n_simul++){

    	cout << "Semilla: "  << seed << endl;
    	cout << endl;
    	//topology << "SIMULACION: " << n_simul << endl;
    	int start_s=clock();

    	cout <<"simulacion: "<< n_simul << endl << endl;

	    /******************************************************************************************************/
	    /*********************************  Declaración de variables   ****************************************/
	    /******************************************************************************************************/

	    vector<particle> system,system_new; //Vector que contiene las particulas y su información relevante.
	    vector<bool>     inter,             //Flag de iteracción para la función evolución. 
	                     inter_old;         //Nos avisa cuándo terminó la interacción.    

	    int num_infected   = 0, 
	        num_refractary = 0, 
	        num_sane       = 0;               //contadores de la evolución de la epidemia.
	    
	    vector<vector<set<int>>> box; //Cuadriculas de set.
	    
	    int num_boxes=floor(L);       //Cantidad de cuadrículas en un lado.
	    
	    /******************************************************************/
	    /*Definimos un vector booleano de interacción*/

	    inter.resize(N);
	    inter_old.resize(N);

	    for(int i; i<N; i++) {
	        inter[i]=false;
	        inter_old[i]=false;
	    }


	    /*========================================================================================================================*/
	    /*Crafting the grid of the square LxL. Each box has side=1 */	    
	    box.resize(num_boxes);
	    for (int i=0; i<box.size(); i++) box[i].resize(num_boxes);

	    /*==========================================*/
	    /* Set para medir R_0 */	
	    /*==========================================*/
		vector<set<int>> R0_set; 
		R0_set.resize(N);
		set<int> first_set, last_set; //ponemos las infectadas iniciales en este set

	    /*========================================================================================================================*/
	    /**************** Condición inicial: todos las partículas se inician de forma que no interactúen **************************/
	    /*========================================================================================================================*/
	    
	    for (int i=0; i<N; i++){
	 		
	 		particle A;           //Definimos una partícula A. 
	        bool accepted=false;  //Variable de aceptación de la partícula. True si no interactúa con otras y false caso contrario.
	        
	        while (!accepted) {
	            
	            accepted=true;
	            A=create_particle();
	   
	            int i_index=floor(A.x), j_index=floor(A.y);
				box[i_index][j_index].insert(i);

	            /* Si interactúa con otra partícula cambiamos la condicion a no aceptada*/
	        
	            for (int l=-2; l<3; l++) {
	                for(int m=-2; m<3; m++){
	        
	                    if (! box[b_condition(i_index+l)][b_condition(j_index+m)].empty() ){
	                    
	                        for(auto elem: box[b_condition(i_index+l)][b_condition(j_index+m)]){
	                        
	                            if (  (elem !=i ) &&  interact(A,system[elem]) ) {
	                                
	                                //cout << i<< " interact with " << elem <<endl;
	                                accepted=false;
	                                box[i_index][j_index].erase(i);

	                            }
	                        }
	                    }    
	                }
	            }        	   
	        } //Si accepted false, el while vuelve a generar una partícula

	        if(i<n_inmunes) A.set_inmune();// Agrego n_inmunes
	        if(i<10) A.set_infected();

	        if (A.is_infected())   num_infected++;
	        if (A.is_sane())       num_sane++;
	        if (A.is_refractary()) num_refractary++;
	        

	        system.push_back(A); //Si fue aceptada la guardamos en el sistema de partículas. 
        
	    } //cierra el for de generacion de partículas

	   // Agregarmos las primeras infectadas
	    forn(i, 0, N){
	    	if (system[i].is_infected()) {
	    		cout << i  << " " << system[i].velocity << endl; 
	    		first_set.insert(i);
	    	}
	    }
	    cout << endl;

	    system_new.resize(system.size());

	    cout << "Infectadas   iniciales = "  << num_infected   << endl;
	    cout << "refractarias iniciales = "  << num_refractary << endl;
	    cout << "Sanas iniciales        = "  << num_sane       << endl;
	    cout << endl;
	  

	    /*========================================================================================================================*/
	    /*Declaración de variables para la evolución del sistema*/
	    /*========================================================================================================================*/

	    double time=time_i;

	    int contador1=0, //Cantidad de interacciones dobles
	        contador2=0, //Cantidad de evoluciones
	        contador3=0, //Cantidad de whiles
	        contador4=0, //Cantidad de interacciones triples
	        contador5=0, //Cantidad de interacciones cuadruples
	        contador6=0; //Cantidad de interaciones quíntuples

	    int int_tot=0; //Cantidad total de interacciones (suponiendo que no hay sextuples)

 
	    /****************************************************************************************************************/
	    /*********************************** Evolución del sistema ******************************************************/
	    /****************************************************************************************************************/

	    double i_max = num_infected,
	    	   t_imax = 0;

	    bool bandera      = false; //para  guardar animación
		bool mips_search  = false; //sarch mips	 
	    int frame_counter = 0;  


		while(num_infected>0){

	            contador3++; //Cantidad de while
				
				if ( (contador3 % 10000) == 0 ){cout << "Time " << time << endl;}	          
	    	    time=time+delta_time;

	            num_sane=0;
	            num_infected=0;
	            num_refractary=0;

	            #pragma omp parallel for
	    	    for (int i = 0; i < N; i++){

	                vector<int> index;
	                index.push_back(i);
	                inter[i]=false;
	                           
	                /*Chequeamos si hay interacciones*/

	                for (int l=-2; l<3; l++) {
	            	   for(int m=-2; m<3; m++){

	                        if(! box[b_condition( floor(system[i].x)+l )][b_condition( floor(system[i].y)+m )].empty() ){

	            		        for(auto elem: box[b_condition( floor(system[i].x)+l )][b_condition( floor(system[i].y)+m )]){
	                 	
	                        	    if (  (elem !=i ) &&  interact(system[i],system[elem]) ) {
	                                
	                    			    inter[i]=true;
	                                    inter_old[i]=true;

	                    			    index.push_back(elem);

	                    			    if (system[i].is_infected() and system[elem].is_sane()) R0_set[i].insert(elem);                                
								               		                             
	                                }//if interaction.
	                            }//for auto elem.
	                        }//if box isn't empty.
	                    }//for m.
	                }//for l.

	                /************************************Fin de chequeo de interacciones******************************************/
	                
	                system_new[i]=evolution( system, index, inter[i] ); //Evolución del sistema.
	                
	                contador2++; // Cantidad de evoluciones.        	
	                
	            } //cierra el for de n particulas

	            
	            forn(i,0,N){

	                if ( contador3 %  1 == 0 ) bandera = false; //si pongo false apago la animación. 
                	if (bandera){

                    	anim << system_new[i].x << " ";
                    	anim << system_new[i].y << " ";                    
                    	anim << time            << " ";
                	}

	            	if(system_new[i].is_infected()) {
                    	num_infected++;
                    	if (bandera) {anim << 1 << endl; bandera=false;}
                	}

                	if(system_new[i].is_sane()){
                    	num_sane++; 
                    	if(bandera) {anim << 0 << endl; bandera=false;}
                	}

                	if(system_new[i].is_refractary()) {
                    	num_refractary++;
						if(bandera) {anim << 2 << endl; bandera=false;}
                	}
	            } //for de particulas-escritura de datos

	      		if (i_max < num_infected) {

	      			frame_counter++;
	      			i_max=num_infected; 
	      			t_imax=time;
	      			if(mips_search){
		      			forn(i,0,N){
	                    	mipsframe << system_new[i].x << " ";
	                    	mipsframe << system_new[i].y << " ";                    
	                    	mipsframe << time            << " ";
	                    	if(system_new[i].is_infected())   mipsframe << 1 << endl;
	                    	if(system_new[i].is_sane())       mipsframe << 0 << endl;
	                    	if(system_new[i].is_refractary()) mipsframe << 2 << endl;
		      			}
	      			}
	      		} //if i_max-frame_mips


	            /*Estabilzamos el set*/

	            for(int j=0; j<N; j++){

	                if ( box[floor(system_new[j].x)][floor(system_new[j].y)].find(j)==box[floor(system_new[j].x)][floor(system_new[j].y)].end() ){

	                    box[floor(system[j].x)][floor(system[j].y)].erase(j);
	                    box[floor(system_new[j].x)][floor(system_new[j].y)].insert(j);

	                }
	            }//cirra el for j.

	            /**********************************************************************************************************************************/

	            if (contador3 % 1 == 0){
	            	epid << num_sane       << " ";
	            	epid << num_infected   << " ";
	            	epid << num_refractary << " ";	            	
	            	epid << time           << " ";
	            	epid << n_simul        << endl;
	            }

	            system=system_new; 

	            forn(i,0,N){
	            	for(auto elem: R0_set[i]){
	            		if (system[elem].is_sane()) R0_set[i].erase(elem);
	            	}// for auto

	            }//forn	          	
		}//cierra el while

	    /*Escritura de nodos y vertices*/

	    //for para sacar repetidas por colisión múltiple.
	    forn(i,0,N){
	    	bool flag=true;
	    	int multiple_inter=0;	    	
	    	forn(j,0,N){
	    		if ( !(R0_set[j].find(i)==R0_set[j].end()) ) {
	    			multiple_inter++; 
	    			flag=false;
	    		}
	    		if (multiple_inter>1 and !flag) 
	    			R0_set[j].erase(i);
	    	}
	    }

	    int epidemy_size=0;
	    double R0_mean=0, n_R0;
	    double R0_first_gen=0;	    

	    //For para agregar las que no contagian a nadie. 
	    forn(i,0,N){
	    	bool flag = R0_set[i].empty();
	    	if (flag){
	    		//Las de la última generación que no producen descendencia. 
	    		forn(j,0,N){
	    			if ( !(R0_set[j].find(i) == R0_set[j].end()) ) {
	    				last_set.insert(i);
	    			}		    		
	    		}
	    		// Las de la primera generación que no producen descendencia. 
	    		if ( !(first_set.find(i) == first_set.end()) )  {
	    			topology << i << endl;	    			
	    		}    		
	    	}
	    }

	    for(auto elem: last_set) {topology << elem << endl;}

		forn(i,0,N){
			if(!R0_set[i].empty()){

				R0_mean+=R0_set[i].size();
				n_R0++;
				
				topology << i << " ";			
	    		for (auto elem: R0_set[i]) {
	    			topology << elem  <<" ";
	    			epidemy_size++;	    		
	    		}
	    		topology << endl;
	    	}	    	
	    }

	    cout << "size of epidemic = " << epidemy_size + first_set.size() << endl;
	    cout << endl;


	    /****************************************************************************/
	    /****************************************************************************/

	    states << num_refractary << " ";
	    states << time << endl;

	    pico << i_max  << " ";
	    pico << t_imax << endl;


	    /**********************************************************************************************************************/
	    /**********************************************************************************************************************/
	    /**********************************************************************************************************************/

	    /*Escritura de Resultados*/
	   	    
	    int stop_s=clock(); //Para el reloj para calcular la duración total del programa. 
	    int cont_inmunes=0, cr=0, ci=0, cs=0 ;

	    for (int i=0; i<N; i++){

	        if(system[i].is_inmune())     cont_inmunes++;
	        if(system[i].is_refractary()) cr++;
	        if(system[i].is_infected())   ci++;
	        if(system[i].is_sane())       cs++;
	    }
	    
	    cout << "Cantidad de agentes= " << cont_inmunes+cr+ci+cs << endl;

	    cout << "Refract = " <<  cr           << endl;
	    cout << "Infect  = " <<  ci           << endl;
	    cout << "Sanos   = " <<  cs           << endl;
	    cout << endl;

	    cout <<"Time in min:" << (((stop_s-start_s)/double(CLOCKS_PER_SEC)*1000)/1000)/60 <<endl;
	    cout << "maximos: " << frame_counter << endl;

	    /*Guardo los resultados finales de cada simulación en el archivo simulaciones.txt*/

	    /*************************************************************************************/

    }//Simulaciones
  
    states.close();
    pico.close();
    epid.close();
    anim.close();
    topology.close();
    mipsframe.close();
        
    /**********************************************************************************************************************/
    /**********************************************************************************************************************/


    return 0;   
}
