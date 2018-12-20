#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define N_BASES 4
#define N_SITIOS 9

typedef struct {
    float min_T, max_T;
    float min_init_h, max_init_h;
    float min_init_J, max_init_J;
} LimitesParametros;

typedef struct {
    float T;
    float * h;
    float * J;
} Parametros;

typedef struct {
    float * p_i;
    float * p_ij;
} Muestras;

long indice_unidimensional_h(int posicion_x, int base_x){
    return N_BASES * posicion_x + base_x;
}

long indice_unidimensional_J(int posicion_x, int posicion_y, int base_x, int base_y){//Funcion para transformar una coordenada en el espacio sitios-bases en un indice unidimensional que voy a usar para todos mis vectores.
        long indice_x = N_BASES * posicion_x + base_x;
        long indice_y = N_BASES * posicion_y + base_y;
        return (N_BASES * N_SITIOS) * indice_y + indice_x;
}

void init_parametros(LimitesParametros * lims, Parametros * parametros){
    //ASIGNACION.
    //Asigno memoria dentro de los campos de parametros.
    int num_entradas = N_BASES * N_SITIOS;
    (parametros -> h) = (float *) malloc(num_entradas * sizeof(float));
    (parametros -> J) = (float *) malloc(num_entradas * num_entradas * sizeof(float));
    
    //INICIALIZACION.
	//Ahora a llenar con defaults.
    //La temperatura es parametro obviamente.
    parametros -> T = lims->min_T + (lims->max_T - lims->min_T) * ((float) rand()) / ((float) RAND_MAX);
    
	//Energia local h.
    long i_b, i_s;
    for(i_b=0; i_b<N_BASES; i_b++){
        for(i_s=0; i_s<N_SITIOS; i_s++){
            *((parametros -> h) + indice_unidimensional_h(i_s, i_b)) = lims->min_init_h + (lims->max_init_h - lims->min_init_h) * ((float) rand()) / ((float) RAND_MAX);
        }
    }

	//Doble interaccion J.
	long j_b, j_s;
    for(i_b=0; i_b<N_BASES; i_b++){
        for(i_s=0; i_s<N_SITIOS; i_s++){
            for(j_b=0; j_b<N_BASES; j_b++){
                for(j_s=0; j_s<=i_s; j_s++){//Hay simetria ante un cambio de posicion. No voy mas alla de i en el contador de j.
                    *((parametros -> J) + indice_unidimensional_J(i_s, j_s, i_b, j_b)) = lims->min_init_J + (lims->max_init_J - lims->min_init_J) * (((float) rand()) / ((float) RAND_MAX));
                    *((parametros -> J) + indice_unidimensional_J(j_s, i_s, j_b, i_b)) = *((parametros -> J) + indice_unidimensional_J(i_s, j_s, i_b, j_b));//Las componentes simetricas.
                }
            }
        }
    }
}

void free_parametros(Parametros * params, int n_parametros){
    int i;
    for(i=0; i<n_parametros; i++){
        free((params + i) -> h);
        free((params + i) -> J);
    }
    free(params);
}

void print_params(Parametros * params){
    int i_b, i_s, j_b, j_s;
    for(i_b=0; i_b<N_BASES; i_b++){
        for(i_s=0; i_s<N_SITIOS; i_s++){
            for(j_b=0; j_b<N_BASES; j_b++){//Como interaccion de dos cuerpos hay simetria en un cambio de bases i y j. No voy mas alla de i en el contador de j.
                for(j_s=0; j_s<N_SITIOS; j_s++){//Mismo tipo de simetria ante un cambio de posicion, y ante un cambio de ambos.
                    printf("Par (%d,%d) (%d,%d) %f \n", i_s, i_b, j_s, j_b, *(params -> J + indice_unidimensional_J(i_s, j_s, i_b, j_b)));
                }
            }
            printf("\n");
        }
    }
}

void print_chain(int ** cadena){
    int idx_nucleotido, i;
    for(i=0; i<N_SITIOS; i++){
        idx_nucleotido = *(*cadena + i);
        switch(idx_nucleotido){//Como esta funcion es solo para debuggear y visualizar, no me importa la eficiencia y uso un switch_case;
            case 1: printf("C"); break;
            case 2: printf("T"); break;
            case 3: printf("G"); break;
            default: printf("A");
        }
    }
    printf("\n");
}

//Para visualizar las probabilidades marginales de orden 1. Las de orden 2 seria confuso visualizar.
void print_muestras(Muestras * muestras){
    printf("\n");
    int i_b, i_s;
    for(i_b=0; i_b<N_BASES; i_b++){//Divido el muestreo por la cantidad de muestras que tome para obtener el promedio.
        
        switch(i_b){
            case 1: printf("C :: "); break;
            case 2: printf("T :: "); break;
            case 3: printf("G :: "); break;
            default: printf("A :: ");
        }
        
        for(i_s=0; i_s<N_SITIOS; i_s++){
            printf("%f ", *((muestras -> p_i) + indice_unidimensional_h(i_s, i_b)));
        }
        printf("\n");
    }
    printf("\n");
}

float step_evolucion(int ** cadena, Parametros * parametros){
    float delta_energia;

    int base_modificar = rand() % N_BASES;
    int sitio_modificar = rand() % N_SITIOS;
    int base_actual = *(*cadena + sitio_modificar);
    if(base_actual != base_modificar){//Si la eleccion me exigiria un cambio.
        float delta_local = ((parametros -> h) + indice_unidimensional_h(sitio_modificar, base_modificar)) - ((parametros -> h) + indice_unidimensional_h(sitio_modificar, base_modificar));//Resto la energia local actual de la posible energia local futura.
        float delta_pares = 0;//Inicializo la energia de a pares a cero.
        int j;
        for(j=0; j<N_SITIOS; j++){//Itero sobre todos los sitios (pares).
            if(j != sitio_modificar){//Obviamente no considero al par del sitio con si mismo. Quiza optimizar esto en el futuro dividiendo el ciclo for() en dos?
                delta_pares += *((parametros -> J) + indice_unidimensional_J(sitio_modificar, base_modificar, j, *(*cadena + j)));
                delta_pares -= *((parametros -> J) + indice_unidimensional_J(sitio_modificar, base_actual, j, *(*cadena + j)));
            }
        }

        float delta_energia = delta_local + delta_pares;
        if(delta_energia > 0){
            float probabilidad_aceptacion = ((float) rand()) / ((float) RAND_MAX);
            if(exp(-delta_energia/(parametros -> T)) > probabilidad_aceptacion){//Muestreo Metropolis por rechazo.
                *(*cadena + sitio_modificar) = base_modificar;//Actualizo la cadena.
            }
        }
        else{//Si la energia no es mayor a cero, el cambio esta garantizado.
            *(*cadena + sitio_modificar) = base_modificar;//Actualizo la cadena.
        }
    }
    else{//Si no me propone un cambio, entonces simplemente actualizo el cambio medio en energia con un cero para esta iteracion.
        delta_energia = 0;
    }
    
    return delta_energia;
}

void evolucion_sistema(Parametros * parametros, Muestras * muestras){
    int n_inicializaciones = 100;//Reinicializaciones para muestrear no solo en tiempo sino en distribucion tambien (para controlar por falta de ergodicidad).
    int n_muestras = 1000;//Muestreo habiendo convergido a un estacionario.
    float criterio_convergencia = 1.01;//Paro de evolucionar cuando mis cambios de energia esten a un 1% de las fluctuaciones termicas.
    
    //Preparo la estructura de muestreo.
    //Asigno memoria dentro de los campos de muestras.
    int num_entradas = N_BASES * N_SITIOS;
    (muestras -> p_i) = (float *) malloc(num_entradas * sizeof(float));
    (muestras -> p_ij) = (float *) malloc(num_entradas * num_entradas * sizeof(float));
    int i_b, i_s, j_b, j_s;
    for(i_b=0; i_b<N_BASES; i_b++){//Inicializo todos los muestreos a cero porque van a ser cumulativos.
        for(i_s=0; i_s<N_SITIOS; i_s++){
            *((muestras -> p_i) + indice_unidimensional_h(i_s, i_b)) = 0;

            for(j_b=0; j_b<=N_BASES; j_b++){
                for(j_s=0; j_s<=i_s; j_s++){//Hay simetria ante un cambio de posicion. No voy mas alla de i en el contador de j.
                    *((muestras -> p_ij) + indice_unidimensional_J(i_s, j_s, i_b, j_b)) = 0;
                    *((muestras -> p_ij) + indice_unidimensional_J(j_s, i_s, j_b, i_b)) = *((muestras -> p_ij) + indice_unidimensional_J(i_s, j_s, i_b, j_b));//Las componentes simetricas.
                }
            }
        }
    }
    
    int r;
    for(r=0; r<n_inicializaciones; r++){//Cantidad de veces que voy a reinicializar aleatoriamente por si el proceso no es ergodico.
        int * cadena = (int *) malloc(N_SITIOS * sizeof(int));
        int i;
        for(i=0; i<N_SITIOS; i++){
            *(cadena + i) = rand() % N_BASES;
        }
        
        int t;
        float fluctuaciones_termicas = (parametros -> T) * (parametros -> T);//La varianza de la distribucion exponencial en dominio positivo (o sea decir que ya no hay descensos de energia). PRO TIP: es mucho mas eficiente multiplicar dos veces que usar la funcion pow(); evitar esa funcion como la peste.
        float accum_delta_energia = 2 * criterio_convergencia * fluctuaciones_termicas;//Los inicializo en valores que me garanticen entrar en el loop while().
        float moving_avg_coef = 0.01;//El coeficiente para el promedio en movimiento.
        
        while(accum_delta_energia > criterio_convergencia * fluctuaciones_termicas){//Iterar hasta estar lo suficientemente cerca a fluctuaciones termicas en el estacionario.
            float delta_energia = step_evolucion(&cadena, parametros);
            accum_delta_energia = (1.0 - moving_avg_coef) * accum_delta_energia + moving_avg_coef * (delta_energia * delta_energia);//Actualizo el delta cumulativo nuevo.
        }

        //Habiendo convergido, es la etapa de muestreo.
        for(t=0; t<n_muestras; t++){
            //Primero solo evoluciono el sistema un paso, como antes.
            step_evolucion(&cadena, parametros);
            
            //Esta segunda parte es distinta y es la de muestreo.
            for(i_b=0; i_b<N_BASES; i_b++){
                for(i_s=0; i_s<N_SITIOS; i_s++){
                    *((muestras -> p_i) + indice_unidimensional_h(i_s, i_b)) += (*(cadena + i_s) == i_b);//Es C entonces los ints son booleanos y puedo hacer hermosuras asi :)

                    for(j_b=0; j_b<=N_BASES; j_b++){
                        for(j_s=0; j_s<=i_s; j_s++){//Hay simetria ante un cambio de posicion. No voy mas alla de i en el contador de j.
                            *((muestras -> p_ij) + indice_unidimensional_J(i_s, j_s, i_b, j_b)) += (*(cadena + i_s) == i_b) && (*(cadena + j_s) == j_b);
                            *((muestras -> p_ij) + indice_unidimensional_J(j_s, i_s, j_b, i_b)) = *((muestras -> p_ij) + indice_unidimensional_J(i_s, j_s, i_b, j_b));//Las componentes simetricas.
                        }
                    }
                }
            }
        }
    }
    
    //Normalizo las muestras ahora.
    for(i_b=0; i_b<N_BASES; i_b++){//Divido el muestreo por la cantidad de muestras que tome para obtener el promedio.
        for(i_s=0; i_s<N_SITIOS; i_s++){
            *((muestras -> p_i) + indice_unidimensional_h(i_s, i_b)) /= (n_muestras * n_inicializaciones);

            for(j_b=0; j_b<=N_BASES; j_b++){
                for(j_s=0; j_s<=i_s; j_s++){//Hay simetria ante un cambio de posicion. No voy mas alla de i en el contador de j.
                    *((muestras -> p_ij) + indice_unidimensional_J(i_s, j_s, i_b, j_b)) /= (n_muestras * n_inicializaciones);
                    *((muestras -> p_ij) + indice_unidimensional_J(j_s, i_s, j_b, i_b)) /= (n_muestras * n_inicializaciones);
                }
            }
        }
    }
}

int main(){
	srand(time(NULL));

    LimitesParametros * limites_parametros = (LimitesParametros *) malloc(sizeof(LimitesParametros));
    limites_parametros -> min_T = 0.1;
    limites_parametros -> max_T = 0.2;
    limites_parametros -> min_init_h = 0.0;
    limites_parametros -> max_init_h = 1.0;
    limites_parametros -> min_init_J = 0.0;
    limites_parametros -> max_init_J = 1.0;

    //Inicializo el barrido que hare sobre el espacio de parametros.
    int n_parametros = 1;//Cantidad de particulas en mi espacio de parametros.
    Parametros * parametros = (Parametros *) malloc(n_parametros * sizeof(Parametros));
    int i;
    for(i=0; i<n_parametros; i++){
        init_parametros(limites_parametros, parametros);
    }

    Muestras * muestras = (Muestras *) malloc(sizeof(Muestras));
    evolucion_sistema(parametros, muestras);
    print_muestras(muestras);
    
    free_parametros(parametros, n_parametros);

    return 0;
}
