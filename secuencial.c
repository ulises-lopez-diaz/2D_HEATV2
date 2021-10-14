////////////////////////////////////////////////////////////////////////
// Ecuación del calor para un espacio 2D de forma secuencial,

// Instrucciones de ejecucion para ubuntu
// 1. Compile usando el comando:
// gcc secuencial.c -o secuencial -lm
// 2. Ejecute usando el comando:
// ./secuencial Nx Ny
// 3. Siga lo que aparece en terminal
//
// Los datos de salida se colocarán en archivos .txt 
// Diferencias finitas -> t.txt


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>


#define MIL_MILLONES 1000000000.0

void malloc_2d_mtx(float **src, int rows, int cols)
{
	int row;
	
	for (row=0; row<rows; row++){
		src[row] = (float *)malloc(cols * sizeof(float));
	}
}

void init_2d_mtx(float **src, int rows, int cols)
{
	int row, col;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			src[col][row] = 0;
		}
	}
}

void print_2d_mtx(float **mat, int rows, int cols)
{
	int row, col;
	
	if (rows<=8 && cols <= 8){
		for (col=0; col<cols*16; col++){
			printf("-");
		}
		
		printf("-\n");
		
		for (row=0; row<rows; row++){
			for (col=0; col<cols; col++){
				if (col==0){
					printf("| ");
				}
				
				printf("%f   \t|", mat[row][col]);
			}
			
			printf("\n");
			
			for (col=0; col<cols*16; col++){
				printf("-");
			}
			printf("-\n");
		}
	}
	else{
		printf("\nMatriz demasiado grande para imprimir !\n");
	}
}

void init_col_temp_ss(float **src, int col_size, int col_idx, float temp)
{
	int row;
	
	for (row=0; row<col_size; row++){
		src[row][col_idx] = temp;
	}
}

void init_row_temp_ss(float **src, int row_size, int row_idx, float temp)
{
	int col;
	
	for (col=0; col<row_size; col++){
		src[row_idx][col] = temp;
	}
}

float find_max(float **src, int rows, int cols)
{
	int row, col;
	float max = src[0][0];
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			if (src[row][col] > max) max = src[row][col];
		}
	}
	
	return max;
}

float find_min(float **src, int rows, int cols)
{
	int row, col;
	float min = src[0][0];
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			if (src[row][col] < min) min = src[row][col];
		}
	}
	
	return min;
}

void malloc_3d_mtx(float ***src, int rows, int cols, int height)
{
	int row, col, step;
	
	for (row=0; row<rows; row++){
		src[row] = (float **) malloc(sizeof(float*)*cols);
		for (col=0; col<cols; col++){
			src[row][col] = (float *) malloc(sizeof(float)*height);
		}
	}
}

void init_3d_mtx(float ***src, int rows, int cols, int height)
{
	int row, col, step;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			for (step=0; step<height; step++){
				src[row][col][step] = 0;
			}
		}
	}
}

void init_col_temp_fd(float ***src, int col_size, int col_idx, int num_steps, float temp)
{
	int row, step;
	
	for (row=0; row<col_size; row++){
		for (step=0; step<num_steps; step++){
			src[row][col_idx][step] = temp;
		}
	}
}

void init_row_temp_fd(float ***src, int row_size, int row_idx, int num_steps, float temp)
{
	int row, col, step;
	
	for (col=0; col<row_size; col++){
		for (step=0; step<num_steps; step++){
			src[row_idx][col][step] = temp;
		}
	}
}

void init_timestep_fd(float ***src, int time_idx, int rows, int cols, float temp)
{
	int row, col, step;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			src[row][col][time_idx] = temp;
		}
	}
}

void print_3d_mtx(float ***src, int rows, int cols, int time)
{
	int row, col;
	
	if (rows<=8 && cols <= 8){		
		for (row=0; row<rows; row++){
			for (col=0; col<cols; col++){
				printf("%.2f    \t", src[row][col][time]);
			}
			printf("\n");
		}
	}
	else{
		printf("\n¡Matrices demasiado grandes para imprimir! \n");
	}
}

void update_Tk(float ***src, int rows, int cols, float **dest, int k, float k1, float dt)
{
	int row, col;
	
	for (row=0; row<rows; row++){
		for (col=0; col<cols; col++){
			dest[row][col] = src[row][col][k] + k1*dt;
		}
	}
}

int main(int argc, char **argv)
{
	////////////////////////////////////////////////////////////////////
	// Seleccion de material
	////////////////////////////////////////////////////////////////////
	
	int material;
	float conductivity;	// conductividad térmica (j/m*C*sec)
	float specific_heat;	// calor especifico  (j/kg*C)
	float density;			// densidad (kg/m^3)
	
	printf("Este programa calcula la distribución de calor en una placa rectangular a lo largo del tiempo. \n\n");
	printf("Los materiales disponibles para la placa son: \n");
	printf("1 - Aluminio\n2 - Cobre\n3 - Plata\n4 - Material Personalizado\n\n");
	printf("Ingrese el número de su material a elegir: ");
	scanf("%d", &material);
	
	// valida la entrada 
	while (material<1 || material>4){
		printf("Ingrese una selección válida(1-4): ");
		scanf("%d", &material);
	}
		
	switch (material){
		case 1: // Aluminio 
			conductivity = 204.3;
			specific_heat = 910;
			density = 2700;
			break;
		case 2: // Cobre
			conductivity = 401;
			specific_heat = 390;
			density = 8940;
			break;
		case 3: // Plata
			conductivity = 629;
			specific_heat = 233;
			density = 10490;
			break;
		case 4: // Material personalizado 
			printf("Ingrese la conductividad de su material:  ");
			scanf("%f", &conductivity);
			printf("Ingrese el calor específico de su material: ");
			scanf("%f", &specific_heat);
			printf("Ingrese la densidad de su material: ");
			scanf("%f", &density);
			break;
		default:
			printf("Error fatal: selección no válida. \n");
	}
	
	////////////////////////////////////////////////////////////////////
	// Declaracion de variables
	////////////////////////////////////////////////////////////////////
	
	float Lx = 1;			// ancho de la placa  (m)
	float Ly = 1;			// longitud de la placa (m)
	int Nx = atoi(argv[1]);	// número de nodos en la dirección x 
	int Ny = atoi(argv[2]);	// número de nodos en la dirección y 
	float T_initial = 0;	// temperatura inicial en todos los nodos 
	float T_up = 50;		// temperatura en los límites (Celsius) 
	float T_down = 100;
	float T_left = 300;
	float T_right = 150; 
	float T_max = 300;
	float T_min = 0;
	float dt = 0.6;		// paso de tiempo (1 seg) 
	float tol = 0.5;		// tolerancia para simulación numérica (Celsius) 
	float tol_ss = 0.001;	// tolerancia para el estado estacionario (Celsius) 
	int i, j, k;				// contadores de iteraciones 
	float err_ss_max = 1;		// error máximo inicial de estado estable 
	float err_ss_min = 1;		// error mínimo inicial en estado estable 
	float dx = Lx/Nx;		// delta x
	float dy = Ly/Ny;		// delta y
	float alpha = conductivity/(specific_heat*density);
	
	// variables de cronometraje 
	struct timespec start, end;
	double ss_time, fd_time, sec, nsec;
	
	// prueba de condición de estabilidad 
	if (dt > 1/(2*alpha*((1/(dx*dx))+(1/(dy*dy))))){
		printf("No se cumple la condición de estabilidad, elija un dt más pequeño. \n");
		return 1;
	}
	
	////////////////////////////////////////////////////////////////////
	// Cálculo del estado estacionario 
	////////////////////////////////////////////////////////////////////	
	
	float *Tss[Nx+2], *Tss2[Nx+2], *diff[Nx+2];
	malloc_2d_mtx(Tss, Nx+2, Ny+2);	 // agregar filas adicionales para la temperatura límite 
	malloc_2d_mtx(Tss2, Nx+2, Ny+2); // agregar filas adicionales para la temperatura límite 	
	malloc_2d_mtx(diff, Nx+2, Ny+2); // agregar filas adicionales para la temperatura límite 
	init_2d_mtx(Tss, Nx+2, Ny+2);
	init_2d_mtx(Tss2, Nx+2, Ny+2);
	init_2d_mtx(diff, Nx+2, Ny+2);
	
	// inicializar los límites de la matriz de temperatura 
	init_col_temp_ss(Tss, Ny+2, 0, T_down);
	init_col_temp_ss(Tss2, Ny+2, 0, T_down);
	init_col_temp_ss(Tss, Ny+2, Ny, T_up);
	init_col_temp_ss(Tss2, Ny+2, Ny, T_up);
	init_col_temp_ss(Tss, Ny+2, Ny+1, T_up);	// utilizado para graficar 
	init_col_temp_ss(Tss2, Ny+2, Ny+1, T_up);	// utilizado para graficar 
	init_row_temp_ss(Tss, Nx+2, 0, T_left);
	init_row_temp_ss(Tss2, Nx+2, 0, T_left);
	init_row_temp_ss(Tss, Nx+2, Nx, T_right);
	init_row_temp_ss(Tss2, Nx+2, Nx, T_right);
	init_row_temp_ss(Tss, Nx+2, Nx+1, T_right);	// utilizado para graficar 
	init_row_temp_ss(Tss2, Nx+2, Nx+1, T_right);// utilizado para graficar 

	// calcular el estado estable 
	printf("Calculo del estado estable ... \n\n");
	float *temp_diff[Nx+2];
	malloc_2d_mtx(temp_diff, Nx+2, Ny+2);
	
	clock_gettime(CLOCK_MONOTONIC, &start); // empezar a contar el tiempo 
	
	while (err_ss_max>=tol_ss || err_ss_min>=tol_ss){		
		for (i=1; i<Nx; i++){
			for (j=1; j<Ny; j++){
				Tss2[i][j] = 0.25*(Tss[i+1][j]+Tss[i][j+1]+Tss[i-1][j]+Tss[i][j-1]);
			}
		}
		
		// Obtiene Tss2(i,j)-Tss(i,j)
		for (i=1; i<Nx+2; i++){
			for (j=1; j<Ny+2; j++){
				temp_diff[i][j] = Tss2[i][j]-Tss[i][j];
			}
		}	
		
		err_ss_max = fabs(find_max(temp_diff, Nx+2, Ny+2));
		err_ss_min = fabs(find_min(temp_diff, Nx+2, Ny+2));
				
		for (i=1; i<Nx+2; i++){
			for (j=1; j<Ny+2; j++){
				Tss[i][j] = Tss2[i][j];
			}
		}
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);
	sec = end.tv_sec-start.tv_sec; 
	nsec = (end.tv_nsec - start.tv_nsec)/MIL_MILLONES; // convertido a s 
	ss_time = sec + nsec;
	
	print_2d_mtx(Tss, Nx+2, Ny+2);

	printf("...finalizó el estado estable. \n\n");
	
	////////////////////////////////////////////////////////////////////
	// Cálculo de la solución en diferencias finitas 
	////////////////////////////////////////////////////////////////////
	printf("Cálculo de la solución de diferencias finitas ... \n\n");
	float k1, k2;
	int num_steps = 75000;			// número de pasos de tiempo para la matriz 
	float err_R_k_max[num_steps];		// error máximo inicial 
	float err_R_k_min[num_steps];		// error mínimo inicial 
	k = 0;							// utilizado para contar mientras que los ciclos de bucle 
	err_R_k_max[k] = 100;
	err_R_k_min[k] = 100;
	
	float **T[num_steps];
	malloc_3d_mtx(T, Nx+2, Ny+2, num_steps);
	init_3d_mtx(T, Nx+2, Ny+2, num_steps);
	
	float *Tk[Nx+2];
	malloc_2d_mtx(Tk, Nx+2, Ny+2);
	
	// inicializar los límites de la matriz de temperatura 
	init_col_temp_fd(T, Ny+2, 0, num_steps, T_down);
	init_col_temp_fd(T, Ny+2, Ny, num_steps, T_up);
	init_col_temp_fd(T, Ny+2, Ny+1, num_steps, T_up);
	init_row_temp_fd(T, Nx+2, 0, num_steps, T_left);
	init_row_temp_fd(T, Nx+2, Ny, num_steps, T_right);
	init_row_temp_fd(T, Nx+2, Ny+1, num_steps, T_right);
	init_timestep_fd(T, 0, Nx+2, Ny+2, T_initial);
	
	clock_gettime(CLOCK_MONOTONIC, &start); // Inicia el contador de tiempo de ejecución
	
	while (err_R_k_max[k]>=tol || err_R_k_min[k]>=tol){
		for (i=1; i<Nx; i++){
			for (j=1; j<Ny; j++){
				k1 = alpha*((T[i-1][j][k]-2*T[i][j][k]+T[i+1][j][k])/(dx*dx)+
							(T[i][j-1][k]-2*T[i][j][k]+T[i][j+1][k])/(dy*dy));
				update_Tk(T, Nx+2, Ny+2, Tk, k, k1, dt);
				k2 = alpha*((Tk[i-1][j]-2*Tk[i][j]+Tk[i+1][j])/(dx*dx)+
							(Tk[i][j-1]-2*Tk[i][j]+Tk[i][j+1])/(dy*dy));
				T[i][j][k+1] = T[i][j][k] + (dt/2)*(k1+k2);
			}
		}
		k++;
		
		// Obtiene T(:,:,k)-Tss
		for (i=1; i<Nx+2; i++){
			for (j=1; j<Ny+2; j++){
				temp_diff[i][j] = T[i][j][k]-Tss[i][j];
			}
		}	
						
		err_R_k_max[k] = fabs(find_max(temp_diff, Nx+2, Ny+2));
		err_R_k_min[k] = fabs(find_min(temp_diff, Nx+2, Ny+2));
				
		// Prueba de convergencia 
		if ((int)(err_R_k_max[k]*10000) == (int)(err_R_k_max[k-1]*10000) && err_R_k_max[k]!=0){
			printf("El error máximo no converge. \n");
			return 1;
		}
			
		if ((int)(err_R_k_min[k]*10000) == (int)(err_R_k_min[k-1]*10000) && err_R_k_min[k]!=0){
			printf("El error mínimo no converge. \n");
			return 1;
		}	
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end); // Finaliza el contado de tiempo de ejecución
	sec = end.tv_sec-start.tv_sec; 
	nsec = (end.tv_nsec - start.tv_nsec)/MIL_MILLONES; // convertido a s 
	fd_time = sec + nsec;
	
	float SStime = k*dt;	// utilizado para el gráfico 
	
	printf("... Cálculo de diferencia finita finalizado.\n\n");
	
	////////////////////////////////////////////////////////////////////
	// Creación de archivos de salida 
	////////////////////////////////////////////////////////////////////
	
	printf("Creando archivos de salida ... \n\n");
	
	// creacion de archivos 
	//FILE *tss;
	FILE *t;
	FILE *int_vars;
	FILE *f_vars;
		
	
	// imprime matriz Tss 
	// tss = fopen("tss.txt", "w+");

	// if (tss == NULL){
	// 	printf("Error: archivo no creado. \n");
	// 	return 1;
	// }
	
	// for (i=0; i<Nx+2; i++){
	// 	for (j=0; j<Ny+2; j++){
	// 		fprintf(tss, "%f\t", Tss[i][j]);
	// 	}
	// 	fprintf(tss, "\n");
	// }
	// fclose(tss);
	
	// imprime matriz T 
	int step;
	
	t = fopen("t.txt", "w+");
	
	for (step=0; step<k; step++){
		for (i=0; i<Nx+2; i++){
			for (j=0; j<Ny+2; j++){
				fprintf(t, "%f\t", T[i][j][step]);
			}
			fprintf(t, "\n");
		}
	}
	
	// imprimir variables enteras asociadas con Tss y T 
	int_vars = fopen("int_vars.txt", "w+");
	
	if (int_vars == NULL){
		printf("Error: File not created.\n");
		return 1;
	}
	
	fprintf(int_vars, "%d\n", Nx);
	fprintf(int_vars, "%d\n", Ny);
	fprintf(int_vars, "%d\n", k);
	
	fclose(int_vars);
	
	// imprimir variables flotantes asociadas con Tss y T 
	f_vars = fopen("f_vars.txt", "w+");
	
	if (f_vars == NULL){
		printf("Error: archivo no creado. \n");
		return 1;
	}
	
	fprintf(f_vars, "%f\n", dx);
	fprintf(f_vars, "%f\n", dy);
	fprintf(f_vars, "%f\n", SStime);
	fprintf(f_vars, "%f\n", Lx);
	fprintf(f_vars, "%f\n", Ly);
	fprintf(f_vars, "%f\n", T_max);	// temperatura máxima en los límites 
	fprintf(f_vars, "%f\n", T_min);	// temperatura mínima en los límites 
	
	printf("... archivos de salida creados con éxito. Programa terminado. \n\n");
	
	printf("Tiempos de ejecución (segundos): \n");
	printf("Diferencias finitas \n");
	printf("%.9f\n", fd_time);
	
	return 0;
}
