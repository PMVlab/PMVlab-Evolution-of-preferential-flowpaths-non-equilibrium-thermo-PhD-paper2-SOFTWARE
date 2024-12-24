#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include "nrutil.h"
#include "nr.h"
#include <math.h>
#include <memory.h>
#include <string.h>

/*====================================================================             */
/*   Routine for the simulation of a reactive transport in 2D                      */
/*     porous media on a Darcy scale using Lagrangian Particle tracking approach   */
/*                                                                                 */
/*					         	                                                   */    
/*                   LAST MODIFIED BY EVGENY SHAVELZON 24-1-2024                   */
/*====================================================================             */

/*GLOBAL SCOPE VARIABLES*/
/*definitions*/
#define GridScale 0.2 /*cell dimension [cm]*/
#define Porosity_lower_bound 0.01 /*porosity limits for reaction*/
#define Porosity_upper_bound 0.99
#define Frequency_of_output 10 /*output / flow calc frequency*/
#define CaCO3_MV 36.93

/*global vars*/
int numcellX, numcellY; /*number of cells in X, Y directions*/
long idum; /*random seed -907 -507139 -39715 -1020307 -97643*/
double dT, delta_s_const, D, dx = GridScale, dy = GridScale, x_loc, y_loc;	/*simulation time step [min], particle spatial step size [cm], diffusion coeff [cm2/min], X and Y cell size, X and Y particle current position */
double leftover_dir, leftover_dist, leftover_time; /*advective velocity direction, "redundant" distance, passed after t=dT, "redundant" time*/
double **K_field, **H_field, **porosity_field; /*K, H, porosity porosity fields*/
//double leftover_dist_diff, double leftover_dir_diff; /*CURRENTLY UNUSED*/

void reactive_particle_tracker(void)
{
	/***************************** DECLARATIONS / INITIALIZATIONS (LOCAL SCOPE VARIABLES) ********************************/   
	
	void move(void);
	
	/*scalars*/
	int i, j, h, max_parts_in_cell, ii, jj, hh; /*indices*/
	int x_cellnum, y_cellnum, x_cellnum_prev, y_cellnum_prev, a_part_num, b_part_num; /*current / previous cell location numbers, particle ID numbers*/
	int num_of_react, output_freq = Frequency_of_output, Na_per_PVT, react_include, mid_point; /*number of time steps / reaction events per PVT, output frequency (number of time steps between successive outputs), total number of injected A parts per PVT, include reaction, inject all at midpoint*/
	int dNa_inj, Na_max, Nb_max, num_of_frames; /*particle injection rate, Max. allowed cumulative number of A/B particles in the field, num of frames (time steps) until termination*/
	int Na_curr = 0, Nb_curr = 0, N_exit_part_a = 0, N_exit_part_b = 0; /*current number of A,B parts, cumulative number of A, B parts that crossed outlet*/
	int write_file_id = 0, time_step_ind = 0, Nsteps, max_por_flag = 0, min_por_flag = 0; /*output frame index, simulation time step index, number of steps to divide dx into, max/min porosity reached flag*/
	int prec_flag, pH_calc_log, calc_dT_straight, particle_track, inj_flag = 1, diss_only; /*precipitation flag, calc pH: (1) log or (2) lin interp, calc dT: (1) straight from mean k or (2) read from file, 0 - do not print out particle trajectories, 1 - print around a selected value, 2 - print the whole path! */

	float temp_head_scan, temp_k_scan, dummy_head; /*temporary vars to read from file*/

	double total_molar_amount, part_molar_amount, caco3_molar_vol = CaCO3_MV; /* number of H+ moles to fill the pore volume, assuming pH_in, how many mols of H+ in 1 A particle (RECALL THAT 1 A PARTICLE BASICALLY CONTAINS 2 H+), [cm3/mol], CaCO3 specific molar volume [cm3/mol]*/
	double qx_cum, dh_dx; /*cumulative Darcy flow at the inlet [cm/min], hyd. head gradient*/
	double porosity_lower = Porosity_lower_bound, porosity_upper = Porosity_upper_bound; /*porosity limits for reaction*/
	double pH_in, pH_res, pH_loc, theta0, k0; /*incoming, resident, local water pH, initial porosity, hyd. conductivty for the case of init. homogeneous field*/
	double field_length, field_height, conc_h2co3_loc, conc_h2co3_eq; /*field length, height [cm], local non-equilibrium fractional amount of H2CO3 in a cell, fractional (molar) amount of dissolved H2CO3 in equilibrium (interpolation from the supplementary plot)*/
	double mean_k, cum_k, mean_vel, head_grad, temp_k, temp_porosity; /* mean / cum hyd. conductivity, mean flow velocity, head gradient over field, temporary hyd. cond. / porosity vars*/
	double Na_ph_in, sim_time = 0.0, react_acc, PVT, D_a_coeff, D_b_coeff, time_watch, epsilon = 0.05; /*num of H+ parts in a cell to obtain pH=pH_in there, simulation time, reaction acceleration parameter, temporal watch value*/

	/*arrays*/
	int    *inlet_part; /*injected particles per inlet cell*/
	double *qx, *BT_times_a, *BT_times_b; /*inlet Darcy flux, A&B parts breakthrough times*/
	double *Xa, *Ya, *Xb, *Yb; /* x, y-coordinates of each particle*/
	double *leftover_dir_array_a, *leftover_dir_array_b; /* "leftover" advective jump direction */
	double *leftover_time_array_a, *leftover_time_array_b; /* "leftover" advective jump time duration */
	double *leftover_dist_array_a, *leftover_dist_array_b; /* "leftover" advective jump distance */
	
	int **cum_cell_vis_a, **cum_cell_vis_b; /*cumulative A,B particle visitation map*/
	int **curr_cell_partnum_a, **curr_cell_partnum_b, **cum_cell_diss_react, **cum_cell_prec_react; /*current num of A,B parts in each cell, cumulative diss/prec reaction counter in cell*/

	double ***id_array_a, ***id_array_b; /*id arrays for A, B parts currently in each cell (should be int, but there's no i3tensor function in NR)*/
	//double *leftover_diff_dist_array_a, *leftover_diff_dir_array_a; /*CURRENTLY UNUSED*/
    //double *leftover_diff_dist_array_b, *leftover_diff_dir_array_b;
	
	/******************************* INPUT PARAMS - READ FROM FILE ***********************************/
	char Input_file_name[50];
	FILE* Input_file; 
	sprintf(Input_file_name, "Input_file.txt");
	
	Input_file = fopen(Input_file_name, "r");
	(void)fscanf(Input_file, "%*s %*s %*s %*s %*s %*s %*s\n"); 
	(void)fscanf(Input_file, "%d %d %lf %ld %d %lf %lf\n", &numcellX, &numcellY, &head_grad, &idum, &Na_per_PVT, &PVT, &react_acc);
	(void)fscanf(Input_file, "%*s %*s %*s %*s %*s %*s %*s\n"); 
	(void)fscanf(Input_file, "%d %d %lf %lf %d %d %d\n", &max_parts_in_cell, &num_of_react, &D_a_coeff, &D_b_coeff, &Nsteps, &pH_calc_log, &calc_dT_straight);
	(void)fscanf(Input_file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s\n");
	(void)fscanf(Input_file, "%lf %lf %d %lf %lf %lf %d %d %d\n", &theta0, &k0, &react_include, &pH_in, &pH_res, &time_watch, &particle_track, &mid_point, &diss_only);
	fclose(Input_file);
    /*********************************************************************************/

	field_length = numcellX * dx;
	field_height = numcellY * dy;

	Nb_max = Na_per_PVT * 10; /*Max. allowed cumulative number of A (H+) particles in the field*/
	Na_max = Na_per_PVT * 10; /*Max. allowed cumulative number of B (H2CO) particles in the field*/
	num_of_frames = num_of_react * 10; /*TERMINATION TIME = 10 PVT*/
	
	total_molar_amount = theta0 * field_height * field_length * 1e-3 * pow(10.0, -pH_in); /* =4.55e-5 [mol], number of H+ moles to obtain pH_in in the available pore volume*/
	part_molar_amount = total_molar_amount / (double)Na_per_PVT; /* [mol/part], how many mols of H+ in 1 A particle, including reaction acceleration (RECALL THAT 1 A PARTICLE BASICALLY CONTAINS 2 H+)*/
	Na_ph_in = (double)Na_per_PVT / ((double)numcellX * (double)numcellY); /* number of A (2*H+) parts in a cell to obtain pH=3.5 there*/
	delta_s_const = GridScale / (double)Nsteps;

	/*array definitions*/
	/*cell numbers go from 1:Nx, 1:Ny, particle numbers go from 1,Na_max*/
	Xa                    = dvector(1, Na_max);
	Ya                    = dvector(1, Na_max);
	leftover_dist_array_a = dvector(1, Na_max);
	leftover_time_array_a = dvector(1, Na_max);
	leftover_dir_array_a  = dvector(1, Na_max);
	BT_times_a            = dvector(1, Na_max);
	
	cum_cell_vis_a        = imatrix(1, numcellX, 1, numcellY);
	curr_cell_partnum_a   = imatrix(1, numcellX, 1, numcellY);
	
	id_array_a            = d3tensor(1, numcellX, 1, numcellY, 1, max_parts_in_cell); 

	//leftover_diff_dist_array_a = dvector(1, Na_max); /*CURRENTLY UNUSED*/
	//leftover_diff_dir_array_a  = dvector(1, Na_max);

	Xb                    = dvector(1, Nb_max);
	Yb                    = dvector(1, Nb_max);
	leftover_dist_array_b = dvector(1, Nb_max);
	leftover_time_array_b = dvector(1, Nb_max);
	leftover_dir_array_b  = dvector(1, Nb_max);
	BT_times_b            = dvector(1, Nb_max);
	
	cum_cell_vis_b        = imatrix(1, numcellX, 1, numcellY);
	curr_cell_partnum_b   = imatrix(1, numcellX, 1, numcellY);
	
	id_array_b            = d3tensor(1, numcellX, 1, numcellY, 1, max_parts_in_cell);

	//leftover_diff_dist_array_b = dvector(1, Na_max); /*CURRENTLY UNUSED*/
	//leftover_diff_dir_array_b  = dvector(1, Na_max);
	
	inlet_part          = ivector(1, numcellY);
	
	qx                  = dvector(1, numcellY);
	
	cum_cell_diss_react = imatrix(1, numcellX, 1, numcellY);
	cum_cell_prec_react = imatrix(1, numcellX, 1, numcellY);
	
	K_field             = dmatrix(1, numcellX, 1, numcellY);
	porosity_field      = dmatrix(1, numcellX, 1, numcellY);  
	H_field             = dmatrix(1, numcellX + 1, 1, numcellY + 1);
	/* OK 24-6-2023 */

	/*array initialization*/
	for (i = 1; i <= Na_max; i++) { BT_times_a[i] = 0.0; }
	for (i = 1; i <= Nb_max; i++) { BT_times_b[i] = 0.0; }

	for (h = 1; h <= numcellX; h++)
	{
		for (i = 1; i <= numcellY; i++)
		{
			cum_cell_vis_b[h][i] = 0;
			cum_cell_vis_a[h][i] = 0;
			cum_cell_diss_react[h][i] = 0;
			cum_cell_prec_react[h][i] = 0;
		}
	}

	/*filename declaration*/
	char total_part_a_filename[30], total_part_b_filename[30];
	char particle_BT_times_a_filename[30], particle_BT_times_b_filename[30];
	char cum_cell_vis_a_filename[30], cum_cell_vis_b_filename[30];
	char curr_cell_partnum_a_filename[30], curr_cell_partnum_b_filename[30];
	char react_diss_file_filename[30], react_prec_file_filename[30];
	char K_filename[30], K_temp_filename[30];
	char H_filename[30], H_temp_filename[30];
	char Por_filename[30], Log_filename[30], inlet_part_filename[30];
	char particle_tracker_apart_filename[30], particle_tracker_bpart_filename[30];  /*particle tracker filenames*/

	/*FILES*/
	FILE *Total_part_a_file, *Total_part_b_file; /*Number of a,b parts in the field in the beginning of each time step dT*/
	FILE *particle_BT_times_a_file, *particle_BT_times_b_file; /*Number of a,b parts that have exited the field each time step dT*/
	FILE *cum_cell_vis_a_file, *cum_cell_vis_b_file; /*particle A,B cumulative cell visitations*/
	FILE *curr_cell_partnum_a_file, *curr_cell_partnum_b_file; /*current particle number in cell*/
	FILE *react_diss_file, *react_prec_file; /*number of diss/prec reactions occured in each cell*/
	FILE *K_file, *K_file_temp; /*K field frames output, dummy K field file for the flow solver*/
	FILE *H_file, *H_file_temp; /*H field frames output, dummy H field file to read the flow solver output*/
	FILE *Por_file, *Log_file, *inlet_part_file; /*porosity field each time step, log file, inlet part file*/
	FILE *particle_tracker_apart_file, *particle_tracker_bpart_file;  /*particle tracker filenames*/

	sprintf(total_part_a_filename, "Total_part_a.dat");
	sprintf(total_part_b_filename, "Total_part_b.dat");
	sprintf(particle_BT_times_a_filename, "particle_BT_times_a.dat");
	sprintf(particle_BT_times_b_filename, "particle_BT_times_b.dat");
	sprintf(K_filename, "OrigTRS0001.txt");
	sprintf(K_temp_filename, "TRS0001.txt");
	sprintf(H_filename, "OrigFIP0001.txt");
	sprintf(H_temp_filename, "FIP0001.txt");
	sprintf(Log_filename, "Logfile.txt");
	
	/******************************************** END DECLARATIONS (LOCAL SCOPE VARIABLES) ********************************/
	
	printf("Program running.\n");

	Log_file = fopen(Log_filename, "w");
	fprintf(Log_file, " Nx Ny head_grad   idum      Na_per_PVT      PVT        react_acc\n");
	fprintf(Log_file, " %d %d %4.3e  %ld       %d      %lf    %lf\n", numcellX, numcellY, head_grad, idum, Na_per_PVT, PVT, react_acc);
	fprintf(Log_file, " max_parts_in_cell   num_of_react   D_a_coeff   D_b_coeff     Nsteps    pH_calc_log    calc_dT_straight\n");
	fprintf(Log_file, " %d                 %d         %4.3e       %4.3e   %d            %d            %d\n", max_parts_in_cell, num_of_react, D_a_coeff, D_b_coeff, Nsteps, pH_calc_log, calc_dT_straight);
	fprintf(Log_file, "Theta0   k0   react_include  pH_in   pH_res   time_watch   particle_track   mid_point    diss_only\n");
	fprintf(Log_file, "%1.2lf    %2.4lf    %d         %1.2lf       %1.2lf     %2.3lf            %d              %d       %d\n", theta0, k0, react_include, pH_in, pH_res, time_watch, particle_track, mid_point, diss_only);
	/********************* DEFINE INITIAL K, H, POROSITY FIELDS (HOMOGENEOUS / HETEROGENEOUS) ***************************************/

	K_file = fopen(K_filename, "r"); /*look for the initial (heterogeneous) field OrigTRS0001.txt*/
	K_file_temp = fopen(K_temp_filename, "w"); /* file to define initial(homogeneous / heterogeneous) conductivity, write to TRS0001.txt for the flow solver*/

	if (K_file == NULL) /*if no OrigTRS0001.txt present*/
	{
		printf("Begin with homogeneous field.\n");
		fprintf(Log_file, "Begin with homogeneous field.\n");
		
		/*create OrigTRS0001.txt, write initial field*/
		K_file = fopen(K_filename, "w");

		for (h = 1; h <= numcellX; h++) /*define uniform K field, write to file*/
		{
			for (i = 1; i <= numcellY; i++)
			{
				K_field[h][i] = k0;
				porosity_field[h][i] = theta0;
				fprintf(K_file_temp, " %8.11f\n ", K_field[h][i]);
				fprintf(K_file, " %8.11f\n ", K_field[h][i]);
			}
		}
		
	}
	else /*scan data from OrigTRS0001.txt*/
	{
		printf("Begin with heterogeneous field.\n");
		fprintf(Log_file, "Begin with heterogeneous field.\n");

		for (h = 1; h <= numcellX; h++) 
		{
			for (i = 1; i <= numcellY; i++)
			{
				(void)fscanf(K_file, "%f\n", &temp_k_scan);
				K_field[h][i] = temp_k_scan;
				fprintf(K_file_temp, " %8.11lf\n ", K_field[h][i]);
				porosity_field[h][i] = theta0; /*assume constant porosity at t=0*/
			}
		}

	}

	fclose(K_file_temp);
	fclose(K_file);
	fclose(Log_file);

	system("H0C100_PAPER2.exe"); /* run flow solver on TRS0001.txt */

	H_file = fopen(H_filename, "w");
	H_file_temp = fopen(H_temp_filename, "r");    
	if ( H_file_temp == NULL) { printf("Error: can't open FIP0001.txt file.\n"); }

	/*read H field from file FIP0001.txt just obtained from the flow solver, print to OrigFIP0001.txt*/
	for (h = 1; h <= numcellX+1; h++)
	{
		for (i = 1; i <= numcellY+1; i++)
		{
			(void)fscanf(H_file_temp, "%f", &temp_head_scan);
			H_field[h][i] = temp_head_scan;

			if (h == 1 && i == 1) { dummy_head = H_field[h][i]; }
			H_field[h][i] = head_grad / dummy_head * H_field[h][i];

			fprintf(H_file, " %8.11f\n ", H_field[h][i]);
		}
	}

	fclose(H_file_temp);
	fclose(H_file);
	
//	dummy_head = H_field[1][1];

	/********************* END OF INITIAL K, H, POROSITY FIELDS INITIALIZATION ***************/
		
	/************************ CALCULATE PARTICLE INJECTION RATE, ETC. *******************************/
	if(calc_dT_straight == 1)
	{
		cum_k = 0;
		for (h = 1; h <= numcellX; h++) { for (i = 1; i <= numcellY; i++) { cum_k += exp(K_field[h][i]); } }
		mean_k = cum_k / ((double)numcellX * (double)numcellY);
		mean_vel = mean_k * (head_grad / field_length) / theta0;
		PVT = field_length / mean_vel;
		dT = PVT / (double)num_of_react; /*time step between reaction events*/
	}
	else if (calc_dT_straight == 0) { dT = PVT / (double)num_of_react; } /*calculate reaction time step dT*/

	dNa_inj = (int)floor((double)Na_per_PVT / (double)num_of_react); /*num of particles to inject in dT, calculated so that N=Na_per_PVT will be reached at Tpv: Ndot = Na_per_PVT/Tpv.*/
	if (mid_point == 1) { dNa_inj = Na_per_PVT; }
	/*for spontaneous injection (Langevin verification), choose midpoint = 1, react_include = 0*/
    /**************************************************************************************************************/

    /********************************** BEGIN SIMULATION (RUN ON TIME AXIS) ***********************************/
	while (time_step_ind < num_of_frames)
    {
		time_step_ind += 1; /*current time step index*/
		sim_time = time_step_ind * dT; /*current time [min]*/

		printf("sim_time %8.11f, %8.5f PVT\n", sim_time, sim_time/PVT);

		if (min_por_flag == 1) { printf("Min. Porosity has been reached.\n"); }
		if (max_por_flag == 1) { printf("Max. Porosity has been reached.\n"); }

		/*INITIALIZE TO ZERO EVERY TIME STEP (AFTER MOVING A,B PARTS THEIR LOCATION IS DEFINED ANEW)*/
	    for (h = 1; h <= numcellX; h++)
	    {
		    for (i = 1; i <= numcellY; i++)
			{
				curr_cell_partnum_a[h][i] = 0;
				curr_cell_partnum_b[h][i] = 0;
				//cum_cell_diss_react[h][i] = 0;
			    //cum_cell_prec_react[h][i] = 0;

				for (j = 1; j <= max_parts_in_cell; j++) /*UPDATE! Previously was not set to zero each time step*/
				{
					id_array_a[h][i][j] = 0.0;
					id_array_b[h][i][j] = 0.0;
				}
			}
		}

		/*flux-weighting for particle injection*/
		qx_cum = 0.0;
		for (i = 1; i <= numcellY; i++) /*calculating vol. flow in inlet cells. Mixed derivative is used for interpolation within cell, not necessary here.*/
		{
			dh_dx = (H_field[2][i] - H_field[1][i]) / dx;
			qx[i] = -exp(K_field[1][i]) * dh_dx * (dy * 1.0);  /*inlet cell vol. flow rate [cm3/min/cm] = Darcy flux * area per unit in-plane depth.*/
			qx_cum = qx_cum + qx[i]; /*total inlet flow*/
		}

			/*INJECT NEW A PARTS AS LONG AS MAX. ALLOWED NUMBER OF A PARTS HAS NOT BEEN REACHED YET*/
			if (Na_curr + dNa_inj < Na_max && inj_flag == 1) /*condition to prevent A arrays overflow*/
			{
				/*distribute each newly injected A particle*/
				for (i = 1; i <= numcellY; i++)
				{
					inlet_part[i] = (int)round(dNa_inj * qx[i] / qx_cum); //(int)round(Na_per_PVT * qx[i] / qx_cum);   
					if (mid_point == 1) { inlet_part[i] = dNa_inj * (i == numcellY/2); }
					/*distribute the next inlet_part[i] particles in the i-th cell*/
					for (j = Na_curr + 1; j <= Na_curr + inlet_part[i]; j++)
					{
						/*init. X location*/
						Xa[j] = 0.0; /*was 1.0e-3 previously, needs check for negative value!*/
						
						if (mid_point == 0) { Ya[j] = ( (double)i - ran2(&idum) ) * dy; } /* flux-weight to distribute in all inlet cells - uniformly distributed*/
						else { Ya[j] = (numcellY / 2) * dy; } /*inject all in the midcell*/
						
						/*out of bounds, reflective BC*/
						if (Ya[j] > field_height) { Ya[j] = field_height + (field_height - Ya[j]); } 
						else if (Ya[j] < 0.0) { Ya[j] = -Ya[j]; }
						
						if (Xa[j] < 0.0) { Xa[j] = -Xa[j]; }

						x_cellnum = (int)ceil(Xa[j] / dx); /*cell location, cell numbers go from 1:Nx, 1:Ny*/
						y_cellnum = (int)ceil(Ya[j] / dy);
						
						/*fix out of bounds*/
						if (y_cellnum > numcellY) { y_cellnum = numcellY; } 
						if (x_cellnum > numcellX) { x_cellnum = numcellX; }

						if (y_cellnum < 1) { y_cellnum = 1; }
						if (x_cellnum < 1) { x_cellnum = 1; }
						
						cum_cell_vis_a[x_cellnum][y_cellnum] += 1;

						leftover_dir_array_a[j] = 0.0;
						leftover_dist_array_a[j] = 0.0;
						leftover_time_array_a[j] = 0.0;

						//leftover_diff_dist_array_a[j] = 0.0;
						//leftover_diff_dir_array_a[j] = 0.0;
					}
					Na_curr = Na_curr + inlet_part[i];
				}
				if (mid_point == 1) { inj_flag = 0; }
			}
			//else { printf("Max. number of A particles has been reached.\n"); }

		Total_part_a_file = fopen(total_part_a_filename, "a"); /*Na_curr IN THE BEGINNING OF A STEP - write to file*/
		fprintf(Total_part_a_file, "%8.11f %i\n", sim_time, Na_curr);
		fclose(Total_part_a_file);

		Total_part_b_file = fopen(total_part_b_filename, "a"); /*Nb_curr IN THE BEGINNING OF A STEP - write to file*/
		fprintf(Total_part_b_file, "%8.11f %i\n", sim_time, Nb_curr);
		fclose(Total_part_b_file);
		/*OK 27-6-2023.*/

		/*move A (H+) particles*/
		D = D_a_coeff;
		i = 1;
		while (i <= Na_curr)
		{
			x_loc = Xa[i];
			y_loc = Ya[i];
			leftover_dist = leftover_dist_array_a[i];
			leftover_time = leftover_time_array_a[i];
			leftover_dir = leftover_dir_array_a[i];

			//leftover_dist_diff = leftover_diff_dist_array_a[i];
			//leftover_dir_diff = leftover_diff_dir_array_a[i];

			x_cellnum_prev = (int)ceil(x_loc / dx); /*current cell, ceil rounds to the nearest integer from above. Cell numbers go from 1:Nx, 1:Ny, thus need to fix for spatial_X=0 (cannot be since all parts are given initial X coord witihn the 1st cell)*/
			y_cellnum_prev = (int)ceil(y_loc / dy);

      		move();

			Xa[i] = x_loc;
			Ya[i] = y_loc;
			leftover_dist_array_a[i] = leftover_dist;
			leftover_time_array_a[i] = leftover_time;
			leftover_dir_array_a[i] = leftover_dir;

			//leftover_diff_dist_array_a[i] = leftover_dist_diff;
			//leftover_diff_dir_array_a[i] = leftover_dir_diff;

			x_cellnum = (int)ceil(x_loc / dx); /*new cell location, cell numbers go from 1:Nx, 1 : Ny*/
			y_cellnum = (int)ceil(y_loc / dy);

			if (x_cellnum <= numcellX && y_cellnum <= numcellY && x_cellnum > 0 && y_cellnum > 0) /*if within bounds*/
			{

				if (x_cellnum != x_cellnum_prev || y_cellnum != y_cellnum_prev) { cum_cell_vis_a[x_cellnum][y_cellnum] += 1; } /*cumulative particle number in the cell*/
				curr_cell_partnum_a[x_cellnum][y_cellnum] += 1; /*current particle number in the cell.*/
				j = curr_cell_partnum_a[x_cellnum][y_cellnum];
				if (j > max_parts_in_cell) { printf("A particle max_parts_in_cell index exceeded!!!\n"); }
				id_array_a[x_cellnum][y_cellnum][j] = (double)i;
			}

			if (Xa[i] > field_length) /*x_cellnum > numcellX, remove the ith particle, stay on the same index as it's replaced by the next one*/
			{
				N_exit_part_a += 1;
				BT_times_a[N_exit_part_a] = sim_time; /*get this more precisely?*/

				memcpy(&(Xa[i]),                     &(Xa[i + 1]),                     (Na_max - i) * sizeof(double));
				memcpy(&(Ya[i]),                     &(Ya[i + 1]),                     (Na_max - i) * sizeof(double));
				memcpy(&(leftover_dist_array_a[i]),  &(leftover_dist_array_a[i + 1]),  (Na_max - i) * sizeof(double));
				memcpy(&(leftover_time_array_a[i]),  &(leftover_time_array_a[i + 1]),  (Na_max - i) * sizeof(double));
				memcpy(&(leftover_dir_array_a[i]),   &(leftover_dir_array_a[i + 1]),   (Na_max - i) * sizeof(double));

			  //memcpy(&(leftover_diff_dist_array_a[i]), &(leftover_diff_dist_array_a[i + 1]), (Na_max - i) * sizeof(double));
			  //memcpy(&(leftover_diff_dir_array_a[i]),  &(leftover_diff_dir_array_a[i + 1]),  (Na_max - i) * sizeof(double));

				Na_curr -= 1; /*no need to update curr_cell_partnum & id arrays here because they didn't include this exited particle from the beginning (they are being redefined every time step)! */
			}
			else { i++; } /*update index - go to the next particle*/

		} /*end of while (i < Na_curr)*/

		particle_BT_times_a_file = fopen(particle_BT_times_a_filename, "a"); /*cumulative number of exited A parts by the current time*/
		fprintf(particle_BT_times_a_file, " %8.11f %i\n", sim_time, N_exit_part_a); /*calc. BT times more accurately???*/
		fclose(particle_BT_times_a_file);

		/*move B (H2CO3) particles (IN THE BEGINNING THERE ARE NO B PARTS BEFORE REACTION TOOK PLACE)*/
		D = D_b_coeff;
		h = 1;
		while (h <= Nb_curr) 
		{
			x_loc = Xb[h]; /*current coords*/
			y_loc = Yb[h];
			leftover_dist = leftover_dist_array_b[h];
			leftover_time = leftover_time_array_b[h];
			leftover_dir = leftover_dir_array_b[h];
			
			//leftover_dist_diff = leftover_diff_dist_array_b[h];
			//leftover_dir_diff = leftover_diff_dir_array_b[h];

			x_cellnum_prev = (int)ceil(x_loc / dx); /*current cell, ceil rounds to the nearest integer from above. Cell numbers go from 1:Nx, 1:Ny, thus need to fix for spatial_X=0 (cannot be since all parts are given initial nonzero X coord witihn the 1st cell)*/
			y_cellnum_prev = (int)ceil(y_loc / dy);
			
			move(); /*advance particles*/
		
			Xb[h] = x_loc; /*new coords*/
			Yb[h] = y_loc;
			leftover_dist_array_b[h] = leftover_dist;
			leftover_time_array_b[h] = leftover_time;
			leftover_dir_array_b[h] = leftover_dir;

			//leftover_diff_dist_array_b[h] = leftover_dist_diff;
			//leftover_diff_dir_array_b[h] = leftover_dir_diff;
			
			x_cellnum = (int)ceil(x_loc / dx); /*new cell, cell numbers go from 1:Nx, 1:Ny*/
			y_cellnum = (int)ceil(y_loc / dy);
			
			if (x_cellnum <= numcellX && y_cellnum <= numcellY && x_cellnum > 0 && y_cellnum > 0) /*if within bounds*/
			{
			   if (x_cellnum != x_cellnum_prev || y_cellnum  != y_cellnum_prev) { cum_cell_vis_b[x_cellnum][y_cellnum] += 1; } /*update particle visitation if particle has jumped to a new cell*/
			   curr_cell_partnum_b[x_cellnum][y_cellnum] += 1; /*REDEFINE curr_cell_partnum_b, NOT WITHIN IF BECAUSE curr_cell_partnum_A,B ARE BEING INITIALIZED TO ZERO EACH STEP!*/ 
			   j = curr_cell_partnum_b[x_cellnum][y_cellnum]; 
			   if (j > max_parts_in_cell) { printf("B particle max_parts_in_cell index exceeded!!!\n"); }
			   id_array_b[x_cellnum][y_cellnum][j] = (double)h;
			}

			if (Xb[h] > field_length) /*if passed the outlet, "h" - B particle index*/
			{ /*in this case the particle is removed from the arrays, the next one adopts its ID number*/
				
				N_exit_part_b += 1;
				BT_times_b[N_exit_part_b] = sim_time; /* BTC time */

				/*Copy block of memory - Copies the values of num bytes from the location pointed to by source directly to the memory block pointed to by destination.*/
				/*void * memcpy ( void * destination, const void * source, size_t num );*/
				/*hth part. has exited -> shift the {h+1,END} array backwards by 1. Effectively remove "h" particle from memory, replace it with "h+1" by moving all array members beginning from "h+1" back by 1*/
				memcpy(&(Xb[h]),                    &(Xb[h + 1]),                    (Nb_max - h) * sizeof(double) );  
				memcpy(&(Yb[h]),                    &(Yb[h + 1]),                    (Nb_max - h) * sizeof(double) );
				memcpy(&(leftover_dist_array_b[h]), &(leftover_dist_array_b[h + 1]), (Nb_max - h) * sizeof(double) );
				memcpy(&(leftover_time_array_b[h]), &(leftover_time_array_b[h + 1]), (Nb_max - h) * sizeof(double) );
				memcpy(&(leftover_dir_array_b[h]),  &(leftover_dir_array_b[h + 1]),  (Nb_max - h) * sizeof(double) );

				//memcpy(&(leftover_diff_dist_array_b[h]), &(leftover_diff_dist_array_b[h + 1]), (Na_max - h) * sizeof(double));
				//memcpy(&(leftover_diff_dir_array_b[h]),  &(leftover_diff_dir_array_b[h + 1]),  (Na_max - h) * sizeof(double));

				Nb_curr -= 1;  /*update current number of B parts*/
			}
			else { h++; } /*advance next particle*/

		} /*end of while (h < Nb_curr)*/

		particle_BT_times_b_file = fopen(particle_BT_times_b_filename, "a");
		fprintf(particle_BT_times_b_file, " %8.11f %i\n", sim_time, N_exit_part_b);
		fclose(particle_BT_times_b_file);

		//------------------------------------------------------------
		//		Reaction part
		//------------------------------------------------------------
		/*After all A,B particles have advanced - perform reaction*/
		if (react_include == 1)
		{

			for (h = 1; h <= numcellX; h++) /*run on all cells*/
			{
				for (i = 1; i <= numcellY; i++)
				{

					if (Na_curr >= Na_max) { printf("Max. number of A particles has been reached.\n"); }
					if (Nb_curr >= Nb_max) { printf("Max. number of B particles has been reached.\n"); }

					if (porosity_field[h][i] >= porosity_upper)
					{
						if (max_por_flag == 0)
						{
							printf("Max. Porosity has been reached.\n");
							Log_file = fopen(Log_filename, "a");
							fprintf(Log_file, "Max. Porosity has been reached at %8.5f PVT\n", sim_time / PVT);
							fclose(Log_file);

							max_por_flag = 1;
						}
					}

					if (porosity_field[h][i] <= porosity_lower) 
					{
						if (min_por_flag == 0)
						{
							printf("Min. Porosity has been reached.\n");
							Log_file = fopen(Log_filename, "a");
							fprintf(Log_file, "Min. Porosity has been reached at %8.5f PVT\n", sim_time / PVT);
							fclose(Log_file);

							min_por_flag = 1;
						}
					}

					/*CALCULATE LOCAL CONCENTRATIONS*/
					/*(molar) fractional amount of dissolved H2CO3 as a function of pH in a cell [H2CO3][mol/cm3] / [total C] (which is [H2CO3]+[HCO3-]+[CO3(2-)]) in equilibrium (interpolation of the supplementary plot)*/
					/*current pH in cell, linear interpolation*/
					
					if (pH_calc_log == 0) { pH_loc = pH_res - (pH_res - pH_in) * (double)curr_cell_partnum_a[h][i] / Na_ph_in; }  
					else if (pH_calc_log == 1) { pH_loc = -log10( part_molar_amount*(double)curr_cell_partnum_a[h][i] / (porosity_field[h][i]*dx*dy*1.0e-3) + pow(10.0, -pH_res) ); }
					
					conc_h2co3_eq = 0.5 * (1.0 - tanh(pH_loc - 6.35)); /*GOOD APPROX. TO THE H2CO3(PH) CURVE FROM THE BJERRUM PLOT*/
					conc_h2co3_loc = (double)curr_cell_partnum_b[h][i] / Na_ph_in; /*fractional amount of dissolved H2CO3 currently (NOT IN EQUILIBRIUM). SHOULD BE MULTIPLIED BY 1/2 BECAUSE OF DIFFERENCE IN MOLAR WEIGHTS!*/
					/*BASED ON ASSUMPTION THAT TOTAL CARBONATE AMOUNT [MOL] IN A CELL EQUALS TO [H+] [MOL] TO OBTAIN PH=3.5 IN A CELL*/
					/*ASSUMING CT CONSTANT THROUGHOUT THE FIELD AND IN TIME (WELL-MIXED), THIS CAN BE ACHIEVED BY NOTING THAT FOR PH=3.5, CT = [H2CO3]. ASSUMING [CA2+], WE THEN HAVE EVERYTHING FROM KEQ*/

					/*DISSOLUTION REACTION - while in the cell current value of [H2CO3]/[total C] is SMALLER than should be in equilibrium, current number of B parts in the field is less than max. allowed,
					porosity is less than max. allowed and at least 1 A part is present*/
					prec_flag = 0;
					//while (conc_h2co3_eq > conc_h2co3_loc && Nb_curr < Nb_max && porosity_field[h][i] < porosity_upper && curr_cell_partnum_a[h][i]>0)
					while (conc_h2co3_eq > conc_h2co3_loc && fabs(conc_h2co3_eq - conc_h2co3_loc)/ conc_h2co3_eq > epsilon && Nb_curr < Nb_max && porosity_field[h][i] < porosity_upper && curr_cell_partnum_a[h][i] > 1)
					{
						/*FIND ID OF THE A PARTICLE TO UNDERGO REACTION*/
						j = curr_cell_partnum_a[h][i]; /*current number of H+ parts in a cell -> particle index in id_array_a,b arrays, begins from 1*/
						a_part_num = (int)floor(id_array_a[h][i][j]); /*ID number of an A part to react - the last part to arrive in cell, begins from 1!*/
						if (a_part_num < 1 || a_part_num > Na_curr) 
						{ printf("a_part_num out of bounds.\n"); }

						/*UPDATE MOMENTOUS COUNTERS*/
						Na_curr -= 1; /*A part disappears & B appears in dissolution reaction (RECALL THAT THEIR MOLAR RATIO is [A]=x, [B]=x/2)*/
						curr_cell_partnum_a[h][i] -= 1;  /*A part disappears & B appears in dissolution reaction*/
						Nb_curr += 1;
						curr_cell_partnum_b[h][i] += 1;

						/*INITIALIZE TRANSPORT ARRAYS OF THE NEWLY APPEARED B PARTICLE*/
						Xb[Nb_curr] = Xa[a_part_num];
						Yb[Nb_curr] = Ya[a_part_num];
						leftover_dist_array_b[Nb_curr] = leftover_dist_array_a[a_part_num];
						leftover_time_array_b[Nb_curr] = leftover_time_array_a[a_part_num];
						leftover_dir_array_b[Nb_curr] = leftover_dir_array_a[a_part_num];
					  //leftover_diff_dist_array_b[Nb_curr] = leftover_diff_dist_array_a[a_part_num];
					  //leftover_diff_dir_array_b[Nb_curr] = leftover_diff_dir_array_a[a_part_num];

						/*UPDATE ID_ARRAY_A TO ACCOUNT FOR THE DISAPPEARED A PARTICLE.*/
						id_array_a[h][i][j] = -1.0;

						for (hh = 1; hh <= numcellX; hh++)
						{
							for (ii = 1; ii <= numcellY; ii++)
							{
								for (jj = 1; jj <= curr_cell_partnum_a[hh][ii]; jj++) /*previously read as for (jj = 1; jj <= max_parts_in_cell; jj++) */
								{
									if (id_array_a[hh][ii][jj] > (double)a_part_num) { id_array_a[hh][ii][jj] -= 1.0; }
								}
							}
						}
						
						/*UPDATE ID_ARRAY_B TO ACCOUNT FOR THE NEWLY APPEARED B PARTICLE*/
						j = curr_cell_partnum_b[h][i];
						if (j >= max_parts_in_cell) { printf("B particle max_parts_in_cell index exceeded!!!\n"); }
						id_array_b[h][i][j] = (double)Nb_curr;

						/*remove the reacted A part from movement arrays by shifting backwards one cell and running over the just reacted part*/
						memcpy(&(Xa[a_part_num]),                    &(Xa[a_part_num + 1]),                    (Na_max - a_part_num) * sizeof(double));
						memcpy(&(Ya[a_part_num]),                    &(Ya[a_part_num + 1]),                    (Na_max - a_part_num) * sizeof(double));
						memcpy(&(leftover_dist_array_a[a_part_num]), &(leftover_dist_array_a[a_part_num + 1]), (Na_max - a_part_num) * sizeof(double));
						memcpy(&(leftover_time_array_a[a_part_num]), &(leftover_time_array_a[a_part_num + 1]), (Na_max - a_part_num) * sizeof(double));
						memcpy(&(leftover_dir_array_a[a_part_num]),  &(leftover_dir_array_a[a_part_num + 1]),  (Na_max - a_part_num) * sizeof(double));
					  //memcpy(&(leftover_diff_dist_array_a[a_part_num]), &(leftover_diff_dist_array_a[a_part_num + 1]), (Na_max - a_part_num) * sizeof(double));
					  //memcpy(&(leftover_diff_dir_array_a[a_part_num]), &(leftover_diff_dir_array_a[a_part_num + 1]), (Na_max - a_part_num) * sizeof(double));

						/*update cumulative counters*/
						cum_cell_vis_b[h][i] += 1; /*new B particle was created - update cum. vis. counter*/
						cum_cell_diss_react[h][i] += 1; /*update reaction counter in the cell*/

						/*update porosity & hyd. conductivity - KC*/
						temp_porosity = porosity_field[h][i]; /*porosity before reaction*/
						porosity_field[h][i] += 0.5 * part_molar_amount * caco3_molar_vol / (dx * dy) * react_acc;  /*porosity after reaction. Based on d[Calcite] = -0.5*d[2*H+]. React. acc. to enhance process.*/
						temp_k = exp(K_field[h][i]) * pow(porosity_field[h][i] / temp_porosity, 3) * pow( (1 - temp_porosity) / (1 - porosity_field[h][i]), 2); /*Kozeny-Carman relation OK.*/
						K_field[h][i] = log(temp_k); /*conductivity after reaction*/

						/*update concentrations*/
						if (pH_calc_log == 0) { pH_loc = pH_res - (pH_res - pH_in) * (double)curr_cell_partnum_a[h][i] / Na_ph_in; }
						else if (pH_calc_log == 1) { pH_loc = -log10( part_molar_amount*(double)curr_cell_partnum_a[h][i] / (porosity_field[h][i] * dx * dy * 1.0e-3) + pow(10.0, -pH_res) ); }

						conc_h2co3_eq = 0.5 * (1.0 - tanh(pH_loc - 6.35)); /*recalculate equilibrium & current H2CO3 concentration*/
						conc_h2co3_loc = (double)curr_cell_partnum_b[h][i] / Na_ph_in; /*should be in precipitation as well!!!*/

						prec_flag = 1; /*do not proceed with precipitation after dissolution finished*/
					} /*dissolution reaction finished*/

					/*PRECIPITATION REACTION - while in the cell current value of [H2CO3]/[total C] is larger than should be in equilibrium, current numbr of A parts is less than max. allowed,
					porosity is more than min. allowed, at least 1 B part is present and prec. reaction hasn't occured before in the current time step*/
					//while (conc_h2co3_eq < conc_h2co3_loc && Na_curr < Na_max && porosity_field[h][i] > porosity_lower && curr_cell_partnum_b[h][i] > 0 && prec_flag == 0)
					while (conc_h2co3_eq < conc_h2co3_loc && fabs(conc_h2co3_eq - conc_h2co3_loc) / conc_h2co3_eq > epsilon && Na_curr < Na_max && porosity_field[h][i] > porosity_lower && curr_cell_partnum_b[h][i] > 1 && prec_flag == 0)
					{
						/*FIND ID OF THE B PARTICLE TO UNDERGO REACTION*/
						j = curr_cell_partnum_b[h][i];
						b_part_num = (int)floor(id_array_b[h][i][j]);
						if (b_part_num < 1 || b_part_num > Nb_curr) 
						{ printf("b_part_num out of bounds.\n"); }

						/*UPDATE MOMENTOUS COUNTERS*/
						Na_curr += 1; 
						curr_cell_partnum_a[h][i] += 1;
						Nb_curr -= 1;
						curr_cell_partnum_b[h][i] -= 1;

						/*INITIALIZE TRANSPORT ARRAYS OF THE NEWLY APPEARED A PARTICLE*/
						Xa[Na_curr] = Xb[b_part_num];
						Ya[Na_curr] = Yb[b_part_num];
						leftover_dist_array_a[Na_curr] = leftover_dist_array_b[b_part_num];
						leftover_time_array_a[Na_curr] = leftover_time_array_b[b_part_num];
						leftover_dir_array_a[Na_curr] = leftover_dir_array_b[b_part_num];
					  //leftover_diff_dist_array_a[Na_curr] = leftover_diff_dist_array_b[b_part_num];
					  //leftover_diff_dir_array_a[Na_curr] = leftover_diff_dir_array_b[b_part_num];

						/*UPDATE ID_ARRAY_B TO ACCOUNT FOR THE DISAPPEARED B PARTICLE. MAKE ID_ARRAY_B 1D???*/
						id_array_b[h][i][j] = -1.0;

						for (hh = 1; hh <= numcellX; hh++)
						{
							for (ii = 1; ii <= numcellY; ii++)
							{
								for (jj = 1; jj <= curr_cell_partnum_b[hh][ii]; jj++) /*previously read as for (jj = 1; jj <= max_parts_in_cell; jj++)*/
								{
									if (id_array_b[hh][ii][jj] > (double)b_part_num) { id_array_b[hh][ii][jj] -= 1.0; }
								}
							}
						}

						/*UPDATE ID_ARRAY_A TO ACCOUNT FOR THE NEWLY APPEARED A PARTICLE*/
						j = curr_cell_partnum_a[h][i];
						if (j > max_parts_in_cell) { printf("A particle max_parts_in_cell index exceeded!!!\n"); }
						id_array_a[h][i][j] = (double)Na_curr;

						/*remove the reacted B part from movement arrays by shifting backwards one cell and running over the just reacted part*/
						memcpy(&(Xb[b_part_num]),                    &(Xb[b_part_num + 1]),                    (Nb_max - b_part_num) * sizeof(double));
						memcpy(&(Yb[b_part_num]),                    &(Yb[b_part_num + 1]),                    (Nb_max - b_part_num) * sizeof(double));
						memcpy(&(leftover_dist_array_b[b_part_num]), &(leftover_dist_array_b[b_part_num + 1]), (Nb_max - b_part_num) * sizeof(double));
						memcpy(&(leftover_time_array_b[b_part_num]), &(leftover_time_array_b[b_part_num + 1]), (Nb_max - b_part_num) * sizeof(double));
						memcpy(&(leftover_dir_array_b[b_part_num]),  &(leftover_dir_array_b[b_part_num + 1]),  (Nb_max - b_part_num) * sizeof(double));
					  //memcpy(&(leftover_diff_dist_array_b[b_part_num]), &(leftover_diff_dist_array_b[b_part_num + 1]), (Na_max - b_part_num) * sizeof(double));
					  //memcpy(&(leftover_diff_dir_array_b[b_part_num]), &(leftover_diff_dir_array_b[b_part_num + 1]), (Na_max - b_part_num) * sizeof(double));

						/*update cumulative counters*/
						cum_cell_vis_a[h][i] += 1;
						
						/*update porosity & hyd. cond.*/
						if (diss_only == 0)
						{
							cum_cell_prec_react[h][i] += 1;

							temp_porosity = porosity_field[h][i];
							porosity_field[h][i] -= 0.5 * part_molar_amount * caco3_molar_vol / (dx * dy) * react_acc; /*current react_acc is equivalent to 0.25 * react_acc as in the paper draft version!!!*/
							temp_k = exp(K_field[h][i]) * pow(porosity_field[h][i] / temp_porosity, 3) * pow((1 - temp_porosity) / (1 - porosity_field[h][i]), 2); /*Kozeny-Carman relation OK.*/
							K_field[h][i] = log(temp_k);
						}

						/*update concentrations*/
						if (pH_calc_log == 0) { pH_loc = pH_res - (pH_res - pH_in) * (double)curr_cell_partnum_a[h][i] / Na_ph_in; }
						else if (pH_calc_log == 1) { pH_loc = -log10(pow(10.0, -pH_res) + part_molar_amount * (double)curr_cell_partnum_a[h][i] / (porosity_field[h][i] * dx * dy * 1.0e-3)); }

						conc_h2co3_eq = 0.5 * (1.0 - tanh(pH_loc - 6.35)); /*recalculate equilibrium & current H2CO3 concentration*/
						conc_h2co3_loc = (double)curr_cell_partnum_b[h][i] / Na_ph_in; 
					
					} /*precipitation reaction finished*/

				} /*end of for (i = 1; i <= Head_y_length; i++)*/

			} /*end of run on all cells*/

		} /*end of if react_include == 1*/
		/*************************** END OF REACTION MODULE *************************************/
		/*OK 1-7-2023.*/

		/*if current frame is output frame - create files, write to files, recalc H field*/
		if (time_step_ind >= output_freq && time_step_ind % output_freq == 0)  
		{
			write_file_id += 1;

			sprintf(Por_filename,                "POR%04d.txt",					  write_file_id);
			sprintf(K_filename,                   "ModifiedTRS%04d.txt",			  write_file_id);
			sprintf(H_filename,                   "ModifiedFIP%04d.txt",			  write_file_id);
			sprintf(curr_cell_partnum_a_filename, "curr_cell_partnum_a_file%04d.dat", write_file_id);
			sprintf(curr_cell_partnum_b_filename, "curr_cell_partnum_b_file%04d.dat", write_file_id); 
			sprintf(react_diss_file_filename,     "React_DISS%04d.txt",				  write_file_id);
			sprintf(react_prec_file_filename,     "React_PREC%04d.txt",				  write_file_id);
			sprintf(cum_cell_vis_a_filename,      "cum_cell_vis_a_file%04d.dat",	  write_file_id);
			sprintf(cum_cell_vis_b_filename,      "cum_cell_vis_b_file%04d.dat",	  write_file_id);
			
			Por_file				 = fopen(Por_filename, "w");
			K_file					 = fopen(K_filename, "w");
			H_file					 = fopen(H_filename, "w");
			curr_cell_partnum_a_file = fopen(curr_cell_partnum_a_filename, "w");
			curr_cell_partnum_b_file = fopen(curr_cell_partnum_b_filename, "w"); 
			react_diss_file			 = fopen(react_diss_file_filename, "w");
			react_prec_file			 = fopen(react_prec_file_filename, "w");
			cum_cell_vis_a_file		 = fopen(cum_cell_vis_a_filename, "w");
			cum_cell_vis_b_file		 = fopen(cum_cell_vis_b_filename, "w");
			
			K_file_temp              = fopen(K_temp_filename, "w");
				
			for (h = 1; h <= numcellX; h++) /*run on all cells*/
			{
				for (i = 1; i <= numcellY; i++)
				{
					fprintf(Por_file,				  " %8.11f\n ", porosity_field[h][i]);
					fprintf(K_file,					  " %8.11f\n ", K_field[h][i]);
					fprintf(curr_cell_partnum_a_file, " %i\n ",     curr_cell_partnum_a[h][i]);
					fprintf(curr_cell_partnum_b_file, " %i\n ",     curr_cell_partnum_b[h][i]); 
					fprintf(react_diss_file,		  " %i\n ",     cum_cell_diss_react[h][i]);
					fprintf(react_prec_file,		  " %i\n ",     cum_cell_prec_react[h][i]);
					fprintf(cum_cell_vis_a_file,	  " %i\n ",     cum_cell_vis_a[h][i]);
					fprintf(cum_cell_vis_b_file,	  " %i\n ",     cum_cell_vis_b[h][i]);
					
					fprintf(K_file_temp,              " %8.11f\n ", K_field[h][i]);
				}
			}

			/*WRITE INLET PART TO FILE*/
			sprintf(inlet_part_filename, "InletPart%04d.txt", write_file_id);
			inlet_part_file = fopen(inlet_part_filename, "w");

			cum_k = 0.0;
			for (i = 1; i <= numcellY; i++)
			{
				cum_k += inlet_part[i];
				fprintf(inlet_part_file, "%8.11f %8.11f %8.11f %i\n", (double)qx[i]/qx_cum, exp(K_field[1][i]), H_field[1][i]- H_field[2][i], inlet_part[i]);
			}
			fprintf(inlet_part_file, "%i\n", (int)cum_k);
			fclose(inlet_part_file);

			/*WRITE SPATIAL PARTICLE TRACKER TO FILE*/
			if ( (particle_track == 2) || (particle_track == 1 && sim_time >= 0.9 * time_watch && sim_time <= 1.1 * time_watch) ) /*if the particle time has just passed the Time watch monitor*/
			{
				sprintf(particle_tracker_apart_filename, "particle_tracker_apart%04d.txt", write_file_id);
				particle_tracker_apart_file = fopen(particle_tracker_apart_filename, "w");
				for (i = 1; i <= Na_curr; i++) { fprintf(particle_tracker_apart_file, " %i %15.8f %15.8f %15.8f\n", i, sim_time, Xa[i], Ya[i]); }

				sprintf(particle_tracker_bpart_filename, "particle_tracker_bpart%04d.txt", write_file_id);
				particle_tracker_bpart_file = fopen(particle_tracker_bpart_filename, "w");
				for (i = 1; i <= Nb_curr; i++) { fprintf(particle_tracker_bpart_file, " %i %15.8f %15.8f %15.8f\n", i, sim_time, Xb[i], Yb[i]); }
				
				fclose(particle_tracker_apart_file);
				fclose(particle_tracker_bpart_file);
			}

			fclose(Por_file);
			fclose(K_file);
			fclose(curr_cell_partnum_a_file);
			fclose(curr_cell_partnum_b_file);
			fclose(react_diss_file);
			fclose(react_prec_file);
			fclose(cum_cell_vis_a_file);
			fclose(cum_cell_vis_b_file);
			fclose(K_file_temp);

			/*run flow solver*/
			system("H0C100_PAPER2.exe");

			/*read H field from file FIP0001.txt just obtained from the flow solver, print to ModifiedFIP_.txt*/
			H_file_temp = fopen(H_temp_filename, "r");
			if (H_file_temp == NULL) { printf("Error: can't open FIP0001.txt file.\n"); }

			for (h = 1; h <= numcellX + 1; h++)
			{
				for (i = 1; i <= numcellY + 1; i++)
				{
					(void)fscanf(H_file_temp, "%f\n", &temp_head_scan);
					H_field[h][i] = head_grad / dummy_head * temp_head_scan;
					fprintf(H_file, " %8.11f\n ", H_field[h][i]);
				}
			}
			fclose(H_file_temp);
			fclose(H_file);

			printf("sim_time %8.11f\n", sim_time);

		} /*end of if (time_step_ind >= output_freq && time_step_ind % output_freq == 0)*/
		
	} /*end of while (time_step_ind < num_of_frames) - OUTER LOOP*/
	/************************************************** END OF SIMULATION *****************************************/
	/*OK 1-7-2023.*/

	/************************************************* RELEASE MEMORY *********************************************/

	free_dvector(Xa, 1, Na_max);
	free_dvector(Ya, 1, Na_max); 
	
	free_dvector(leftover_dist_array_a, 1, Na_max);
	free_dvector(leftover_time_array_a, 1, Na_max);
	free_dvector(leftover_dir_array_a, 1, Na_max);

	free_dvector(BT_times_a, 1, Na_max);

	free_imatrix(cum_cell_vis_a, 1, numcellX, 1, numcellY);
	free_imatrix(curr_cell_partnum_a, 1, numcellX, 1, numcellY);

	free_d3tensor(id_array_a, 1, numcellX, 1, numcellY, 1, max_parts_in_cell);

	//free_dvector(leftover_diff_dist_array_a, 1, Na_max);
    //free_dvector(leftover_diff_dir_array_a, 1, Na_max);

	/********************************************************************************/
	free_dvector(Xb, 1, Nb_max);
	free_dvector(Yb, 1, Nb_max);

	free_dvector(leftover_dist_array_b, 1, Nb_max);
	free_dvector(leftover_time_array_b, 1, Nb_max);
	free_dvector(leftover_dir_array_b, 1, Nb_max);

	free_dvector(BT_times_b, 1, Nb_max);

	free_imatrix(cum_cell_vis_b, 1, numcellX, 1, numcellY);
	free_imatrix(curr_cell_partnum_b, 1, numcellX, 1, numcellY);

	free_d3tensor(id_array_b, 1, numcellX, 1, numcellY, 1, max_parts_in_cell);
	
	//free_dvector(leftover_diff_dist_array_b, 1, Na_max);
	//free_dvector(leftover_diff_dir_array_b, 1, Na_max);

	/********************************************************************************/
	
	free_ivector(inlet_part, 1, numcellY);
	free_dvector(qx, 1, numcellY);

	free_imatrix(cum_cell_diss_react, 1, numcellX, 1, numcellY);
	free_imatrix(cum_cell_prec_react, 1, numcellX, 1, numcellY);

	free_dmatrix(K_field, 1, numcellX, 1, numcellY);
	free_dmatrix(porosity_field, 1, numcellX, 1, numcellY);
	free_dmatrix(H_field, 1, numcellX + 1, 1, numcellY + 1);

	/*********************** END OF MEMORY RELEASE *******************************/
}