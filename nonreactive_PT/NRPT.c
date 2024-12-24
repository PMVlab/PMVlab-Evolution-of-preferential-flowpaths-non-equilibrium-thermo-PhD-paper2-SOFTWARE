#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nr.h"
#include <math.h>
#include <memory.h>
#include <string.h>
#include <time.h>

/*====================================================================*/
/*                                                                    */
/*          NON-REACTIVE 2D PARTICLE TRACKER                          */
/*                                                                    */  
/*       LAST UPDATE 21-2-2024 BY EVGENY SHAVELZON                    */  
/*                                                                    */  
/*====================================================================*/

extern int time_frame; /*current frame number*/

void NRPT(void)
{
/*====================================================================*/
/*                     Declarations                                   */
/*====================================================================*/
	
    /*Numerical*/
	int Nx, Ny, Nsteps, Npart_TOT; /*cell numbers in longitudinal / transverse direction, how many steps to divide dx by to get ds*/
	int mid_point, particle_track; /*0 - distribute particles along inlet, else - inject all in the mid inlet Y location; 0 - do not print out particle trajectories, 1 - print around a selected value, 2 - print the whole path!; number of injected particles*/
	int i, j, k, x_cellnum, y_cellnum; /*array indices, x/y axis cell number*/
	int part_count = 0, monitor_flag, Npart_REAL = 0; /*auxilliary variable to distribute flux-weighted particles in the inlet cells, BTC monitor flag, exact number of injected particles*/

	long idum; // (-10709/-9456); /*Random seed, has to be a negative integer!*/
	
	float temp_K, temp_H; /*dummy vars for scanf*/

	double dx=0.2, dy=0.2, delta_s_const; /*cell dims, constant spatial step size [cm]*/
	double theta0, Ddiff; /*mean field porosity, should be a field itself (!!!), tracer diffusion coefficient [cm2/min]*/
	double BTC_watch_loc, time_watch; /*X coordinate for the BTC monitor location; time watch monitor value to print particle trajectories, [min]*/
	double A, B, C, Vx, Vy, V, phi, x_loc, y_loc, delta_t, tot_time, Qtot = 0.0; /*hyd. head gradients (inlcuded mixed deriv.???), advection velocity components, mag. & direction angle, current particle coords, dt to pass const. ds, total particle time, total inlet vol. flow rate*/
	double tot_delta_s, temp_Dx, temp_Dy; /*total particle path, diffusive contribution*/
	
	int * inlet_part; /*number of particles injected in an inlet cell array*/
	double * qcell, * y_loc_0, * BTC; /*inlet Darcy flux array, initial particle Y location array, BT times for all paricles*/
	
	int ** flag_cell_vis, ** cum_cell_vis; /*matrix that indicates which cells were visited by particle, total number of particle visitations in cell*/
	double ** K_field, ** H_field,  ** cum_cell_time, ** cum_cell_dist; /*Hyd. Cond. & Head field matrices, total time spent in cell, total distance passed in cell by visiting particles */

	/*Chars*/
	char cum_cell_vis_filename[30];
	char cum_cell_time_filename[30];
	char cum_cell_dist_filename[30]; 
	char BTC_filename[30];
	char K_filename[30];
	char H_filename[30];
	char particle_tracker_time_filename[50];  /*particle tracker filename*/
	char particle_tracker_x_loc_filename[50];  /*particle tracker filename*/
	char particle_tracker_y_loc_filename[50];  /*particle tracker filename*/

	/*Files*/
	FILE *cum_cell_vis_file; /*total num of particles in cell*/
	FILE *cum_cell_time_file; /*total time in cell*/
	FILE *cum_cell_dist_file; /*total dist in cell*/
	FILE *BTC_file; /*BTC times*/
    FILE *K_file;
	FILE *H_file;
	FILE* particle_tracker_time_file; /*particle tracker file*/
	FILE* particle_tracker_x_loc_file; /*particle tracker file*/
	FILE* particle_tracker_y_loc_file; /*particle tracker file*/
	
	/******************************* INPUT PARAMS - READ FROM FILE ***********************************/
	char Input_filename[50];
	FILE* Input_file;
	sprintf(Input_filename, "NRPT_input_file.txt");

	Input_file = fopen(Input_filename, "r");
	(void)fscanf(Input_file, "%*s %*s %*s %*s\n");
	(void)fscanf(Input_file, "%ld %d %lf %lf\n", &idum, &Npart_TOT, &theta0, &Ddiff);
	(void)fscanf(Input_file, "%*s %*s %*s %*s %*s %*s %*s\n");
	(void)fscanf(Input_file, "%d %d %d %lf %lf %d %d\n", &Nx, &Ny, &Nsteps, &BTC_watch_loc, &time_watch, &mid_point, &particle_track);
	fclose(Input_file);
	/*********************************************************************************/
	
	delta_s_const = dx / Nsteps;

	/*Arrays*/
	inlet_part    = ivector(1, Ny);
	qcell         = dvector(1, Ny);
	BTC           = dvector(1, Npart_TOT);
	y_loc_0       = dvector(1, Npart_TOT);

	cum_cell_vis  = imatrix(1, Nx,     1, Ny);
	flag_cell_vis = imatrix(1, Nx,     1, Ny);
	cum_cell_time = dmatrix(1, Nx,     1, Ny);
	cum_cell_dist = dmatrix(1, Nx,     1, Ny);
	K_field       = dmatrix(1, Nx,     1, Ny);
	H_field       = dmatrix(1, Nx + 1, 1, Ny + 1);

	/*Set file names for the current frame*/				 
	sprintf(cum_cell_vis_filename, "%d cum_cell_vis_file.dat", time_frame); // particle density matrix
	sprintf(cum_cell_time_filename, "%d cum_cell_time_file.dat", time_frame); // time spent by all particles in each cell
	sprintf(cum_cell_dist_filename, "%d cum_cell_dist_file.dat", time_frame); // distance passed by all particles in each cell
	sprintf(BTC_filename, "%d BTC.dat", time_frame); // particles BT times
	sprintf(particle_tracker_time_filename, "%d particle_tracker_time_file.dat", time_frame);
	sprintf(particle_tracker_x_loc_filename, "%d particle_tracker_x_loc_file.dat", time_frame);
	sprintf(particle_tracker_y_loc_filename, "%d particle_tracker_y_loc_file.dat", time_frame);

	/*Create output files*/
	cum_cell_vis_file = fopen(cum_cell_vis_filename, "w");
	cum_cell_time_file = fopen(cum_cell_time_filename, "w");
	cum_cell_dist_file = fopen(cum_cell_dist_filename, "w");
	BTC_file = fopen(BTC_filename,"w");
	
	if (particle_track > 0)
	{
		particle_tracker_time_file = fopen(particle_tracker_time_filename, "w");
		particle_tracker_x_loc_file = fopen(particle_tracker_x_loc_filename, "w");
		particle_tracker_y_loc_file = fopen(particle_tracker_y_loc_filename, "w");
	}

	/*Set input filenames*/
    if (time_frame == 0)
	{
	    sprintf(K_filename, "OrigTRS0001.txt");
	    sprintf(H_filename, "OrigFIP0001.txt");
	}
	else 
	{
		sprintf(K_filename, "TRS%04d.txt", time_frame);
		sprintf(H_filename, "FIP%04d.txt", time_frame);
	}
				
	/*Open input files*/
	K_file = fopen(K_filename,"r"); 
	H_file = fopen(H_filename,"r");

    if(K_file==NULL || H_file==NULL)
	{
	printf(K_filename);
	printf("Error: can't open file.\n");
    }
			
	/*Scan input files, initialize arrays.*/
	for (i=1; i<=Nx; i++)
	{
	     for (j=1; j<=Ny; j++)
	     {
			 fscanf(K_file, "%f", &temp_K);
			 K_field[i][j] = temp_K;
	
			 fscanf(H_file,"%f", &temp_H);
			 H_field[i][j] = temp_H;

			 cum_cell_time[i][j] = 0.0;
			 cum_cell_dist[i][j] = 0.0;
			 cum_cell_vis[i][j]  = 0;
		 }
		 fscanf(H_file, "%f", &temp_H);
		 H_field[i][j] = temp_H;
	} 

	for (j=1; j<=Ny+1; j++)
	{
		fscanf(H_file,"%f",&temp_H);
		H_field[i][j]=temp_H;
    } 

	fclose(K_file);
	fclose(H_file);

	/*Initialize*/
	for (k = 1; k <= Npart_TOT; k++) { BTC[k] = 0.0; }

	printf("Doing the job.\n");

	/*inlet vol. flow rate (~ to flux) per cell for flux-weighting*/
	for (j=1; j<=Ny; j++)
	{
		x_cellnum = 1;
		y_cellnum = j;
		A = ( H_field[x_cellnum+1][y_cellnum]-H_field[x_cellnum][y_cellnum] ) / dx; /*dh/dx*/
		qcell[j] = -exp(K_field[x_cellnum][y_cellnum]) * A * (dy * 1.0); /*inlet cell vol. flow rate [cm3/min/cm] = flux * area per unit in plane depth, ~ to flux.*/
		Qtot = Qtot + qcell[j]; /*total inlet vol. flow rate*/
    }
	
	/*distribute particles in inlet cells according to flux-weighting*/
	for (j = 1; j <= Ny; j++)
	{
	    inlet_part[j] = (int)floor(Npart_TOT*qcell[j]/Qtot); /*number of injected particles in each cell scales with relative vol. flow rate in that cell*/
		Npart_REAL = Npart_REAL + inlet_part[j];
		
		for (k = 1; k <= inlet_part[j]; k++)
		{
			part_count = part_count+1;
			y_loc_0[part_count] = ( (double)j-ran2(&idum) ) * dy; /*initial y location of particle within the current inlet cell, governed by uniform pdf*/
	    }
	} 
	
	/*****************************************************************/
	/*Run on each injected particle*/
	for (k = 1; k <= Npart_REAL; k++)
	{
		
		/*initialize vars for current particle*/
		tot_time = 0.0; /*total particle time in the field*/
		tot_delta_s = 0.0; /*total particle path in the field*/
		monitor_flag = 0; /*BTC monitor flag*/

		for (i = 1; i <= Nx; i++) { for (j = 1; j <= Ny; j++) { flag_cell_vis[i][j] = 0; } /*initialize first cell visitation flag*/ }
		
		/*init. X location*/
		x_loc = 0.0; 

		/*init. Y location*/
		if (mid_point == 0) {y_loc = y_loc_0[k];} /* flux-weight to distribute in all inlet cells - use init. Y location*/
		else {y_loc = (Ny/2)*dy;} /*inject all in the midcell*/ 
		
		fprintf(BTC_file, " %d ", k); /*BT times output file - make row for the current particle*/
		
		if (particle_track == 1)
		{
			fprintf(particle_tracker_time_file, " %d ", k); /*particle tracker time output - make file row for the current particle*/
			fprintf(particle_tracker_x_loc_file, " %d ", k); /*particle tracker x_loc output - make file row for the current particle*/
			fprintf(particle_tracker_y_loc_file, " %d ", k); /*particle tracker y_loc output - make file row for the current particle*/
		}
					
		while (x_loc < Nx*dx) // while particle hasn't reached outlet
		 {
 			  x_cellnum= (int)ceil(x_loc/dx); /*X, Y current cell indices*/
			  y_cellnum= (int)ceil(y_loc/dy);
		  
			  if (y_cellnum > Ny) {y_cellnum = Ny;} /*fix out of bounds*/
		      if (x_cellnum > Nx) {x_cellnum = Nx;}
		  
			  if (y_cellnum < 1){y_cellnum = 1;}
		      if (x_cellnum < 1){x_cellnum = 1;}
		  
		      A = (H_field[x_cellnum+1][y_cellnum] - H_field[x_cellnum][y_cellnum])/dx; /*dh/dx*/
		      B = (H_field[x_cellnum][y_cellnum+1] - H_field[x_cellnum][y_cellnum])/dy; /*dy/dx*/
		      C = (H_field[x_cellnum+1][y_cellnum+1] - H_field[x_cellnum][y_cellnum+1] - H_field[x_cellnum + 1][y_cellnum] + H_field[x_cellnum][y_cellnum])/(dx*dy);/*d2h/dxdy*/
			
             /*calculate current particle velocity*/
		     Vx = -exp(K_field[x_cellnum][y_cellnum]) / theta0 * ( A+C*( y_loc-((double)y_cellnum-1.0)*dy ) ); /*mixed derivative used to interpolate within cell, OK*/
		     Vy = -exp(K_field[x_cellnum][y_cellnum]) / theta0 * ( B+C*( x_loc-((double)x_cellnum-1.0)*dx ) );
		  
			 V = sqrt(Vx*Vx+Vy*Vy);  
		     phi = atan(Vy/Vx); /*flow direction angle wrt horizontal axis, [-pi/2, +pi/2]*/
		  
			 if (fabs(phi) < 1.0e-9) { phi = 0.0; }
			 if (fabs(phi) > 1.57079) { phi =(phi/ fabs(phi))*1.57079; }

		     delta_t = delta_s_const / V; /*current jump time step*/
		     tot_time = tot_time + delta_t; /*total particle time (neglect diffusive contribution)*/

			 temp_Dx = sqrt(2.0 * Ddiff * delta_t) * (gasdev(&idum));
			 x_loc = x_loc + delta_s_const * cos(phi) + temp_Dx; /*current location*/
				
			 temp_Dy = sqrt(2.0 * Ddiff * delta_t) * (gasdev(&idum));
			 y_loc = y_loc + delta_s_const * sin(phi) + temp_Dy;
			
			 tot_delta_s = tot_delta_s + sqrt( pow(delta_s_const*cos(phi) + temp_Dx, 2) + pow(delta_s_const*sin(phi) + temp_Dy, 2) );

		     if (y_loc > Ny*dy) {y_loc = Ny*dy + (Ny*dy-y_loc);} /*out of bounds, reflective BC*/
 		     if (y_loc < 0.0 ) { y_loc = -y_loc; }
		     if (x_loc < 0.0 ) { x_loc = -x_loc; }  
		  
   		     if (flag_cell_vis[x_cellnum][y_cellnum] == 0) /*if the current cell is being visited by the particle for the first time*/
		     {
				 flag_cell_vis[x_cellnum][y_cellnum] = 1;
				 cum_cell_vis[x_cellnum][y_cellnum]  = cum_cell_vis[x_cellnum][y_cellnum] + 1; /*assuming particle can visit cell only once*/

				 cum_cell_time[x_cellnum][y_cellnum] = cum_cell_time[x_cellnum][y_cellnum] + delta_t; /*assuming particle flight always begins and ends in the same cell???*/
			     cum_cell_dist[x_cellnum][y_cellnum] = cum_cell_dist[x_cellnum][y_cellnum] + sqrt( pow(delta_s_const*cos(phi) + temp_Dx, 2) + pow(delta_s_const*sin(phi) + temp_Dy, 2) );
  		     }  
		     else
		     {
			     cum_cell_time[x_cellnum][y_cellnum] = cum_cell_time[x_cellnum][y_cellnum] + delta_t;
			     cum_cell_dist[x_cellnum][y_cellnum] = cum_cell_dist[x_cellnum][y_cellnum] + sqrt( pow(delta_s_const*cos(phi) + temp_Dx, 2) + pow(delta_s_const*sin(phi) + temp_Dy, 2) );
	  	     }

			 if (x_loc > BTC_watch_loc && monitor_flag == 0) /*if the particle has just passed the BTC monitor*/
			 {
				 monitor_flag = 1;
				 BTC[k] = tot_time - (x_loc - BTC_watch_loc) / (V*cos(phi)); /*calc and print to file BT time past the BTC monitor*/
				 fprintf(BTC_file, " %15.8f ", BTC[k]);
			 }

			 if ( (particle_track == 2) || (particle_track == 1 && tot_time >= 0.9 * time_watch && tot_time <= 1.1 * time_watch) ) /*if the particle time has just passed the Time watch monitor*/
			 {
				 fprintf(particle_tracker_time_file, " %15.8f ", tot_time);
				 fprintf(particle_tracker_x_loc_file, " %15.8f ", x_loc);
				 fprintf(particle_tracker_y_loc_file, " %15.8f ", y_loc);
			 }

		 } /*end of while particle has not reached the outlet*/
				   
		 /*BT time for the outlet*/
	 	 BTC[k] = tot_time - (x_loc - Nx*dx) / (V*cos(phi));

		 fprintf(BTC_file," %15.8f %15.8f\n",BTC[k], tot_delta_s);
		 /*BTC file: num. part, BTC time monitor, BTC time outlet, total part. path*/

		 if (particle_track == 1)
		 {
			 fprintf(particle_tracker_time_file, "\n");
			 fprintf(particle_tracker_x_loc_file, "\n");
			 fprintf(particle_tracker_y_loc_file, "\n");
		 }

    } /*end of particle number loop*/
	/*************************************************************/
	
	/*print output*/
	for (i = 1; i <= Nx; i++)
	{
		for (j = 1; j <= Ny; j++) 
		{
	        fprintf(cum_cell_time_file, " %5.8f ", cum_cell_time[i][j]);
			fprintf(cum_cell_dist_file, " %5.8f ", cum_cell_dist[i][j]);
			fprintf(cum_cell_vis_file,  " %d ",     cum_cell_vis[i][j]);
		}
		fprintf(cum_cell_time_file, " \n");
		fprintf(cum_cell_dist_file, " \n");
		fprintf(cum_cell_vis_file," \n");
	}

 	/*close output files*/
	fclose(cum_cell_vis_file);
	fclose(cum_cell_time_file);
	fclose(cum_cell_dist_file);
    fclose(BTC_file);
	
	if (particle_track == 1)
	{
		fclose(particle_tracker_time_file);
		fclose(particle_tracker_x_loc_file);
		fclose(particle_tracker_y_loc_file);
	}
}
