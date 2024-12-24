#include <stdlib.h>
#include "nrutil.h"
#include "nr.h"
#include <math.h>
#include <memory.h>
#include <string.h>
//#include <stdio.h>

/*external vars*/
extern long idum; /*random seed*/
extern int numcellY, numcellX; /*number of cells in X, Y directions*/
extern double dT, delta_s_const, D, dx, dy, x_loc, y_loc; /*simulation time increment, constant particle spatial step, diffusion coefficient [cm2/min], cell dimensions in X, Y, particle X,Y coordinates*/
extern double leftover_dist, leftover_time, leftover_dir; /*"redundant" parameters from the previous time step: leftover distance, time and advective velocity direction*/
extern double **H_field, **K_field, **porosity_field; /*H, K, porosity fields*/
//extern double leftover_dist_diff, leftover_dir_diff; /*CURRENTLY UNUSED "redundant" parameters from the previous time step: leftover diffusive distance and diffusive movement direction*/

/*====================================================================  */
/*   Routine for advancing particles in the domain using 2D Langevin eq.*/
/*          (GASDEV.C)                                                  */
/*       LAST MODIFIED BY EVGENY SHAVELZON 16-2-2024                    */
/*======================================================================*/

void move(void)
{
	int x_cell, y_cell; /*particle location cell numbers*/
	double A, B, C, Vx, Vy, V; /*head gradients, flow velocities*/
	double t = 0;  /*local jump time*/
	double dist_adv, dist_diff, phi_adv; /*magnitude of advective/diffusive jump (partial), direction of advective jump*/
	double delta_s, delta_t;  /*full leftover advective jump length, leftover distance from previous time step, jump duration*/

	/*FIRST expend the remainder of the previous jump, then proceed to the current if time is left*/
	delta_s = leftover_dist; /*assign "redundant" distance, to be expended first, at the beginning of the jump to a local var*/
	phi_adv = leftover_dir; /*assign "redundant" advective velocity direction*/
	delta_t = leftover_time; /*assign "redundant" time after dT*/

	while (t < dT) /*advance particle until the end of the current time step is reached*/
	{
		if (delta_s == 0) /*if no jump remainder left from the previous time step*/
		{
			/*current cell indices, ceil number rounds to the nearest integer from above. Cell numbers go from 1:Nx, 1:Ny, thus need to fix for x_loc=0  (cannot be since all parts are given initial nonzero X coord within the 1st cell)*/
			x_cell = (int)ceil(x_loc / dx); 
			y_cell = (int)ceil(y_loc / dy);

			/*check bounds just in case*/
			if (y_cell > numcellY) { y_cell = numcellY; } /*out of bounds*/
			if (x_cell > numcellX) { x_cell = numcellX; }
			if (y_cell < 1) { y_cell = 1; } /*Cell numbers go from 1:Nx, 1:Ny, thus need to fix for x_loc=0 (cannot be since all parts are given initial nonzero X coord witihn the 1st cell)*/
			if (x_cell < 1) { x_cell = 1; }

			/*head gradient components. MIXED DERIVATIVE TO INTERPOLATE WITHIN CELL*/
			A = ( H_field[x_cell + 1][y_cell] - H_field[x_cell][y_cell] ) / dx; 
			B = ( H_field[x_cell][y_cell + 1] - H_field[x_cell][y_cell] ) / dy;
			C = ( H_field[x_cell + 1][y_cell + 1] - H_field[x_cell][y_cell + 1] - H_field[x_cell + 1][y_cell] + H_field[x_cell][y_cell] ) / (dx*dy);
			/* check for positive C???*/

			/*get flow velocity components [cm/min] in the current cell from Darcy*/
			Vx = -exp( K_field[x_cell][y_cell] ) / porosity_field[x_cell][y_cell] * ( A + C * (y_loc - ((double)y_cell-1.0)*dy) ); // PAY ATTENTION THIS READ BEFORE EXP(K/THETA)! 
			Vy = -exp( K_field[x_cell][y_cell] ) / porosity_field[x_cell][y_cell] * ( B + C * (x_loc - ((double)x_cell-1.0)*dx) );
			V = sqrt(Vx*Vx + Vy*Vy); 
			phi_adv = atan(Vy / Vx);
			
			/*duration of a jump - equals to a time for a particle at the current location to travel delta_s_const*/
			delta_t = delta_s_const / V; 

			if (t + delta_t < dT) /*if the whole delta_t jump fits into dT interval - full flight duration!*/
			{
				t = t + delta_t; /*update current time*/
				dist_adv = delta_s_const; /*advective contribution*/
				
				dist_diff = sqrt(2.0 * D * delta_t) * gasdev(&idum);  /* diffusive contribution, Gaussian dist*/
				x_loc = x_loc + dist_adv * cos(phi_adv) + dist_diff; /*current location*/
				
				dist_diff = sqrt(2.0 * D * delta_t) * gasdev(&idum);  /* diffusive contribution, Gaussian dist*/
				y_loc = y_loc + dist_adv * sin(phi_adv) + dist_diff;

				/*if out of bounds, apply reflective BC. Exit field checked in the main program*/
				if (y_loc > (double)numcellY * dy) { y_loc = (double)numcellY * dy + ((double)numcellY * dy - y_loc); } 
				if (y_loc < 0.0) { y_loc = -y_loc; }
				if (x_loc < 0.0) { x_loc = -x_loc; }  

				delta_s = 0.0; /*just to make order*/

			} /*end of if (t + delta_t < dT)*/

			else /*if the whole delta_t jump doesn't fit into dT interval - partial flight duration (leave remainder for the next time step)*/
			{
				t = t + delta_t; /*update current time by delta_t*/
				dist_adv = ( (delta_t - (t - dT)) / delta_t ) * delta_s_const; /*relative displacement part to reach dT. (delta_t - (t - dT)) / delta_t is just the relative jump part to reach dT.*/

				dist_diff = sqrt(2.0 * D * (delta_t - (t - dT))) * gasdev(&idum); /*diffusion contribution, relative part to reach dT*/
	 		    x_loc = x_loc + dist_adv * cos(phi_adv) + dist_diff; /*update particle location*/
				
				dist_diff = sqrt(2.0 * D * (delta_t - (t - dT))) * gasdev(&idum); /*diffusion contribution, relative part to reach dT*/
				y_loc = y_loc + dist_adv * sin(phi_adv) + dist_diff;

				/*if out of bounds, apply reflective BC. Exit field checked in the main program*/
				if (y_loc > (double)numcellY * dy) { y_loc = (double)numcellY * dy + ((double)numcellY * dy - y_loc); }
				if (y_loc < 0.0) { y_loc = -y_loc; }
				if (x_loc < 0.0) { x_loc = -x_loc; }  

				delta_s = ( (t - dT) / delta_t ) * delta_s_const; /*"redundant" distance, the remainder to be passed during the next time step*/
				delta_t = t - dT; /*"redundant" time after dT, to be saved for the next time step*/

			}
		} /*end of if (delta_s == 0) - no jump remainder left from the previous time step. OK 19-6-2023.*/

		else /*if (delta_s != 0), if some jump remainder is left from the previous time step. relevant when the jump has continued into the next time interval -> the particle begins with delta_s>0*/
		{
			if (t + delta_t < dT) /*if the whole remainder delta_t jump fits into dT interval - full flight duration for the jump remainder delta_s!*/
			{
				t = t + delta_t; /*update jump time*/
				dist_adv = delta_s;  

				dist_diff = sqrt(2.0 * D * delta_t) * gasdev(&idum);
				x_loc = x_loc + dist_adv * cos(phi_adv) + dist_diff;

				dist_diff = sqrt(2.0 * D * delta_t) * gasdev(&idum);
				y_loc = y_loc + dist_adv * sin(phi_adv) + dist_diff;

				/*if out of bounds, apply reflective BC. Exit field checked in the main program*/
				if (y_loc > (double)numcellY * dy) { y_loc = (double)numcellY * dy + ((double)numcellY * dy - y_loc); }
				if (y_loc < 0.0) { y_loc = -y_loc; }
				if (x_loc < 0.0) { x_loc = -x_loc; }  
				 
				delta_s = 0.0; /*set jump remainder to zero, can proceed with the current time step jump now!*/
			}

			else /*if the whole remainder delta_t jump doesn't fit into dT interval - partial flight duration (leave remainder of the remainder for the next time step...)*/
			{   /*BETTER NOT TO ARRIVE AT THIS!*/
				t = t + delta_t; 
				dist_adv = ( (delta_t - (t - dT)) / delta_t ) * delta_s;

				dist_diff = sqrt(2.0 * D * (delta_t - (t - dT))) * gasdev(&idum);
				x_loc = x_loc + dist_adv * cos(phi_adv) + dist_diff;

				dist_diff = sqrt(2.0 * D * (delta_t - (t - dT))) * gasdev(&idum);
				y_loc = y_loc + dist_adv * sin(phi_adv) + dist_diff;

				/*if out of bounds, apply reflective BC. Exit field checked in the main program*/
				if (y_loc > (double)numcellY * dy) { y_loc = (double)numcellY * dy + ((double)numcellY * dy - y_loc); }
				if (y_loc < 0.0) { y_loc = -y_loc; }
				if (x_loc < 0.0) { x_loc = -x_loc; } 

				delta_s = ( (t - dT) / delta_t ) * delta_s; /*updating the remainder of the remainder...*/
				delta_t = t - dT;
			}

		} /*end of if (delta_s != 0) - some jump remainder is left from the previous time step.*/
	
	} /*while (t < dT)*/

	/*update "leftovers"*/
	leftover_dir = phi_adv; /*advective velocity direction*/
	leftover_dist = delta_s; /*"redundant" distance after t=dT*/
	leftover_time = delta_t; /*"redundant" time after dT*/
	//leftover_dist_diff = dist_diff;
	//leftover_dir_diff = phi_diff;
	/*OK 16-2-2024.*/
} 
