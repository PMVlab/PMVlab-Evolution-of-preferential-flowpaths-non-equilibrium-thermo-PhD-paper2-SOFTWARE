#include <stdio.h>
#include <stdlib.h>
//#include "nrutil.h"
//#include "nr.h"
#include <math.h>
#include <memory.h>
#include <string.h>

extern double **K_field_mat, **Head_field_mat;
extern double bin_x, bin_y,Ux, Uy;
double q1, q2, q3, q4, q1ce, q2ce, q3ce, q4ce;
int ix, iy, Add_X, Add_Y;
	double ddx,ddy,dz,ddz, Qx, Qy;
	double fac,rsq,v1,v2,csi,ic,ir,eta,detj;
	double q1f,q2f,q3f,q4f,q1,q2,q3,q4,q1int,q2int,q3int,q4int;
	double q1ce,q2ce,q3ce,q4ce,q1n,q2n,q3n,q4n,q1p,q2p,q3p,q4p;

void Kinz_corr(int loc_x,int loc_y,double spatial_X,double spatial_Y)
{
   
	
	

	void faiQiEl(int ix,int iy);
	void faiQce(double q1,double q2,double q3,double q4);
	ix=loc_x;
	iy=loc_y;

		/*
	  identify sub-element containing spatial_X 
	  define distance between spatial_X and the left hand and bottom sides of the element ---DDX=spatial_X-bin_x, ddy=spatial_Y-bin_y
	  */
		ddx=spatial_X-(ix)*bin_x;
		ddy=spatial_Y-(iy)*bin_y;
		/*ddx=spatial_X-(ix-1)*bin_x;
		ddy=spatial_Y-(iy-1)*bin_y;*/
		dz=bin_x/2;
		ddz=bin_x/4;

		
		faiQiEl(ix,iy);
		q1p=q1;
		q2p=q2;
		q3p=q3;
		q4p=q4;
	/*
	 Sub 'faiQce' calculates the fluxes 
     along the sides (horiz and vertical) connecting the mean points of the element where spatial_X is,
    i.e., along the internal sides of the sub-elements constructed in the current element
	 (eq.10-12 paper)
	 */
	 faiQce(q1p,q2p,q3p,q4p);
//	   saves fluxes on internal sides
	q1int= q1ce;
	q2int= q2ce;
	q3int= q3ce;
	q4int=q4ce;
	
	
/*      
	  Define the fluxes (qif) along all sides of the current subelement
       interpolation is performed according to Eq.(17-18)
       rendering velocity in the given point

	WE CAN NOW HAVE FOUR SITUATIONS, AS DETAILED IN WHAT FOLLOWS
	 1- IF spatial_X belongs to the top-right sub-element 
	 */
	if (ddx>=dz && ddy>=dz )
	{
	q1f=-q2int;
	q4f=q3int;
/*
      to define the fluxes along the the external sides of the subelements containing spatial_X (q2f and q3f)
      - we need to consider the top-right node of the element  
      - we need to define (employing sub 'faiQiEl') the contributions to it
        from the elements which communicate with the right top node 
      Once these fluxes are known, THEN calculate q2f and q3f through 'faiQce'
	  
 	  q1n is the contribution to the node coming from the element containing spatial_X
	  */
	q1n=-q3p;
	/*
	 look up for the other contributions to the node by interrogating all other elements 
      communicating with the node
	  */
	faiQiEl(ix,iy+1);
	q4n=-q2;
	faiQiEl(ix+1,iy);
	q2n=-q4;
	faiQiEl(ix+1,iy+1);
	q3n=-q1;

//	 calculate the fluxes q2f and q3f along the sides of sub-element containing spatial_X
      faiQce(q1n,q2n,q3n,q4n);
	  q2f=q1ce;
	q3f=-q4ce;
	}

//	 2- IF spatial_X belongs to the bottom-right sub-element 

	if (ddx >= dz && ddy<dz )
	{
	q3f=q2int;
	q4f=-q1int; 

//	proceed with calculating fluxes along external sides of the sub-element
      q4n=-q2p;

	faiQiEl(ix+1,iy);
	q3n=-q1;

	faiQiEl(ix,iy-1);
	q1n=-q3;

	faiQiEl(ix+1,iy-1);
	q2n=-q4;

	 faiQce(q1n,q2n,q3n,q4n);
	q1f=q4ce;
	q2f= -q3ce;
	}

//	 3- IF spatial_X belongs to the bottom-left sub-element 

	if (ddx<dz && ddy<dz )
	{

	
//	 define fluxes along the internal sides
	 
	q3f=-q4int;
 	q2f=q1int;
//	proceed with calculating fluxes along external sides of the sub-element

	q3n=-q1p;
	faiQiEl(ix-1,iy);
	q4n=-q2;

	faiQiEl(ix,iy-1);
	q2n=-q4;

	faiQiEl(ix-1,iy-1);
	 q1n=-q3;
	faiQce(q1n,q2n,q3n,q4n);
	q1f=-q2ce;
	q4f= q3ce;
	
	}

//4- IF spatial_X belongs to the top-left sub-element 

	if (ddx<dz && ddy>=dz ) 
	{
	
//	 define fluxes along the internal sides
		
  	q1f=q4int;
 	q2f=-q3int;
//	proceed with calculating fluxes along external sides of the sub-element

	q2n=-q4p;

	
	faiQiEl(ix,iy+1);
	q3n=-q1;
	faiQiEl(ix-1,iy);
	q1n=-q3;
	faiQiEl(ix-1,iy+1);
	q4n=-q2;
	faiQce(q1n,q2n,q3n,q4n);
	q4f=-q1ce;
	q3f= q2ce;
	
	}
/*	
	Once we know the fluxes along the sides of the sub-element, we can calculate the velocity
	at a desired point upon performing linear interpolation of the fluxes calculated
	calculated along the sides of the sub-elements
	 Note that (Eq.(17), (18) are slightly different, but they should be equivalent 
	to the formulation here adopted
      
	 ic and ir identify the numbering of the half-column / half-rows
	 */
	 ic=floor(spatial_X/dz);
	 ir=floor(spatial_Y/dz);
	 
//	csi and eta vary between -1 and 1
	csi= ( (spatial_X-(ic*dz))-ddz ) /ddz;
      eta= (spatial_Y-(ir*dz) -ddz)	 /ddz;
	  /*
c 	 the denominator includes the area of passage of the flux qif, i.e., (Dx/2), 
c      to go from fluxes to velocity
*/
	detj=2*(bin_x/2);

	Qx= ( (q2f-q4f) + csi*(q2f+q4f))/detj;
	Qy= ((q3f-q1f)+eta*(q3f+q1f))/detj;
	Ux=Qx;
	Uy=Qy;  

	
	
	
	
	
	
	
	/*
	
	
	faiQiEl(loc_x, loc_y);
		  if (loc_x>0 && loc_y>0)
		  {
			  Add_X=loc_x+floor(ddx/dz)-1;
			  Add_Y=loc_y+floor(ddy/dz)-1;
			  for (iy=0;iy<=1;iy++)
			  {
				  for (ix=0;ix<=1;ix++)
				  {
					  Add_Y=Add_Y+iy;
					  Add_X=Add_X+ix;
					 if (Add_Y!=loc_y
					faiQiEl( ix, iy);
				    
					 
				  }
			}
		  faiQce( q1, q2, q3, q4);

		}


		*/
	}

