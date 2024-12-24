#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int i, time_frame;
 
int main(void)
{
 void NRPT(void);

 printf("Where to start? ");
 scanf("%d", &time_frame);

 for (i= time_frame;i<1001;i++)
 {
 	 time_frame = i;
 	 printf("Doing work. %i\n",time_frame);
	 NRPT();		
 }

return 0;
}