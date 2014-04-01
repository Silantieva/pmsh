#include "defns5.h"

void outgeomview(
int        mypid,
int        numprocs,
nglib::Ng_Mesh  *submesh,
char     *file)
{
   FILE   *fp ; 
   char filename[100] ; 
   char file2name[100] ;
   vector<float>  rgbvals ; 
   float rgb[3] ; 
   
   if(mypid == 0) {
      sprintf(filename,"%s.off",file) ; 
      fp = fopen(filename,"w") ; 
      if (fp == NULL) {
        printf("Cannot open file: %s\n",filename) ; 
        return ; 
      }
 
      fprintf(fp,"LIST\n") ;
      for(int p=0 ; p < numprocs ; p++) {
        sprintf(file2name,"%s.%d.off",file,p) ;
        fprintf(fp,"{ < %s }\n",file2name) ; 
      }
      fclose(fp) ; 
   }
   sprintf(file2name,"%s.%d.off",file,mypid) ;
   ncolors(numprocs,rgbvals) ;
   rgb[0] = rgbvals[3*mypid] ;  
   rgb[1] = rgbvals[3*mypid+1] ;
   rgb[2] = rgbvals[3*mypid+2] ;
   outgeomview(submesh,rgb,file2name) ; 

}

void outgeomview(
nglib::Ng_Mesh  *mesh,
float    *rgb,
char     *filename)
{
   FILE   *fp ; 
   double  xyz[3] ;
   int     verts[4] ; 
   int     np = Ng_GetNP(mesh) ;
   int     nse = Ng_GetNSE(mesh) ;  
    
   fp = fopen(filename,"w") ; 
   if (fp == NULL) {
     printf("Cannot open file: %s\n",filename) ; 
     return ; 
   }
   fprintf(fp,"OFF\n%d %d 0\n",np,nse) ; 
   
   for(int i=1 ; i <= np ; i++) {
     Ng_GetPoint (mesh,i,xyz) ;
     fprintf(fp,"%lf %lf %lf\n",xyz[0],xyz[1],xyz[2]) ; 
   }
   
   for(int i=1 ; i <= nse ; i++) {
       Ng_GetSurfaceElement(mesh,i,verts);
       fprintf(fp,"3 %d %d %d %f %f %f\n",verts[0]-1,verts[1]-1,verts[2]-1,
                                           rgb[0],rgb[1],rgb[2]) ; 
   }
   fclose(fp) ; 
}

void outsmy(
nglib::Ng_Mesh  *mesh,
char     *filename)
{
   FILE   *fp ; 
   double  xyz[3] ;
   int     verts[4] ; 
   int     np = Ng_GetNP(mesh) ;
   int     nse = Ng_GetNSE(mesh) ;  
    
   fp = fopen(filename,"w") ; 
   if (fp == NULL) {
     printf("Cannot open file: %s\n",filename) ; 
     return ; 
   }
   fprintf(fp,"%d %d\n",np,nse) ; 
   
   for(int i=1 ; i <= np ; i++) {
     Ng_GetPoint (mesh,i,xyz) ;
     fprintf(fp,"%lf %lf %lf\n",xyz[0],xyz[1],xyz[2]) ; 
   }
   
   for(int i=1 ; i <= nse ; i++) {
       Ng_GetSurfaceElement(mesh,i,verts);
       fprintf(fp,"%d %d %d\n",verts[0],verts[1],verts[2]) ; 
   }
   fclose(fp) ; 
}

void ncolors(
int n,
vector<float> & rgbvals)
{
   int i ; 
   float hsv[3],rgb[3]  ; 
   float adjbright, incrhue ;
   float brightbase = 0.6 ;
   float brightvar = 0.3 ; 
   float satur = 0.9 ; 
   
   adjbright = brightvar ; 
   hsv[0] = 0.0 ;
   incrhue = 1.0/n ; 
   hsv[1] = satur ; 
   
   for(i=0 ; i < n ; i++) {
      hsv[2] = brightbase + adjbright ; 
      hsv2rgb(rgb,hsv) ; 
      rgbvals.push_back(rgb[0]) ; 
      rgbvals.push_back(rgb[1]) ; 
      rgbvals.push_back(rgb[2]) ; 
      hsv[0] += incrhue ; 
      adjbright *= -1 ; 
   }
}

void hsv2rgb(
float *rgb, 
float *ihsv) 
{
	int i;
	float f, p, q, t;
	float hsv[3] ; 
	
	hsv[0] = ihsv[0]*360.0 ; hsv[1] = ihsv[1] ; hsv[2] = ihsv[2] ;
	if( hsv[1] < 1.0e-6 ) {
		rgb[0] = rgb[1] = rgb[2] = hsv[2] ;
		return;
	}
	
	hsv[0] /= 60.0 ;			
	i = floor( hsv[0] );
	f = hsv[0] - i;			
	p = hsv[2] * ( 1 - hsv[1] );
	q = hsv[2] * ( 1 - hsv[1] * f );
	t = hsv[2] * ( 1 - hsv[1] * ( 1 - f ) );
	switch( i ) {
		case 0:    rgb[0] = hsv[2]; rgb[1] = t;       rgb[2] = p;      break; 		
		case 1:    rgb[0] = q;      rgb[1] = hsv[2];  rgb[2] = p;      break; 
		case 2:    rgb[0] = p;      rgb[1] = hsv[2];  rgb[2] = t;      break;
		case 3:    rgb[0] = p;      rgb[1] = q;       rgb[2] = hsv[2]; break; 
		case 4:    rgb[0] = t;      rgb[1] = p;       rgb[2] = hsv[2]; break;
	    default:   rgb[0] = hsv[2]; rgb[1] = p;       rgb[2] = q;      break;			
	}
}
