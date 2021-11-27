/**********************************************************************
***********************************************************************

   Version:        1.1
   Last modified:  July 25, 201
   Authors:        
     Olaf Schenk
  
***********************************************************************

     laplace.c is the driver routine that computes the selected
     inversion of a 3D Laplacian operator 
     -\partial^2/\partial x^2-\partial^2/\partial y^2-\partial ^2/\partial z^2
     with Dirichlet boundary condition, discretized on a 3D square domain 
     [0, 1] * [0, 1] * [0, 1] with grid size nx * ny * nz.

     
     Usage: laplace.x -nx=<nx> -ny=<ny> -nz=<nz> -chkerr=<0|1> -complx=<0|1> -cfull=<0|1> -cstore=<0|1> -printa=<0|1>\n");
     -nx        :   The number of grid points in x direction. (default: nx=10 )
     -ny        :   The number of grid points in y direction. (default: ny=10)
     -nz        :   The number of grid points in z direction. (default: nz=10)
     -nb        :   The number of dense rhs. (default: nb=10)

     -chkerr    :   Whether full inversion is to be computed for
		    efficiency comparison.  
		    chkerr = 0 (default). Only perform selected inversion.
		    chkerr = 1. Both selected inversion and full
    		      inversion are computed.

     -complx    :   Generate symmetric matrix that is  
                    complx = 0 (default). real and indefinite.
		    complx = 1. Complex symmetric

     -cstore    :   Storage for inverse factor
                    cstore = 0 (default) overwrite factor
		    cstore = 1  do not overwrite factor
		    cstore = 2  Store values on disc

     -cfull     :   Return symmetric entries in A=[ia, ja, a]
                    cfull = 0 upper symmetric CSR format (default) 
                    cfull = 1 full  symmetric CSR format

     -printa    :   Output the sparse matrix A. 
                    printa = 0 (default). Do not output A.
		    printa = 1 . Output A.
***********************************************************************
**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include <unistd.h>
#include <string.h>

#define mesh(i,j,k) mesh[nx*ny*((k)-1) + nx*((j)-1)+(i)-1]
#define nzvals(i) nzvals[(i)-1]
#define rowind(i) rowind[(i)-1]
#define colptr(i) colptr[(i)-1]
#define xsol(i)   xsol[(i)-1]
#define rhs(i)    rhs[(i)-1]
#define diag(i)   diag[(i)-1]
#define diag2(i)  diag2[(i)-1]


#ifdef TIMING
extern double getime(void);
#endif

int main(int argc, char ** argv)
{
   int i, j, k, nnodes,nedges, ia, nnz, ibeg, iend ;
   int nx= -1, ny = -1, nb = -1, nz = -1, complx = 0, cfull = 0, count, node;
   int *mesh;
   int *rowind, *colptr;
   double *nzvals, *nzvals2, *rhs, *xsol, *diag, *diag2;
   int token, order=-1;
   int Lnnz, nnz_full;
   double t0,t1, errmax; 
   long myclktck;
   double dval, errabs;
   int chkerr = 0, cdiag = 0, printa = 0, dumpL=0, cstore=0;
   int ierr = 0;
   double hx, hy, hz;  /* mesh size */


   ia = 1;
   while (ia < argc) {
      if ( !strncmp(argv[ia],"-nx",3) ) {
 	 nx = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-ny",3) ) {
 	 ny = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-nz",3) ) {
    	 nz = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-nb",3) ) {
    	 nb = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-chkerr",7) ) {
    	 chkerr = atoi(&argv[ia][8]);
         if (chkerr != 0) chkerr = 1;
      }
      else if ( !strncmp(argv[ia],"-complx",7) ) {
    	 complx = atoi(&argv[ia][8]);
         if (complx != 0) complx = 1;
      }
      else if ( !strncmp(argv[ia],"-cdiag",6) ) {
    	 cdiag = atoi(&argv[ia][7]);
         if (cdiag != 0) cdiag = 1;
      }
      else if ( !strncmp(argv[ia],"-cstore",7) ) {
    	 cstore = atoi(&argv[ia][8]);
         if (cstore > 2  || cstore < 0 ) cstore = 0;
      }
      else if ( !strncmp(argv[ia],"-cfull",6) ) {
    	 cfull = atoi(&argv[ia][7]);
         if (cfull != 0) cfull = 1;
      }
      else if ( !strncmp(argv[ia],"-printa",7) ) {
    	 printa = atoi(&argv[ia][8]);
         if (printa != 0) printa = 1;
      }
      else {
 	 fprintf(stderr, "invalide argument!\n");
	 fprintf(stderr, "Usage: laplace -nx=<nx> -ny=<ny> -nz=<nz> -nb=<nb> -chkerr=<0|1> -complx=<0|1> -cdiag=<0|1> -cstore=<0|1> -cfull=<0|1>  -printa=<0|1>\n");
         return 1;
      }
      ia++;
   }

   if (nx == -1 && ny > 0) nx = ny;
   if (ny == -1 && nx > 0) ny = nx;
   if (nz == -1 && nz > 0) nz = nx;

   if (nx == -1) nx = 10;
   if (ny == -1) ny = 10;
   if (nz == -1) nz = 10;

   hx = 1.0 / (nx+1);
   hy = 1.0 / (ny+1);
   hz = 1.0 / (nz+1);

   nnodes = nx*ny*nz;
   mesh = (int*)malloc(nnodes*sizeof(int));

   for (i = 0; i<nnodes; i++) mesh[i]=i+1;
   
   /* first pass to count the number of edges */
   /* Dirichlet BC */ 
   nedges = 0; 
  
   for (k = 1; k <= nz; k++) {
     for (j = 1; j <= ny; j++) {
       for (i = 1; i <= nx; i++) {
       if (k < nz) nedges++;
       if (k > 1)  nedges++;
       if (j < ny) nedges++;
       if (j > 1)  nedges++;
       if (i < nx) nedges++;
       if (i > 1)  nedges++;   
     }
   }
   }

   /* print the matrix dimension and number of nonzeros */
   nnz = nedges/2 + nnodes;
   printf("%d  %d\n", nnodes, nnz);
   nnz_full = nnz;
   if (cfull == 1) 
      nnz_full = 2*nnz;
   
   colptr = (int*)malloc((nnodes+1)*sizeof(int));
   rowind = (int*)malloc(nnz_full*sizeof(int));

   if (!complx)
   {
   	nzvals  = (double*)malloc(nnz_full*sizeof(double));
   	nzvals2 = (double*)malloc(nnz_full*sizeof(double));
	rhs     = (double*)calloc(nnodes,sizeof(double));
        xsol    = (double*)calloc(nnodes,sizeof(double));
   }
   else
   {
   	nzvals  = (double*)malloc(2*nnz_full*sizeof(double));
   	nzvals2 = (double*)malloc(2*nnz_full*sizeof(double));
	rhs     = (double*)calloc(2*nnodes,sizeof(double));
        xsol    = (double*)calloc(2*nnodes,sizeof(double));
   }

   colptr(1) = 1;
   count = 0;
   node  = 0;

   dval = 2.0/(hx*hx) + 2.0/(hy*hy) + 2.0/(hz*hz);

   dval = 20.0;

   printf(" Dirichlet boundary condition\n");

   for (k = 1; k <= nz; k++) {
     for (j = 1; j <= ny; j++) {
        for (i = 1; i <= nx; i++) {
	 /* diagonal */
         if (printa)
            printf("%d %d  %8.2e\n", mesh(i,j,k), mesh(i,j,k), dval); 

         rowind[count] = mesh(i,j,k);

	 if (!complx)
         {
		 nzvals[count] = (i + (j-1)*nx + (k-1)*nx*ny); /* no pivots */
		 nzvals[count] = nx*nx*ny; /* no pivots */
	 }
	 else
         {
		 nzvals[2*count  ] = (i + (j-1)*nx + (k-1)*nx*ny); /* no pivots */
		 nzvals[2*count+1] = -0.0; /* complex entry */

	 }        
         count++;

         /* lower */
         if (i < nx) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i+1,j,k), mesh(i,j,k));

            rowind[count] = mesh(i+1,j,k);
            dval = -1.0;
	    if (!complx)
            {
            	nzvals[count] = dval;
	    }
	    else
            {
            	nzvals[2*count  ] = dval;
            	nzvals[2*count+1] = -0.1*((i + (j-1)*nx + (k-1)*nx*ny));    /* complex entry */
	    }
            count++;
         }


         /* right */
         if (j < ny) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i,j+1,k), mesh(i,j,k)); 

            rowind[count] = mesh(i,j+1,k);
            dval = -1.0;
	    if (!complx)
            {
            	nzvals[count] = dval;
	    }
	    else
            {
            	nzvals[2*count  ] = dval;
            	nzvals[2*count+1] = -0.1*((i + (j-1)*nx + (k-1)*nx*ny));    /* complex entry */
	    }
            count++;
         }       

         /* lower */
         if (k < nz) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i,j,k+1), mesh(i,j,k));

            rowind[count] = mesh(i,j,k+1);
            dval = -1.0;
	    if (!complx)
            {
            	nzvals[count] = dval;
	    }
	    else
            {
            	nzvals[2*count  ] = dval;
            	nzvals[2*count+1] = -0.1*((i + (j-1)*nx + (k-1)*nx*ny));    /* complex entry */
	    }
            count++; 
         }

         node++;
         colptr(node+1) = count+1; 
         
      } 
   }
   }
   if (count != nnz) {
       printf(" is this correct? count = %d, nnz = %d\n", count, nnz);  
       return 1;
   }

   nnz = count;  
   for (i = 0; i < nnz; i++) 
   { 
	    if (!complx)
            {
            	nzvals2[i] = nzvals[i];
	    }
	    else
            {
            	nzvals2[2*i  ] = nzvals[2*i  ];
            	nzvals2[2*i+1] = nzvals[2*i+1];
	    }
   }


   token = 0;
   for (j = 1; j <= nnodes; j++)
	if (!complx)
	   rhs(j) = j;
	else
	{	
	   rhs(2*j-1) = j;
	   rhs(2*j)   = -0.000*j;
	}	

   if (printa) 
   { 
	printf("%d \n", nnodes );
	printf("%d \n", colptr(nnodes) );
        for (j = 1; j <= nnodes+1; j++) 
		printf("%d \n", colptr(j) );
        for (j = 1; j <= colptr(nnodes); j++) 
		printf("%d \n", rowind(j) );
        for (j = 1; j <= colptr(nnodes); j++)
        { 
		if (!complx)
			printf("%e \n", nzvals(j) );
		else
			printf("%e %e \n", nzvals(2*j-1), nzvals(2*j) );			
        }
   }

   /* add additional dense row of size nb */
 
   if (nb==-1) 
   {
       exit(9);
   }
   else
   {
      int NDENSE;
      int *rowindD, *colptrD;
      double *nzvalsD,*rhsD, *xsolD;
      int nnzD = 0; 
      int nnzDrow; 

      colptrD = (int*)malloc((nnodes+1+nb)*sizeof(int));
      rowindD = (int*)malloc((nnz_full+nb*nnodes +nb*nb/2+nb)*sizeof(int));
      nzvalsD  = (double*)malloc((nnz_full+nb*nnodes +nb*nb/2+nb)*sizeof(double));
      rhsD     = (double*)calloc(nnodes+nb,sizeof(double));
      xsolD    = (double*)calloc(nnodes+nb,sizeof(double));
      NDENSE = nnodes + nb; 


      colptrD[0] = 1;

      for (i = 0; i < nnodes; i++) {
       nnzDrow = 0; 
       for (j = colptr[i]-1; j <colptr[i+1]-1; j++) {
         rowindD[nnzD] = rowind[j]; 
         nzvalsD[nnzD]= nzvals[j];
         nnzD += 1;
         nnzDrow += 1; 
       }      
       for (k = 1; k <= nb; k++) {
         rowindD[nnzD] = nnodes + k; 
         nzvalsD[nnzD]=  -1.0/nb;
         nnzD +=1;
         nnzDrow += 1; 
       }
       colptrD[i+1] = colptrD[i] + nnzDrow;
     }

      for (i = 0; i < nb; i++) {
       nnzDrow = 0; 
       for (j = i; j < nb; j++) {
         rowindD[nnzD] = nnodes + j + 1; 
         if (i == j ) 
           nzvalsD[nnzD]=  1.0*nb;
         else
           nzvalsD[nnzD]=  -1.0/nb;
 
         nnzD += 1;
         nnzDrow += 1; 
       }      
       colptrD[nnodes+i + 1] = colptrD[nnodes+i] + nnzDrow;
     }


     if (0) 
     { 

       printf(" Sparse Matrix \n" );

	printf("%d \n", NDENSE );
	printf("%d \n", colptrD[NDENSE] );
        for (j = 1; j <= NDENSE+1; j++) 
		printf("%d \n", colptrD[j-1] -1);
        for (j = 1; j <= colptrD[NDENSE]; j++) 
		printf("%d \n", rowindD[j-1]-1 );

	printf(" Sparse & Dense Matrix \n" );

        for (i = 1; i <= NDENSE+1; i++) {
        for (j = colptrD[i-1]; j < colptrD[i]; j++) {
          printf("%d %d %e \n", i, rowindD[j-1], nzvalsD[j-1]);
         }      
       }
     }


     if (1) 
     { 
        FILE *mat_file;
        mat_file = fopen ("matrixRGF.mat", "w");

        fprintf (mat_file, "%d \n", NDENSE);
        fprintf (mat_file, "%d \n", NDENSE);
        fprintf (mat_file, "%d \n", colptrD[NDENSE+1-1]-1);
        for (j = 1; j <= NDENSE+1; j++) 
		fprintf(mat_file,"%d \n", colptrD[j-1]-1 );
        for (j = 1; j < colptrD[NDENSE]; j++)
		fprintf(mat_file,"%d \n", rowindD[j-1]-1 );
        for (j = 1; j <= colptrD[NDENSE]; j++) 
		fprintf(mat_file,"%24.16e \n", nzvalsD[j-1] );

         exit(9);
       
     }


   free(colptrD); 
   free(rowindD); 
   free(nzvalsD); 
   free(rhsD); 
   free(xsolD); 

   }





   free(mesh);
   free(colptr); 
   free(rowind); 
   free(nzvals); 
   free(nzvals2); 
   free(rhs); 
   free(xsol); 
   return 0;
}
