gfortran -c multigrid.f90     
gfortran -c relaxation.f90
gfortran -c gauss_seidel.f90  
gfortran -c residual.f90
gfortran -c l2norm.f90        
gfortran -c prolongation.f90  
gfortran -c restriction.f90
gfortran -c poisson3d.f90     

gfortran -o multigrid *.o
rm *.o
