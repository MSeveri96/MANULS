 Briefly speaking, the code works as follows:
 
 1) The code reads the data from the files supplied by the user and reshapes the data in matrices, in this way it is easier to implement the derivaties.    
 
2) The code generates several points on the unit sphere. The points are generated using the Fibonacci sphere algorithm. Each point, along with the origin, defines a unit vector. The unit vectors generated are stored and used later.

3) The code picks a direction computed in point 2 and calculates the "ingredients" needed for the anaylis.  It computes the first and second deriatives of the dipole, of the fe vector, the hessian and the gradient of the original PES. (see eqs. XXX of YYY)

4) Maintaining the direction chosen in point 2 the code solves eq XXX of YYY for each point of the grid. The solution of this equation is an amplitude of the electric field. Now we have the amplitude of the field and the direction of the field. The code checks if in this point the pertubed gradient is "close to zero", that is its norm is below a treshold. If this requirement is met, the code proceeds, otherwise it goes back to step 3 and uses a different direction from the set of directions computed in point 2.

5) Now we have a direction of the field, the associated amplitude and we know that the 
pertubed gradient in this point is reasonably close to zero. 
