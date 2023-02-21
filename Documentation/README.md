# Documentation

## 


## Further Details

This section is aimed to the interested user. This section gives just the essential theoretical details and outlines the functioning of the code without loosing too much time in technicalities.

For further details about the theory the user is referred to the papers (xxx yyy) and is encouraged to write to the corresponding author for any further clarifications.  

Regarding how the code works "under the hood" the following sections outiline the basic functioning and the comments in the code explain the remaining details.

### Theoretical Background
 
### Code Overview 
 The code works as follows:
 
1) The code reads the data from the files supplied by the user and reshapes the data in matrices, in this way it is easier to implement the derivaties.    

2) The code generates several points on the unit sphere. The points are generated using the Fibonacci sphere algorithm. Each point, along with the origin, defines a unit vector. The unit vectors generated are stored and used later.

3) The code picks a direction computed in point 2 and calculates the "ingredients" needed for the anaylis.  It computes the first and second deriatives of the dipole, of the fe vector, the hessian and the gradient of the original PES. (see eqs. XXX of YYY)

4) Maintaining the direction chosen in point 2 the code solves eq XXX of YYY for each point of the grid. The solution of this equation is an amplitude of the electric field. Now we have the amplitude of the field and the direction of the field. The code checks if in this point the pertubed gradient is "close to zero", that is its norm is below a treshold (defined by the user). If this requirement is met, the code proceeds, otherwise it goes back to step 3 and uses a different direction from the set of directions computed in point 2.

5) Now we have a direction of the field, the associated amplitude and we know that the pertubed gradient in this point is reasonably close to zero. A point must satisfy this condition to be a bond-breaking-point, but it is not a sufficient condition. To be in the **optimal** bond-breaking-point means that the perturbed gradient is zero, that the sigma function is zero and that the difference between the first and the third term of the sigma function is zero (see eq. xxx in yyy). Now the code checks if the two other conditions are satisfied. 

6) In the point in which the pertubed gradient is zero the code calculates the terms of the sigma function. (see eq. xxx in yyyy) then checks if the sigma function is below a treshold (used defined) and if the difference between the first and the third term of the sigma function is below the same treshold. If these criterion is met the coordinates of this point, the field and the sigma function terms are stored.

7) The code repeats points 3-6 for all the directions generated in point 2. Usually the code finds several candidates for the role of **optimal** bond-breaking point. The code now finds the point that is closest to the optimal bond breaking point between the candidates found in the previous points. The optimal bond-breaking-point is the point in which the sigma function is the smallest.

8) The code outputs two files: a) "fields_and_bbps.txt" contains all bond-breaking-points found by the code and b) "optimal_field_and_optimal_bbp.txt" contains the optimal field and the optimal bond-breaking-point.

9) Lastly the code plots the original PES and the PES perturbed by the optimal external electric field.



