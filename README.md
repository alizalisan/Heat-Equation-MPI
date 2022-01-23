# 2D Steady State Heat Equation using Iterative Method

Objective:

Write an MPI program to simulate the 2-dimensional (2D) steady-state "heat" equation using an iterative method.

Problem:

Consider a 2D "plate" of some material (e.g., steel) with the edges being held at fixed, constant temperatures (i.e., boundary values): Vleft, Vright, Vtop, Vbottom. Assume the rest of the plate starts with an initial temperature (at T=0): Vinit. The steady state heat equation is used to determine the temperature at any point on the plate at any future time.

One computational approach to solving for the heat values in the plate at a particular time is to use an iterative simulation method where time is advanced in steps (delta) and the current heat values at T are used to estimate the next heat values at T+delta. The equations used would be approximating the physics partial differential equation for heat transfer in the physical material. Because a physical plate is continuous, it is also necessary to discretize the plate "region" for purposes of the iterative computation.

Consider a 2D plate with dimensions M x N, where each 1 x 1 element represents a temparture value on the region at a location (i,j). There are MxN total elements that will need to be updated at each time step. Let an array V[M][N] how these values. The initial temperatures values are set as follows:

  V[0:M-1,0]     = Vleft
  V[0:M-1,N-1]   = Vright
  V[0,0:N-1]     = Vtop
  V[M-1,0:N-1]   = Vbottom
  V[1:M-2,1:N-2] = Vinit
Given plate temperatures at time T, V(T), the following equation is used to compute the temperature at time T+delta at element (i,j):
  V[i,j] = ( V[i-1,j] V[i+1,j] V[i,j-1] V[i,j+1] ) / 4
This equation uses a "4-point stencil" since it 4 neighboring elements to (i,j) determine its next value. A "5-point stencil" would include V[i][j] in the sum and would divide by 5.

To parallel the 2D heat equation using iteration as described above using message passing, it is necessary to "decompose" the region into subregions and assign those subregions to processes. For instance, suppose that the plate is discretized into a 1024x1024 region. If there are 16 processes (i.e., ranks)in the message passing program, the region (known as the "domain" of the simulation) could be partitioned into 16 subregions, each of size 256x256 elements (see below).

  +----------+----------+----------+----------+
  |          |          |          |          |
  |   P0     |   P1     |   P2     |   P3     |
  |          |          |          |          |
  +----------+----------+----------+----------+
  |          |          |          |          |
  |   P4     |   P5     |   P6     |   P7     |
  |          |          |          |          |
  +----------+----------+----------+----------+
  |          |          |          |          |
  |   P8     |   P9     |   P10    |   P11    |
  |          |          |          |          |
  +----------+----------+----------+----------+
  |          |          |          |          |
  |   P12    |   P13    |   P14    |   P15    |
  |          |          |          |          |
  +----------+----------+----------+----------+
At each time step, each process would update the elements in its subregion. (The process is considered to "own" these this elements and its responsibility to update them follows the "owner's compute" rule.) However, there is a minor difficulty to attend to. Let the following figure represent all of the elements (e) in a subregion that a process owns. For the elements at the edges, 1 or 2 elements necessary in the update equation are not owned by the process, but are owned by another process either "above" (A), "below" (B), "left" (L), or "right" (R) of the subregion. (Notice the boundary elements are special cases.) It is necessary for these elements to be provided by the neighboring processes.

   A A A A A A A A
  +-+-+-+-+-+-+-+-+
L |e|e|e|e|e|e|e|e| R
  +-+-+-+-+-+-+-+-+
L |e|e|e|e|e|e|e|e| R
  +-+-+-+-+-+-+-+-+
L |e|e|e|e|e|e|e|e| R
  +-+-+-+-+-+-+-+-+
L |e|e|e|e|e|e|e|e| R
  +-+-+-+-+-+-+-+-+
   B B B B B B B B
The elements that need to be communicated are called "ghost cells" and are essentially copies of the edge cells that are sent to the neighbors. The most straightforward way to implement the ghost cells is as extra rows and columns in the subregion data structure (see below). (Notice that the corner ghost cells are not neede for a 4-point stencil.)

  +-+-+-+-+-+-+-+-+-+-+
  | |G|G|G|G|G|G|G|G| |
  +-+-+-+-+-+-+-+-+-+-+
  |G|e|e|e|e|e|e|e|e|G|
  +-+-+-+-+-+-+-+-+-+-+
  |G|e|e|e|e|e|e|e|e|G|
  +-+-+-+-+-+-+-+-+-+-+
  |G|e|e|e|e|e|e|e|e|G|
  +-+-+-+-+-+-+-+-+-+-+
  |G|e|e|e|e|e|e|e|e|G|
  +-+-+-+-+-+-+-+-+-+-+
  | |G|G|G|G|G|G|G|G| |
  +-+-+-+-+-+-+-+-+-+-+
Thus, the procedure for each process at each update step is:

Send the edge elements to the appropriate neighbor processes.
Receive edge elements from neighbor processes and put them in the ghost cells.
Update the elements in the subregion.
This procedure can be executed for a fixed # timesteps or until a convergence criteria is reach. Convergence is typically determined by comparing how much the updated state differs from the previous state. This is done element by element. Given criteria C, if any element differs more then C, convergence is not reached. Convergence testing can be done per process, but determining global convergence across all processes requires finding out if any subregion is not converged.

Assignment:

Write an MPI program in C or C++ to implement the 2D heat equation using an iterative algorithm as described above. Demonstrate the program on a 1024x2048 region for 1, 2, 4, 8, 16, and 32 processes, running for 100 timesteps. Let all temperatures be floating point values in the range 0.0 to 1.0. Let the boundry values be 1.0. (You can change these values to get different outcomes.) Assume the rest of the plate starts with an initial temperature of 0.0. Print out the maximum element update difference for every subregion at every 5 timesteps.

A good way to test your code is to first write a sequential solution and make sure it is working correctly. Then this version will be your "gold standard" for comparing the parallel results at any timestep.
