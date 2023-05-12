This is a bunch of images that works without any modify of codes in this repisitory.

Images are provided by Katie Bunde in Matrix Mechanotransduction Laboratory at Boston University. 

Dot Tracking (MATLAB)
-Put all files under the same file holder.

-Rename all images as "...01", "...02" and so on.

-If the images have a drift, run "registration.m" to get rid of the drift.

-Run "pattern_generation.m" to create a reference image.

-Run "Dot_Tracking_Part1.m".

-Run "Dot_Tracking_Part2.m".

-Run "traction.m".

-Run "run_distmesh.m" to create a boundary for a cell.

-Run "distmesh_2d" to generate a mesh for a cell.

Traction forces and modulus distribution (FEniCS) 
-Put "cell.xml" and "celldata.txt" created in Dot tracking part in MATLAB under the same file holder.

-Run "cell_homogeneous.py" in FEniCS to compute the traction forces and modulus distributions for a homogeneous cell.

-Run "cell_heterogeneous.py" in FEniCS to compute the traction forces and modulus distributions for a heterogeneous cell.
