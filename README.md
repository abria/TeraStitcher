TeraStitcher
===========================================================

A tool for fast automatic 3D-stitching of teravoxel-sized 
microscopy images (BMC Bioinformatics 2012, 13:316)

Exploiting multi-level parallelism for stitching very large 
microscopy images (Frontiers in Neuroinformatics, 13, 2019)

===========================================================

Before using this software, you MUST accept the LICENSE.txt

Documentation,  help and  other info  are available on  our 
Github wiki at http://abria.github.io/TeraStitcher/.

===========================================================
Contributors

- Alessandro Bria (email: a.bria@unicas.it).
  Post-doctoral Fellow at University of Cassino (Italy).
  Main developer.

- Giulio Iannello (email: g.iannello@unicampus.it).
  Full Professor at University Campus Bio-Medico of Rome (italy).
  Supervisor and co-developer.
  
===========================================================
Main features

- designed for images exceeding the TeraByte size
- fast and reliable 3D stitching based on a multi-MIP approach
- typical memory requirement below 4 GB (8 at most)
- 2D stitching (single slice images) supported
- regular expression based matching for image file names
- data subset selection
- sparse data support
- i/o plugin-based architecture
- stitching of multi-channel images
- support for big tiff files (> 4 GB)
- HDF5-based formats
- parallelization on multi-core platform
- fast alignment computation on NVIDIA GPUs

===========================================================
