multitaperDevelopment
======================

Development code for multitaper package. I do not plan to put this on CRAN as a package but functionality from this may be moved to the multitaper package at some point.

No documentation provided, but this code should work when the multiaper package is installed. 

Build instructions for general Linux machine.

1) download and unzip to a folder called multitaperDevelopment-master
2) from the parent folder: 
 a) R CMD build multitaperDevelopment-master/
 b) R CMD check multitaperDevelopment_0.1-1.tar.gz (optional)
 c) R CMD INSTALL multitaperDevelopment_0.1-1.tar.gz

Note: The version number may change, and you will likely have to set the PATH variable for other operating systems.



