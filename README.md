# standing_genetic_variation_and_de_novo_mutations
This repository contains one .cpp file:

(1) The code standing-denovo-github.cpp simulates the rescue process. Standing genetic variation is considered by assuming that the initial
phenotypic distribution for each trait is taken from a Gaussian distribution. 


Requires the GSL scientific library:

Ubuntu system installation through terminal:

(C++)	~ sudo apt-get install build-essential
(GSL)	~ sudo apt-get install gsl-bin
(GSL)	~ sudo apt-get install libgsl-dev

Ubuntu systems compilation through terminal:

~ c++ -O3 name.cpp -o [executable_name] -lm -lgsl -lgslcblas

To run the code, we must provide input data.
                
The input data for code (1) are:

     i- carrying capacity K
     ii - number of traits n
     iii - mutation probability 
     iv - number of independent runs   
     v - the strength of selection, which we set at 1
     vi - mutation effect size
     vii - stress level, delta
     viii - Maximum fitness W_max
    
     Here is an example of how to run the code in an Ubuntu terminal

     ./[executable_name] 10000 3 0.005 10000 1.0 0.1 0.3 1.5

    The output of the code is created in a.DAT file. The name of the file includes the values of the parameter 
    used. The first column of the DAT file is delta (drop in fitness), whereas the second column provides the 
    extinction probability, and the third column provides the rescue probability


