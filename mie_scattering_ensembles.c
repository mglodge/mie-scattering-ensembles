/*

This code uses Mie theory to calculate scattering phase functions for ensembles of particles in air.

    INPUTS:

        - refractive index of material

        - wavelength of light used (units: um)

        - "size_distribution.txt" 

            -  Size distribution of particles, using the direct output file from ExCESS (Exoplanet Cloud Ensemble Scattering System)
            -  Lines 1-16: Header information
            -  Line 17: Number of bins
            -  Line 18: onwards: bin_start_diameter, bin_end_diameter, number of particles in bin

    OUTPUTS:

        - "S1_and_S2_vs_angle.txt"

            -  .csv file containing: 
                - angle (degrees)
                - average ensemble scattering amplitude S1^2
                - average ensemble scattering amplitude S2^2
                - total unpolarised average scattering magnitude
                - degree of linear polarisation

        - Tables of data for each size bin, including:

            - Q_sca (scattering efficiency)
            - Q_ext (extinction efficiency)
            - Q_abs (absorption efficiency)
            - S1 (cumulative scattering amplitude at an angle of 0 degrees)
            - S2 (cumulative scattering amplitude at an angle of 0 degrees)

        - Phase function graph: Total unpolarised scattering intensity of ensemble vs angle

        - Phase function graph: Polarised components S1 and S2, as well as the average 

        - Graph of 'Degree of linear polarisation vs angle'


                                                                    Written by Matt Lodge (31/01/24)

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


/* =============== FILE HANDLES =============== */
FILE *sizedistributionfile = NULL;          /* size distribution input file */
FILE *outfile = NULL;                     /* angular scattering output file */

/* =============== LOOP COUNTERS / INDICES =============== */
int i = 0;                                  /* generic loop index */
int j = 0;                                  /* generic loop index */
int k = 0;                                  /* generic loop index */
int n = 0;                                  /* Mie series index */
int n_max = 0;                              /* max Mie order */
int number_of_wavelengths = 0;              /* number of wavelengths */
int N_bins = 0;                             /* number of size distribution bins */
int r_index = 0;                            /* radius/bin index */

/* =============== CONSTANTS =============== */
const double pi = 3.14159265358979323846264;

/* =============== PHYSICAL PARAMETERS =============== */
double wavelength = 0.0;                    /* wavelength (µm) */
double max_wavelength = 0.0;                /* max wavelength */
double min_wavelength = 0.0;                /* min wavelength */
double wavelength_increment = 0.0;          /* wavelength increment (um) */
double sphere_radius = 0.0;                 /* particle radius */
double magnetic_perm_medium = 0.0;          /* magnetic permeability of the medium */
double magnetic_perm_sphere = 0.0;          /* magnetic permeability of the medium */
double x_size_parameter = 0.0;              /* size parameter (2*pi*R/lambda) */
double theta = 0.0;                         /* scattering angle (rad) */
double N_particles = 0.0;                   /* number of particles in a given bin */

/* =============== MIE EFFICIENCIES =============== */
double Q_scatter = 0.0;                     /* scattering efficiency */
double Q_extinction = 0.0;                  /* extinction efficiency */
double Q_absorption = 0.0;                  /* absorption efficiency */

/* =============== OUTPUTS (each with 181 values) =============== */
double total_magnitudes[181] = {0};        /* intensity */
double DOLP_magnitudes[181] = {0};         /* degree of linear polarisation */
double average_S1_squared[181] = {0};      /* <|S1|^2> */
double average_S2_squared[181] = {0};      /* <|S2|^2> */

/* =============== VARIABLES FOR USE IN BESSEL FUNCTION/SCATTERING CALCULATIONS =============== */
double sum_scatter = 0.0; // cumulative sum of scattering terms
double sum_extinction = 0.0; // cumulative sum of extinction terms
double magnitude_S1 = 0.0; // magnitude of scattering amplitude S1
double magnitude_S2 = 0.0; // magnitude of scattering amplitude S2

/* =============== COMPLEX VARIABLES (in Mie calculations) =============== */
double _Complex ref_ind_sphere = 0.0;       /* refractive index sphere */
double _Complex diff_xjx = 0.0;             /* differential of x * j_n(x) */
double _Complex diff_xhx = 0.0;             /* differential of x * h_n(x) */
double _Complex diff_mxjmx = 0.0;           /* differential of mx * j_n(mx) */
double _Complex a_coeff = 0.0;              /* Mie a_n */
double _Complex b_coeff = 0.0;              /* Mie b_n */
double _Complex sum_S1 = 0.0;               /* S1 sum */
double _Complex sum_S2 = 0.0;               /* S2 sum */
double _Complex S1[181] = {0};              /* scattering amplitude S1 */
double _Complex S2[181] = {0};              /* scattering amplitude S2 */

/* =============== DYNAMIC ARRAYS (size determined during runtime) =============== */
double *bessel_first_kind_x = NULL;                 /* j_n(x) */
double *bessel_second_kind_x = NULL;                /* y_n(x) */
double *pi_function = NULL;                         /* pi_n */
double *tau_function = NULL;                        /* tau_n */
double _Complex *bessel_first_kind_mx = NULL;       /* j_n(mx) */
double _Complex *complex_hankel_function = NULL;    /* h_n(x) */

/* ================= SIZE DISTRIBUTION (size determined during runtime)  ================= */
double **input_size_distribution = NULL;    /* raw input */
double **size_distribution = NULL;          /* processed */
double ***S1_S2_squared_values = NULL;      /* |S1|^2, |S2|^2 storage */


int main()
{

    /* ------------------ Set initial parameters ------------------ */

    ref_ind_sphere = 1.49 + 0*I; // remember -- these are complex
    wavelength= 0.532; // um

    /* ------------------------------------------------------------ */


    /* these values are usually constant, but you can alter them if required */

    magnetic_perm_medium= 0.0000012566; // 1.26 uH/m
    magnetic_perm_sphere= 0.0000012566; // 1.26 uH/m

    /* read in size distribution file */
    sizedistributionfile = fopen("size_distribution.txt","r");
    if (sizedistributionfile == NULL) {
        printf("\n\nError- the shape file cannot be found!! \n\n\n");
        return 1;
    }

    for(int i=0;i<16;i++){
        fscanf(sizedistributionfile, "%*[^\n]\n"); // skip the 16 header lines
    }
    fscanf(sizedistributionfile, " %d", &N_bins); // record the number of bins

    /* initialise matrices to store the size distribution and computed quantities for S1/S2, now that we know how much data there is */

    input_size_distribution=(double**)malloc(N_bins*sizeof(double*));
    for (i=0; i<N_bins; i++) {
        input_size_distribution[i]=(double*)malloc((3)*sizeof(double));   // matrix will record: d_min, d_max (both diameters), number for each row of the size distribution
    }

    size_distribution=(double**)malloc(N_bins*sizeof(double*));
    for (i=0; i<N_bins; i++) {
        size_distribution[i]=(double*)malloc((2)*sizeof(double));   // matrix will record: r_average (average radius), number of particles at this radius (for each row of the size distribution), S1, S2, S_total (all magnitudes of the complex numbers for the intensities)
    }

    S1_S2_squared_values=(double***)malloc(N_bins*sizeof(double**));
    for (i=0; i<N_bins; i++) {
        S1_S2_squared_values[i]=(double**)malloc(181*sizeof(double*));
        for (j=0; j<181; j++) {
            S1_S2_squared_values[i][j]=(double*)malloc((2)*sizeof(double));  //makes a 3D matrix with index [a][b][c] where [a]=size distribution radius index, [b]=angle index (0 -> 180 degrees) and [c]= 2 elements [0,1] = S1, S2
        }
    }

    /* read in the actual size distribution data from the file */
    for(i=0;i<N_bins;i++){
        fscanf(sizedistributionfile, " %lf %lf %lf", &input_size_distribution[i][0], &input_size_distribution[i][1], &input_size_distribution[i][2]); // read in d_min, d_max, number for each row
    }
    fclose(sizedistributionfile);

    printf("\n SIZE DISTRIBUTION DATA IMPORTED:\n");
    printf("\n	 Diameter start (um) 	 Diameter end (um)     Number of particles ");
    for(i=0;i<N_bins;i++){
        printf("\n\t      %lf                 %lf              %lf", input_size_distribution[i][0], input_size_distribution[i][1], input_size_distribution[i][2]);
    }

    /* rearrange the data so that we only store the middle value of the bin size and convert it to radius */

    printf("\n\n SIZE DISTRIBUTION USING THE MEAN RADIUS OF EACH BIN:\n");
    printf("\n	 Mean radius (um) 	 Number of particles");
    N_particles=0;

    for(i=0;i<N_bins;i++){
        size_distribution[i][0]= (input_size_distribution[i][0]+input_size_distribution[i][1])/2; //find average bin diameter
        size_distribution[i][0]= size_distribution[i][0]/2; // convert diameter to radius
        size_distribution[i][1]= input_size_distribution[i][2]; //number stays the same
        printf("\n\t     %lf                  %lf ", size_distribution[i][0], size_distribution[i][1]);

        N_particles += input_size_distribution[i][2]; //add together to find the total number of particles in the distribution
    }

    printf("\n\n A total of %lf particles were found in %d bins.\n", N_particles, N_bins);
    
            // -------------------- MIE CODE BEGINS HERE -------------------- 

    // print header for table of key results
    printf("\n %15s %15s %20s %8s %15s %15s %15s %25s %30s\n","Wavelength (um)","Radius (um)","x (Size Parameter)","nmax","Qscatter","Qextinction","Qabsorption","Sum[S1(theta=0)]", "Sum[S2(theta=0)]");

    for(r_index=0; r_index<N_bins; r_index++){  /* loop Mie code over each radius in the size distribution */
 
        sphere_radius=size_distribution[r_index][0]; //set sphere radius as each value from our size distribution

        /* calculate size parameter and maximum number of terms in summation equation for this radius and wavelength */

        x_size_parameter= 2*pi*sphere_radius/wavelength;
        n_max= ceil(x_size_parameter + 4*pow(x_size_parameter, (1.0/3.0)) +2);

        /* Declare dynamic arrays */

        bessel_first_kind_x = (double*)malloc(n_max*sizeof(double));
        bessel_first_kind_mx = (double _Complex*)malloc(n_max*sizeof(double _Complex));
        bessel_second_kind_x = (double*)malloc(n_max*sizeof(double));
        pi_function = (double*)malloc(n_max*sizeof(double));
        tau_function = (double*)malloc(n_max*sizeof(double));
        complex_hankel_function = (double _Complex*)malloc(n_max*sizeof(double _Complex));

        for(k=0; k<=180; k++){ //run between all angles 0 -> +180

            theta=k*pi/180; // change theta every cycle (and convert from deg -> radians)

            /* Set initial values for computing Bessel, Hankel, pi and tau functions as a series (we compute them from the two intial terms in each case) */

            bessel_first_kind_x[0]= sin(x_size_parameter)/x_size_parameter;
            bessel_first_kind_x[1]= sin(x_size_parameter)/pow(x_size_parameter,2) - cos(x_size_parameter)/x_size_parameter;

            bessel_second_kind_x[0]= -cos(x_size_parameter)/x_size_parameter;
            bessel_second_kind_x[1]= -cos(x_size_parameter)/pow(x_size_parameter,2) - sin(x_size_parameter)/x_size_parameter;

            bessel_first_kind_mx[0]= csin(ref_ind_sphere*x_size_parameter)/(ref_ind_sphere*x_size_parameter); //these functions are complex too
            bessel_first_kind_mx[1]= csin(ref_ind_sphere*x_size_parameter)/cpow(ref_ind_sphere*x_size_parameter,2) - ccos(ref_ind_sphere*x_size_parameter)/(ref_ind_sphere*x_size_parameter);

            complex_hankel_function[0]= bessel_first_kind_x[0] + bessel_second_kind_x[0]*I;
            complex_hankel_function[1]= bessel_first_kind_x[1] + bessel_second_kind_x[1]*I;

            pi_function[0]= 0;
            pi_function[1]= 1;

            tau_function[0]= 0;
            tau_function[1]= cos(theta);

            /* initialise summation terms */
            
            sum_scatter=0;
            sum_extinction=0;
            sum_S1=0;
            sum_S2=0;
            for(i=1;i<n_max;i++){ // Calculate summation terms for Mie coefficients

                if(i>1){ // Calculate new function terms (only after the first loop -- we already have terms 0 and 1 from initialisation above)
                    bessel_first_kind_x[i]= ((2*(i-1)+1)/x_size_parameter)*bessel_first_kind_x[i-1] - bessel_first_kind_x[i-2];
                    bessel_first_kind_mx[i]= ((2*(i-1)+1)/(ref_ind_sphere*x_size_parameter))*bessel_first_kind_mx[i-1] - bessel_first_kind_mx[i-2];
                    bessel_second_kind_x[i]= ((2*(i-1)+1)/x_size_parameter)*bessel_second_kind_x[i-1] - bessel_second_kind_x[i-2];
                    complex_hankel_function[i]=((2*(i-1)+1)/x_size_parameter)*complex_hankel_function[i-1] - complex_hankel_function[i-2];

                    pi_function[i]= ((2*i-1.0)/(i-1.0))*cos(theta)*pi_function[i-1] - (i/(i-1.0))*pi_function[i-2];
                    tau_function[i]= i*cos(theta)*pi_function[i] - (i+1)*pi_function[i-1];
                }

                /* Calculate coefficients a_n and b_n */

                diff_xjx= x_size_parameter*bessel_first_kind_x[i-1] - i*bessel_first_kind_x[i];
                diff_mxjmx= ref_ind_sphere*x_size_parameter*bessel_first_kind_mx[i-1] - i*bessel_first_kind_mx[i];
                diff_xhx = x_size_parameter*complex_hankel_function[i-1] - i*complex_hankel_function[i];

                a_coeff= (magnetic_perm_medium*cpow(ref_ind_sphere,2)*bessel_first_kind_mx[i]*diff_xjx - magnetic_perm_sphere*bessel_first_kind_x[i]*diff_mxjmx) / (magnetic_perm_medium*cpow(ref_ind_sphere,2)*bessel_first_kind_mx[i]*diff_xhx - magnetic_perm_sphere*complex_hankel_function[i]*diff_mxjmx);
                b_coeff= (magnetic_perm_sphere*bessel_first_kind_mx[i]*diff_xjx - magnetic_perm_medium*bessel_first_kind_x[i]*diff_mxjmx) / (magnetic_perm_sphere*bessel_first_kind_mx[i]*diff_xhx - magnetic_perm_medium*complex_hankel_function[i]*diff_mxjmx);
                
                // calculate summation terms for scattering and extinction
                sum_scatter= sum_scatter + (2.0*i+1.0)*(pow(creal(a_coeff),2) + pow(cimag(a_coeff),2) + pow(creal(b_coeff),2) + pow(cimag(b_coeff),2));
                sum_extinction= sum_extinction + (2.0*i+1.0)*(creal(a_coeff+b_coeff));
                
                // calculate summation terms for S1 and S2 (scattering amplitudes)
                sum_S1= sum_S1 + ((2.0*i+1)/(i*(i+1.0)))*(a_coeff*pi_function[i] + b_coeff*tau_function[i]);
                sum_S2= sum_S2 + ((2.0*i+1)/(i*(i+1.0)))*(a_coeff*tau_function[i] + b_coeff*pi_function[i]);

            }
            
            //determine Mie Q_efficiencies
            Q_scatter= (2.0/pow(x_size_parameter,2))*sum_scatter;      //Q_scatter
            Q_extinction = (2.0/pow(x_size_parameter,2))*sum_extinction;  //Q_extinction
            Q_absorption = Q_extinction - Q_scatter;  //Q_absorption

            if(k==0){ // for angle = 0 degrees only, print key terms
                printf("%10.3f  %16.3f   %16.5f  %12d    %14.5e %14.5e %15.5e   %12.6f + %12.6fi     %12.6f + %12.6fi\n",  wavelength, sphere_radius, x_size_parameter, n_max, Q_scatter, Q_extinction, Q_absorption, creal(sum_S1), cimag(sum_S1), creal(sum_S2), cimag(sum_S2));
            }

            S1[k]=sum_S1; // store perpendicular polarisation intensity
            S2[k]=sum_S2; // store parallel polarisation intensity

        }

        // save intensities for all angles at this radius into S1_S2_squared_values matrix
    
        for(k=0; k<=180; k++){ //run between all angles 0 -> +180
            
            // find the magitudes |S_1| and |S_2| for each angle and particle size
            magnitude_S1 = sqrt(pow(creal(S1[k]),2)+pow(cimag(S1[k]),2));
            magnitude_S2 = sqrt(pow(creal(S2[k]),2)+pow(cimag(S2[k]),2));

            // square the magnitudes to find |S_1|^2 and |S_2|^2 and store in matrix
            S1_S2_squared_values[r_index][k][0] = pow(magnitude_S1,2); 
            S1_S2_squared_values[r_index][k][1] = pow(magnitude_S2,2);
        }

        /*free memory for arrays for each radius/wavelength combination*/
        free((void*) bessel_first_kind_x);
        free((void*) bessel_first_kind_mx);
        free((void*) bessel_second_kind_x);
        free((void*) pi_function);
        free((void*) tau_function);
        free((void*) complex_hankel_function);

        //return to top of radius loop and use next particle radius in the distribution    
    }


    /* calculate average intensities for each angle using the size distribution */

    for(k=0; k<=180; k++){ // initialise arrays for all angles between 0 -> +180
        average_S1_squared[k]=0;
        average_S2_squared[k]=0;
    }

    for(k=0; k<=180; k++){ //calculate bin-weighted averages for all angles between 0 -> +180
        for(r_index=0; r_index<N_bins; r_index++){
            average_S1_squared[k] += S1_S2_squared_values[r_index][k][0]*size_distribution[r_index][1]; // add together all of the |S_1|^2 values, multiplied by the number of particles in this bin, for each radius
            average_S2_squared[k] += S1_S2_squared_values[r_index][k][1]*size_distribution[r_index][1]; // add together all of the |S_2|^2 values, multiplied by the number of particles in this bin, for each radius
        }

        //divide by the total number of particles to find the average
        average_S1_squared[k] = average_S1_squared[k] / N_particles;
        average_S2_squared[k] = average_S2_squared[k] / N_particles;
    }

    // find total unpolarised magnitudes and DOLP of the averages for each angle
    for(k=0; k<=180; k++){
        total_magnitudes[k]= (average_S1_squared[k]+average_S2_squared[k])/2; //finding the unpolarised average of the S1+S2 averages right at the end (the mean resulting magnitude)
        DOLP_magnitudes[k]= 100.0*(average_S1_squared[k]-average_S2_squared[k])/(average_S1_squared[k]+average_S2_squared[k]); // x 100.0 to express as a percentage: finds DOLP from the final averaged values of S1 and S2
    }

    /* save results to outfile to print a graph of the averages */
    outfile=fopen("S1_and_S2_vs_angle.txt","w");
    for(k=0; k<=180; k++){ //run between all angles 0 -> +180
        fprintf(outfile,"%d,%f,%f,%f,%f\n", k, average_S1_squared[k], average_S2_squared[k], total_magnitudes[k], DOLP_magnitudes[k]);
    }    
    fclose(outfile);
    
    /* load python plotting tool and graph results */
    printf("\n\n Loading mie_plot.py...\n");
    system("python mie_plot.py"); // run matplotlib to produce a graph of the extracted results

    printf("\n Code complete!\n\n");

    /* free data stored in arrays */

    for (i=0; i<N_bins; i++) {
        free((void*)input_size_distribution[i]);
    }
    free((void*)input_size_distribution);

    for (i=0; i<N_bins; i++) {
        free((void*)size_distribution[i]);
    }
    free((void*)size_distribution);

    for (i=0;i<N_bins;i++) {
        for(j=0;j<181;j++){
            free((void*)S1_S2_squared_values[i][j]);
        }
        free((void*)S1_S2_squared_values[i]);
    }
    free((void*)S1_S2_squared_values);

    return 0;
}

