//  Reconstruction for 2D and 3D
//  J. Nelson & J. Steinbrener

// Currently this program: 
// Opens the given script file, reads h5 filename and reconstruction 
// parameters, opens and reads h5 file, creates an object with random 
// values as itn_array, performs difference map, and updates itn_array.
// Note that for difference map:
//  - opt_itn_cas = Support estimate
//  - final_itn_cas = Modulus estimate. This will be saved to the file 
//    in the end or used for averaging!

// To run program with "x" number of processes:
// mpirun -np x dm_recon script_name


#include "dm_recon.h"

#define STRLEN 128

int main(int argc, char *argv[]) {
  char this_arg[STRLEN], error_string[STRLEN];
  char scriptname[STRLEN], filename[STRLEN];
  char png_name[STRLEN];
  reconParams recon_params;
  allocatedArrays allocated_arrays;
  dm_itn_struct itn_struct;
  dm_spt_struct spt_struct;
  /* Adi error array is converted to absolute error for magnitudes in the 
   * routine dm_recon_get_adm
   */
  dm_array_real_struct mag_error_ras; 
  dm_array_real_struct magnitude_ras;
  dm_array_complex_struct current_itn_cas, avg_final_itn_cas;
  dm_array_complex_struct final_itn_cas, opt_itn_cas;
  dm_array_byte_struct spt_bas;
  dm_array_real_struct recon_errors;
  int i_arg,it,avg_count;
  int debugWait;
  int my_rank, p;
  float error_i;
 
  dm_init(&p,&my_rank);

  /* define hardwired defaults */
  i_arg = 1;
  avg_count = 0;
  allocated_arrays.adi_error_allocated = 0;
  allocated_arrays.magnitude_allocated = 0;
  allocated_arrays.spt_allocated = 0;
  allocated_arrays.current_itn_allocated = 0;
  allocated_arrays.avg_final_itn_allocated = 0;
  allocated_arrays.final_itn_allocated = 0;
  allocated_arrays.final_itn_allocated = 0;
  error_i = 0.;
  
  /* define defaults that user can change through CLA */
  debugWait = 0;
  
  while (i_arg < argc) {
    strcpy( this_arg, argv[i_arg] );
    if ((strncasecmp("-?",this_arg,2) == 0) || 
	(strncasecmp("-H",this_arg,2) == 0)) {
      exit(1);
    } else if (strncasecmp("-W",this_arg,2) == 0) {
      debugWait = 1;
      i_arg++;
    } else {
      strcpy(scriptname,this_arg);  // reads script name
      i_arg++;
    }
  }
  
  while (debugWait);
  
  /* --------------------------------------------------------------- 
     ----------------- reading out script file ---------------------
     --------------------------------------------------------------- */ 
    if (dm_recon_read_script(scriptname, filename, &recon_params, 
			     error_string) < 0) {
      printf("%s", error_string);
      exit(1);
    } 
  
  /* --------------------------------------------------------------- 
     -------------------- Initialize errors ------------------------
     --------------------------------------------------------------- */ 
    if (my_rank == 0) {
      recon_errors.nx = recon_params.iterate_count;
      recon_errors.ny = 1;
      recon_errors.nz = 1;
      recon_errors.npix = recon_errors.nx;
      recon_errors.real_array = (dm_array_real *)
	malloc(recon_params.iterate_count*sizeof(dm_array_real));
    }

  
  /* ---------------------------------------------------------------  
     ----------- Get adi and adi_error arrays from file ------------
     ------------------- check total memory impact -----------------
     ------ crop 2D files if necessary and return magnitudes -------
     ------------- set nx, ny, and nz in recon_params --------------
     --------------------------------------------------------------- */
  if (dm_recon_get_adm(filename, &allocated_arrays, 
		       &mag_error_ras,
		       &magnitude_ras, 
		       &recon_params, error_string, 
		       my_rank, p) < 0) {
    printf("%s", error_string);
    exit(1);
  }


  /* --------------------------------------------------------------- 
     ------------- Get the existing spt_struct --------------------- 
     --------------------------------------------------------------- */
  if (dm_recon_get_spt(filename, &recon_params,
		       &spt_struct,
		       &spt_bas, 
		       &allocated_arrays, 
		       error_string, my_rank, p) < 0) {
    printf("%s", error_string);
    exit(1);
  }

  
  /* --------------------------------------------------------------- 
     ------------ Allocate memory for all complex arrays ----------- 
     --------------------------------------------------------------- */
  dm_recon_init_complex_arrays(&current_itn_cas, 
			       &final_itn_cas,
			       &opt_itn_cas,
			       &avg_final_itn_cas,
			       &allocated_arrays, &recon_params,
			       error_string,p);
  
  /* --------------------------------------------------------------- 
     ----- Create new initial itn array or use previously saved----- 
     --------------------------------------------------------------- */
  dm_recon_init_itn(filename,error_string,
                    &recon_params,&current_itn_cas, &spt_bas, 
		    &magnitude_ras,p,my_rank);
  
  /* ---------------------------------------------------------------
     ------------------ Main reconstruction loop -------------------
     ---------------------------------------------------------------  */
  dm_array_fft(&final_itn_cas,p,DM_ARRAY_CREATE_FFT_PLAN 
	       | DM_ARRAY_FFT_MEASURE, my_rank);
  printf("Created FFT plan.\n");

  /* current_itn needs plan for shrinkwrap! */
  dm_array_fft(&current_itn_cas,p,DM_ARRAY_CREATE_FFT_PLAN 
	       | DM_ARRAY_FFT_MEASURE, my_rank);
  printf("Created 2nd FFT plan.\n");

  if ((recon_params.diffmap) || (recon_params.shrink_every != 0)) {
    dm_array_fft(&opt_itn_cas,p,DM_ARRAY_CREATE_FFT_PLAN 
		 | DM_ARRAY_FFT_MEASURE, my_rank);  
    printf("Created 3rd FFT plan.\n");
  }

  for (it=0; it<recon_params.iterate_count; it++){
    
    /* -------------------------------------------------------------
       ---------- Perform reconstruction & calulate error ----------
       -------------------------------------------------------------  */ 
    if (recon_params.diffmap) {
      dm_recon_do_diffmap(&recon_params, &magnitude_ras,
			  &mag_error_ras, &current_itn_cas,
			  &final_itn_cas, &opt_itn_cas, 
			  &avg_final_itn_cas, &spt_bas, 
			  &recon_errors,p,my_rank,it,&avg_count,
                          filename);
    } else if (recon_params.hio) {
      dm_recon_do_hio(&recon_params, &magnitude_ras,
		      &mag_error_ras, &current_itn_cas,
		      &final_itn_cas, &opt_itn_cas, 
		      &avg_final_itn_cas, &spt_bas, 
		      &recon_errors, p,my_rank,it,&avg_count,
                      filename);
    }
  } /* endfor */
    
  itn_struct.iterate_count = it;
  printf("Finished %u iteration FFT loop. \n", it);
  
  if (recon_params.avg_every != 0) {
      dm_array_multiply_real_scalar(&avg_final_itn_cas, 
				    1.0/(dm_array_real)avg_count);
  }

  /* destroy the plans */
  dm_array_fft(&final_itn_cas,p,DM_ARRAY_DESTROY_FFT_PLAN,
               my_rank);
  dm_array_fft(&current_itn_cas,p,DM_ARRAY_DESTROY_FFT_PLAN,
               my_rank);
  if (recon_params.diffmap) {
    dm_array_fft(&opt_itn_cas,p,DM_ARRAY_DESTROY_FFT_PLAN,
		 my_rank);
  }

  /* ---------------------------------------------------------------  
     ------------------------ Update hdf5 --------------------------
     --------------------------------------------------------------- */ 
  if (recon_params.avg_every == 0) {
    if (dm_recon_save_recon(filename, &recon_params, &itn_struct, 
			    &final_itn_cas, &spt_struct, &spt_bas, 
			    &recon_errors,error_string,
			    my_rank, p) < 0) {
      printf(error_string);
      exit(1);
    }
  } else {
    if (dm_recon_save_recon(filename, &recon_params, &itn_struct, 
			    &avg_final_itn_cas,&spt_struct, &spt_bas,
			    &recon_errors, error_string,
			    my_rank, p) < 0) {
      printf(error_string);
      exit(1);
    }
  }
  
  /* ---------------------------------------------------------------
     --------------------------- free memory ----------------------- 
     --------------------------------------------------------------- */
  dm_recon_free_mem(&allocated_arrays,&mag_error_ras,
		    &magnitude_ras,
		    &spt_bas, &current_itn_cas,
		    &avg_final_itn_cas, &final_itn_cas,
		    &opt_itn_cas,&recon_errors,my_rank);
  
  
  dm_exit();
  exit(0);
}     

/*------------------------------------------------------------*/
int dm_recon_read_script(char *scriptname, 
			 char *filename,
			 reconParams *these_params, 
			 char *error_string)
{

  char current_char;
  FILE *pFile;

  if (strlen(scriptname) == 0) {
    sprintf(error_string,"dm_recon_read_script: No scriptfile given!\n");
    return(-1);
  }
  
  pFile = fopen(scriptname, "r");
  if (!pFile) {
    sprintf(error_string,"dm_recon_read_script: Error opening file %s. \n",
	    scriptname);
    return(-1);
  }
  
  /*    filename     */
  fscanf(pFile, "%s", filename);
  if (strlen(filename) == 0) {
    sprintf(error_string,"dm_recon_read_script: No filename specified!");
    return(-1);
  }
  printf("File to use \"%s\". \n", filename);
  
  /*      HIO       */
  current_char = fgetc(pFile);
  while (current_char != '=') {
    current_char = fgetc(pFile);
  }
  fscanf(pFile,"%d",&(these_params->hio));     // also a pointer
  if (these_params->hio) {
    these_params->diffmap = 0;
    printf("Using HIO\n");  
  } else {
    these_params->diffmap = 1;
    printf("Using DIFFMAP\n");  
  }


  /*      beta       */
  current_char = fgetc(pFile);
  while (current_char != '=') {
    current_char = fgetc(pFile);
  }
  fscanf(pFile,"%f",&(these_params->beta));     // also a pointer
  printf("Beta is equal to %f \n", these_params->beta);  
  
  these_params->gamma_mod = -1/(these_params->beta);
  these_params->gamma_spt = 1/(these_params->beta);
  printf("The modulus and support gamma are equal to %f and %f. \n",
	 these_params->gamma_mod, these_params->gamma_spt);
  
  /*    iterations     */
  current_char = fgetc(pFile);
  while (current_char != '=') {
    current_char = fgetc(pFile);
  }
  fscanf(pFile,"%d",&(these_params->iterate_count));
  printf("Number of iterations to perform: %u \n", these_params->iterate_count);

  /*    use saved iterate as starting guess?     */
  current_char = fgetc(pFile);
  while (current_char != '=') {
    current_char = fgetc(pFile);
  }
  fscanf(pFile,"%d",&(these_params->use_saved_iterate));
  if (these_params->use_saved_iterate == 1) {
      printf("Use previously saved iterate as starting guess if available.\n");
  } else {
      printf("Start with new starting guess.\n");
  }

  /*    use saved error array for modulus constraint?     */
  current_char = fgetc(pFile);
  while (current_char != '=') {
    current_char = fgetc(pFile);
  }
  fscanf(pFile,"%d",&(these_params->use_error_array));
  if (these_params->use_error_array == 1) {
      printf("Use magnitudes with uncertainties if available.\n");
  } else {
      printf("Assume exact magnitudes.\n");
  }

  /*     positivity    */
  current_char = fgetc(pFile);
  while (current_char != ':') {
    current_char = fgetc(pFile);
  }
  current_char = fgetc(pFile);
  while (current_char != 'y' && current_char != 'Y' &&
         current_char != 'n' && current_char != 'N') {
    current_char = fgetc(pFile);
  }

  if (current_char == 'y' || current_char == 'Y') {
    these_params->positivity = 1;
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    
    fscanf(pFile,"%d",&(these_params->imaginary_only));
    if (these_params->imaginary_only == 1) {
      printf("Use positivity constraint on imaginary part only. \n");
    } else {
      printf("Use positivity constraint. \n");
    }
  } else {
    these_params->positivity = 0;
    printf("No positivity constraint. \n");
  }
  
  /*     averaging     */
  current_char = fgetc(pFile);
  while (current_char != ':') {
    current_char = fgetc(pFile);
  }
  current_char = fgetc(pFile);
  while (current_char != 'y' && current_char != 'Y' && 
	 current_char != 'n' && current_char != 'N') {
    current_char = fgetc(pFile);
  }
  if (current_char == 'y' || current_char == 'Y') {
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%d",&(these_params->avg_after));
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%d",&(these_params->avg_every));
    printf("Average every %u iterations after the first %u.\n", 
	   these_params->avg_every, these_params->avg_after);
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%d",&(these_params->avg_save_iterates));
    if (these_params->avg_save_iterates == 1) {
        printf("Save each iterate before averaging.\n");
    }
  } else {
    these_params->avg_after = these_params->iterate_count;
    these_params->avg_every = 0;
    printf("Don't average.\n");
  }
  /* Make sure the values make sense */
  if (these_params->avg_after >= these_params->iterate_count) {
    these_params->avg_after = these_params->iterate_count;
    these_params->avg_every = 0;
  }

  /*    shrinkwrap     */
  current_char = fgetc(pFile);
  while (current_char != ':') {
    current_char = fgetc(pFile);
  }
  current_char = fgetc(pFile);
  while (current_char != 'y' && current_char != 'Y' && 
	 current_char != 'n' && current_char != 'N') {
    current_char = fgetc(pFile);
  }
  if (current_char == 'y' || current_char == 'Y') {
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%d",&(these_params->shrink_first));
    /*these_params->shrink_first -= 1;*/

    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%d",&(these_params->shrink_every));  
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%d",&(these_params->shrink_last));  
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%f",&(these_params->sigma));   
    current_char = fgetc(pFile);
    while (current_char != '=') {
      current_char = fgetc(pFile);
    }
    fscanf(pFile,"%f",&(these_params->threshold));  
    printf("Do shrinkwrap every %u iterations starting with iteration %u \n using a sigma of %f and threshold of %f.\n Finish shrinking with iteration %u.\n", 
	   these_params->shrink_every, (these_params->shrink_first), 
	   these_params->sigma, these_params->threshold, 
	   these_params->shrink_last);
  } else {
    these_params->shrink_first = these_params->iterate_count;
    these_params->shrink_every = 0;
    these_params->shrink_last = these_params->iterate_count;
printf("Don't perform shrinkwrap.\n");
  }

  /* Make sure the value make sense */
  if ((these_params->shrink_first >= these_params->iterate_count) || 
      (these_params->shrink_last < these_params->shrink_first)){
    these_params->shrink_first = these_params->iterate_count;
    these_params->shrink_every = 0;
    these_params->shrink_last = these_params->iterate_count;
  }


  fclose (pFile);
  return(0);
}

/*------------------------------------------------------------*/
int dm_recon_get_adm(char *filename, 
		     allocatedArrays *ptr_allocated_arrays, 
		     dm_array_real_struct *ptr_mag_error_ras,
		     dm_array_real_struct *ptr_magnitude_ras,
		     reconParams *ptr_recon_params,
		     char *error_string,
		     int my_rank,
		     int p)
{
 
  hid_t h5_file_id;
  int nx, ny, nz, n_dims;
  dm_array_index_t npix,ipix;
  int error_is_present = 0;
  dm_adi_struct this_adi_struct;
  dm_array_real_struct temp_adi_ras;
  dm_array_real_struct temp_mag_error_ras;
  

  if (dm_h5_openwrite(filename,&h5_file_id,error_string,my_rank) 
      != DM_FILEIO_SUCCESS) {
    return(-1);
  }
  printf("Opening file \"%s\".\n",filename);
  
  /* Get the existing adi_struct and allocate memory */
  if (dm_h5_read_adi_info(h5_file_id,&nx,&ny,&nz,&error_is_present,
			  &this_adi_struct,
			  error_string,my_rank) == DM_FILEIO_FAILURE) {
    dm_h5_close(h5_file_id,my_rank);
    return(-1);
  } 
 
  n_dims = 0;
  if (nx > 1) n_dims = 1;
  if (ny > 1) n_dims = 2;
  if (nz > 1) n_dims = 3;

  /* no MPI support for 2D or 1D arrays! */
#ifdef __MPI__
  /*if (n_dims < 3) {
    sprintf(error_string,
	    "dm_recon: Currently only MPI support for 3D arrays!\n");
    return(-1);
    }*/
#endif 

  npix = nx*ny*nz;
  /* check if over available memory for 3D only and disable saving iterates*/
  if (n_dims == 3) {
      if (ptr_recon_params->avg_save_iterates) {
          printf("[WARNING] DM_RECON: saving iterates disabled for 3D!\n");
          ptr_recon_params->avg_save_iterates = 0;
      }
      
#if defined DM_ARRAY_DOUBLE
      /*if ((ptr_recon_params->avg_every != 0) && (npix > 256*256*256)) {
      sprintf(error_string,
	      "Cannot use averaging with arrays larger than 256^3. \n");
      dm_h5_close(h5_file_id,my_rank);
      return(-1);
      }*/
    if (npix > 512*512*512) {
      sprintf(error_string,"With DM array cannot be larger than 512^3. \n");
      dm_h5_close(h5_file_id,my_rank);
      return(-1);
    }
#else
    if ((ptr_recon_params->avg_every != 0) && (npix > 512*512*512)) {
      sprintf(error_string,
	      "Cannot use averaging with arrays larger than 512^3. \n");
      dm_h5_close(h5_file_id,my_rank);
      return(-1);
    }
    if (npix > 1024*1024*1024) {
      sprintf(error_string,"With DM array cannot be larger than 1024^3. \n");
      dm_h5_close(h5_file_id,my_rank);
      return(-1);
    }
#endif /* DM_ARRAY_DOUBLE */
  }

  temp_adi_ras.nx = nx;
  temp_adi_ras.ny = ny;
  temp_adi_ras.nz = nz;
  temp_adi_ras.npix = (dm_array_index_t)nx*
    (dm_array_index_t)ny*(dm_array_index_t)nz;

  /* Allocate memory for temp adi array */
  DM_ARRAY_REAL_STRUCT_INIT((&temp_adi_ras),
			    temp_adi_ras.npix,p);
  

  if (error_is_present == 1) {
    temp_mag_error_ras.nx = temp_adi_ras.nx;
    temp_mag_error_ras.ny = temp_adi_ras.ny;
    temp_mag_error_ras.nz = temp_adi_ras.nz;
    temp_mag_error_ras.npix = temp_adi_ras.npix;
    DM_ARRAY_REAL_STRUCT_INIT((&temp_mag_error_ras),
			      temp_mag_error_ras.npix,p);  
  } else {
    temp_mag_error_ras.nx = 0;
    temp_mag_error_ras.ny = 0;
    temp_mag_error_ras.nz = 0;
    temp_mag_error_ras.npix = 0;
  }
  
  /* Read the adi */
  if (dm_h5_read_adi(h5_file_id,&temp_adi_ras,
		     &temp_mag_error_ras,
		     error_string,my_rank,p) != DM_FILEIO_SUCCESS) {
    dm_h5_close(h5_file_id,my_rank);
    return(-1);
  }
  
  /* Crop the adi if it is 2D (assume it is cropped in 3D already) 
  if (n_dims == 2) {
    dm_recon_crop_2d_adi(&this_adi_struct, &temp_adi_ras,
			 &temp_mag_error_ras,my_rank,p);
			 } */

  /* Set the nx,ny, nz tags in recon_params */
  ptr_recon_params->nx = temp_adi_ras.nx;
  ptr_recon_params->ny = temp_adi_ras.ny;
  ptr_recon_params->nz = temp_adi_ras.nz;
  
  /* Now convert to magnitudes */
  ptr_magnitude_ras->nx = temp_adi_ras.nx;
  ptr_magnitude_ras->ny = temp_adi_ras.ny;
  ptr_magnitude_ras->nz = temp_adi_ras.nz;
  ptr_magnitude_ras->npix = temp_adi_ras.npix;
  
  DM_ARRAY_REAL_STRUCT_INIT(ptr_magnitude_ras, 
			    ptr_magnitude_ras->npix, p);
  ptr_allocated_arrays->magnitude_allocated = 1;
    
  dm_array_magnitude_real(ptr_magnitude_ras, 
			  &temp_adi_ras,1);

  if (error_is_present == 1) {
    dm_array_reassign_real(ptr_mag_error_ras,
			   &temp_mag_error_ras,
			   ptr_allocated_arrays->adi_error_allocated);    
    ptr_allocated_arrays->adi_error_allocated = 1;
    
    /* Convert to absolute error (Factor 0.5 because it is the absolute error 
     * of the magnitudes not intensities)!! */
    for (ipix=0;ipix<ptr_mag_error_ras->local_npix;ipix++) {
      *(ptr_mag_error_ras->real_array+ipix) *= 
	0.5*(*(ptr_magnitude_ras->real_array+ipix));
    }

  } else {
    ptr_mag_error_ras->nx = 0;
    ptr_mag_error_ras->ny = 0;
    ptr_mag_error_ras->nz = 0;
    ptr_mag_error_ras->npix = 0;
  }
  
  /* Now free the temp_adi_array */
     free(temp_adi_ras.real_array);
  
  dm_h5_close(h5_file_id, my_rank);
  return(0);
}

/*------------------------------------------------------------*/
void dm_recon_crop_2d_adi(dm_adi_struct *ptr_adi_struct,
			  dm_array_real_struct *ptr_adi_ras,
			  dm_array_real_struct *ptr_mag_error_ras,
			  int my_rank,
			  int p)
{
  /* crop 2d adi */
  dm_array_crop_2d_real(ptr_adi_ras,ptr_adi_struct->xcenter_offset_pixels,
			ptr_adi_struct->ycenter_offset_pixels,my_rank,p);
  
  /* Do we have an error array? */
  if (ptr_mag_error_ras->nx > 0) {
    dm_array_crop_2d_real(ptr_mag_error_ras,
			  ptr_adi_struct->xcenter_offset_pixels,
			  ptr_adi_struct->ycenter_offset_pixels,
			  my_rank,p);
  } 
}

/*------------------------------------------------------------*/
int dm_recon_get_spt(char *filename,
		     reconParams *ptr_recon_params,
		     dm_spt_struct *ptr_spt_struct,
		     dm_array_byte_struct *ptr_spt_bas,
		     allocatedArrays *ptr_allocated_arrays,
		     char *error_string, int my_rank, int p)
{
  hid_t h5_id;
  int nx,ny,nz,iz,iy,ix;
  int local_nx,local_ny,local_nz;
  int global_ix,global_iy,global_iz;
  int xmax,ymax,zmax,i;
  dm_adi_struct this_adi_struct;

  local_nx = ptr_recon_params->nx;
  local_ny = ptr_recon_params->ny;
  local_nz = ptr_recon_params->nz;

  if (dm_h5_openwrite(filename,&h5_id,error_string,my_rank) 
      != DM_FILEIO_SUCCESS) {
    return(-1);
  }

  if (dm_h5_spt_group_exists(h5_id,my_rank) == 1) {
    printf("Using stored support array.\n");
 
    if (dm_h5_read_spt_info(h5_id,&nx,&ny,&nz,
			    ptr_spt_struct,
			    error_string,my_rank) == DM_FILEIO_FAILURE) {
      return(-1);
    }

    if ((local_nx != nx) || (local_ny != ny) || (local_nz != nz)) {
      sprintf(error_string,
     "dm_recon_get_spt: SPT array should be [%d,%d,%d] instead of [%d.%d,%d]\n",
	      local_nx,local_ny,local_nz,nx,ny,nz);
      return(-1);
    }

    ptr_spt_bas->nx = local_nx;
    ptr_spt_bas->ny = local_ny;
    ptr_spt_bas->nz = local_nz;
    ptr_spt_bas->npix = 
      (dm_array_index_t)local_nx*(dm_array_index_t)local_ny*
      (dm_array_index_t)local_nz;
    DM_ARRAY_BYTE_STRUCT_INIT(ptr_spt_bas,
			      ptr_spt_bas->npix,p);
    
    ptr_allocated_arrays->spt_allocated = 1;

    
    if (dm_h5_read_spt(h5_id,ptr_spt_bas,
		       error_string,my_rank,p) == DM_FILEIO_FAILURE) {
      dm_h5_close(h5_id,my_rank);
      return(-1);
    }
 
    /* ------------In case there is no spt_array -----------*/   
  } else {
    printf("Creating a support array.\n");
    
    ptr_spt_bas->nx = local_nx;
    ptr_spt_bas->ny = local_ny;
    ptr_spt_bas->nz = local_nz;
    ptr_spt_bas->npix = 
      (dm_array_index_t)local_nx*(dm_array_index_t)local_ny*
      (dm_array_index_t)local_nz;
    DM_ARRAY_BYTE_STRUCT_INIT(ptr_spt_bas,
			      ptr_spt_bas->npix,p);
    ptr_allocated_arrays->spt_allocated = 1;

    ptr_spt_struct->support_scaling = 1.0;

    /* We need this in case we use MPI */
    if (local_ny == 1) {
      xmax = local_nx/p;
      ymax = local_ny;
      zmax = local_nz;
    } else if (local_nz == 1) {
      xmax = local_nx;
      ymax = local_ny/p;
      zmax = local_nz;
    } else {
      xmax = local_nx;
      ymax = local_ny;
      zmax = local_nz/p;
    }
    /* Create support to be a cube half the size of array */
    if (local_ny == 1) {
        for (ix=0; ix<local_nx/p; ix++) {
            i = ix;
            global_ix = my_rank*local_nx/p + ix;
            if ((global_ix >= local_nx/4) && (global_ix < 3*local_nx/4)) {
                *(ptr_spt_bas->byte_array+i) = (u_int8_t)1;
            } else {
                *(ptr_spt_bas->byte_array+i) = (u_int8_t)0;
            }
        }
    } else if (local_nz == 1) {
        for (iy=0; iy<local_ny/p; iy++) {
            global_iy = my_rank*local_ny/p + iy;
            for (ix=0; ix<local_nx; ix++) {
                i = (ix+iy*local_nx);
                
                if ((global_iy >= local_ny/4) && (global_iy < 3*local_ny/4)
                    && (ix >= local_nx/4) && (ix < 3*local_nx/4)) {
                    *(ptr_spt_bas->byte_array+i) = (u_int8_t)1;
                } else {
                    *(ptr_spt_bas->byte_array+i) = (u_int8_t)0;
                }
            }
        }
    } else {
        for (iz=0; iz<local_nz/p; iz++) {
            global_iz = my_rank*local_nz/p + iz;

            for (iy=0; iy<local_ny; iy++) {
                for (ix=0; ix<local_nx; ix++) {
                    i = (ix+iy*local_nx+iz*local_nx*local_ny);

                    if ((global_iz >= local_nz/4) && (global_iz < 3*local_nz/4)
                        && (iy >= local_ny/4) && (iy < 3*local_ny/4)
                        && (ix >= local_nx/4) && (ix < 3*local_nx/4)) {
                        *(ptr_spt_bas->byte_array+i) = (u_int8_t)1;
                    } else {
                        *(ptr_spt_bas->byte_array+i) = (u_int8_t)0;
                    }
                }
            }
        }
    }
  }
  
  dm_h5_close(h5_id,my_rank);
  return(0);
}

/*------------------------------------------------------------*/
void dm_recon_init_complex_arrays(dm_array_complex_struct *ptr_current_itn_cas, 
				  dm_array_complex_struct *ptr_final_itn_cas,
				  dm_array_complex_struct *ptr_opt_itn_cas,
				  dm_array_complex_struct *ptr_avg_final_itn_cas,
				  allocatedArrays *ptr_allocated_arrays,
				  reconParams *ptr_recon_params,
				  char *error_string, 
				  int p)
{
  ptr_current_itn_cas->nx = ptr_recon_params->nx;   
  ptr_current_itn_cas->ny = ptr_recon_params->ny;
  ptr_current_itn_cas->nz = ptr_recon_params->nz; 
  ptr_current_itn_cas->npix = (dm_array_index_t)ptr_recon_params->nx*
    (dm_array_index_t)ptr_recon_params->ny*
    (dm_array_index_t)ptr_recon_params->nz;

  ptr_final_itn_cas->nx = ptr_recon_params->nx;
  ptr_final_itn_cas->ny = ptr_recon_params->ny;
  ptr_final_itn_cas->nz = ptr_recon_params->nz;
  ptr_final_itn_cas->npix = (dm_array_index_t)ptr_recon_params->nx*
    (dm_array_index_t)ptr_recon_params->ny*
    (dm_array_index_t)ptr_recon_params->nz;


  DM_ARRAY_COMPLEX_STRUCT_INIT(ptr_current_itn_cas,
			       ptr_current_itn_cas->npix,p);
  DM_ARRAY_COMPLEX_STRUCT_INIT(ptr_final_itn_cas,
			       ptr_final_itn_cas->npix,p);
  ptr_allocated_arrays->current_itn_allocated = 1;
  ptr_allocated_arrays->final_itn_allocated = 1;

  /* Allocate if we need it */
  if ((ptr_recon_params->diffmap) || (ptr_recon_params->shrink_every != 0)) {
    ptr_opt_itn_cas->nx = ptr_recon_params->nx;
    ptr_opt_itn_cas->ny = ptr_recon_params->ny;
    ptr_opt_itn_cas->nz = ptr_recon_params->nz;
    ptr_opt_itn_cas->npix = (dm_array_index_t)ptr_recon_params->nx*
      (dm_array_index_t)ptr_recon_params->ny*
      (dm_array_index_t)ptr_recon_params->nz;
    
    DM_ARRAY_COMPLEX_STRUCT_INIT(ptr_opt_itn_cas,
				 ptr_opt_itn_cas->npix,p);
    ptr_allocated_arrays->opt_itn_allocated = 1;
  }

  /* Allocate if we need it */
  if (ptr_recon_params->avg_every != 0) {
    ptr_avg_final_itn_cas->nx = ptr_recon_params->nx;
    ptr_avg_final_itn_cas->ny = ptr_recon_params->ny;
    ptr_avg_final_itn_cas->nz = ptr_recon_params->nz;
    ptr_avg_final_itn_cas->npix = (dm_array_index_t)ptr_recon_params->nx*
      (dm_array_index_t)ptr_recon_params->ny*
      (dm_array_index_t)ptr_recon_params->nz;

    DM_ARRAY_COMPLEX_STRUCT_INIT(ptr_avg_final_itn_cas,
				 ptr_avg_final_itn_cas->npix,p);
    /* initialize the array to zeros */
    dm_array_zero_complex(ptr_avg_final_itn_cas);
    ptr_allocated_arrays->avg_final_itn_allocated = 1;
  }
}

/*------------------------------------------------------------*/
void dm_recon_init_itn(char *filename,
                       char *error_string,
                       reconParams *ptr_recon_params,
                       dm_array_complex_struct *ptr_current_itn_cas, 
		       dm_array_byte_struct *ptr_spt_bas, 
		       dm_array_real_struct *ptr_magnitude_ras,
		       int p, int my_rank)
{

  dm_array_real power_mags;
  dm_array_real this_power;
  hid_t h5_id;
  int nx,ny,nz, recon_errors_npix;
  int recon_errors_allocated, got_it;
  dm_itn_struct my_itn_struct;
  dm_array_real_struct recon_errors;
  

  if (ptr_recon_params->use_saved_iterate == 1) {
      if (dm_h5_openwrite(filename,&h5_id,error_string,my_rank) 
          != DM_FILEIO_SUCCESS) {
          return;
      }
      
      if (dm_h5_itn_group_exists(h5_id,my_rank) == 1) {
          
          if (dm_h5_read_itn_info(h5_id,&nx,&ny,&nz,
                                  &recon_errors_npix, &my_itn_struct,
                                  error_string,my_rank) == DM_FILEIO_FAILURE) {
              printf("Error reading \"/itn\" group info\n");
              dm_h5_close(h5_id,my_rank);
              exit(1);
          }

          /* Need to allocate errors to be able to read itn group */
          if (my_rank == 0) {
              if (recon_errors_npix > 0) {
                  recon_errors.nx = recon_errors_npix;
                  recon_errors.ny = 1;
                  recon_errors.nz = 1;
                  recon_errors.npix = recon_errors.nx;
                  /* Recon errors will not be split up over the nodes 
                   * so p = 1 in this case.
                   */
                  recon_errors.real_array = 
                      (dm_array_real *)malloc(recon_errors.npix*
                                              sizeof(dm_array_real));
                  recon_errors_allocated = 1;
              } else {
                  recon_errors_allocated = 0;
              }
          }
          
          if (dm_h5_read_itn(h5_id,&recon_errors,ptr_current_itn_cas,
                             error_string,my_rank,p) == DM_FILEIO_FAILURE) {
              printf("Error reading \"/itn\" array\n");
              dm_h5_close(h5_id,my_rank);
              exit(1);
          }

          /* Close the file */
          dm_h5_close(h5_id,my_rank);
          
          /* We don't need the errors */
          if (recon_errors_allocated == 1)
              free(recon_errors.real_array);

          /* We have an initial guess */
          got_it = 1;
      } else {
          /* There was no saved itn */
          got_it = 0;
      }
  }

  /* Either the user did not want to or there was no saved itn */
  if ((ptr_recon_params->use_saved_iterate == 0) || got_it == 0) {
      
      dm_array_zero_complex(ptr_current_itn_cas);
      
      /* fill with random numbers in area of support */
      dm_array_rand(ptr_current_itn_cas,0);
      dm_array_multiply_complex_byte(ptr_current_itn_cas,ptr_spt_bas);   
      
      
      power_mags = dm_array_total_power_real(ptr_magnitude_ras,0);
      this_power = dm_array_total_power_complex(ptr_current_itn_cas,NULL,0);
      
      dm_array_multiply_real_scalar(ptr_current_itn_cas,
                                    sqrt(power_mags/this_power));
  }

}

/*------------------------------------------------------------*/
void dm_recon_do_diffmap(reconParams *ptr_recon_params,
			 dm_array_real_struct *ptr_magnitude_ras,
			 dm_array_real_struct *ptr_mag_error_ras,
			 dm_array_complex_struct *ptr_current_itn_cas,
			 dm_array_complex_struct *ptr_final_itn_cas,
			 dm_array_complex_struct *ptr_opt_itn_cas,
			 dm_array_complex_struct *ptr_avg_final_itn_cas,
			 dm_array_byte_struct *ptr_spt_bas,
			 dm_array_real_struct *ptr_recon_errors,
			 int p, int my_rank, int it, int *ptr_avg_count,
                         char *filename)

{
  dm_array_real total_power;

  /* NOTE:
   *  - Support estimate = ptr_opt_itn_cas
   *  - Modulus estimate = ptr_final_itn_cas, this will be saved to 
   *    the file or used for averaging!
   */

  /* support estimate */
  dm_array_copy_complex(ptr_opt_itn_cas, ptr_current_itn_cas);

  dm_array_fft(ptr_opt_itn_cas,p,DM_ARRAY_FORWARD_FFT,
	       my_rank);

  /* Pass with error array if desired */
  if ((ptr_recon_params->use_error_array) && (ptr_mag_error_ras->npix != 0)) {
    dm_array_transfer_magnitudes(ptr_opt_itn_cas,ptr_magnitude_ras,
				 ptr_mag_error_ras,0);
  } else {
    dm_array_transfer_magnitudes(ptr_opt_itn_cas,ptr_magnitude_ras,NULL,0);
  }

  dm_array_fft(ptr_opt_itn_cas,p,DM_ARRAY_INVERSE_FFT,
               my_rank);


  dm_array_multiply_real_scalar(ptr_opt_itn_cas,
				(dm_array_real)(ptr_recon_params->gamma_mod+1));
  
  dm_array_multiply_real_scalar(ptr_current_itn_cas,
				(dm_array_real)(ptr_recon_params->gamma_mod));

  dm_array_subtract_complex(ptr_opt_itn_cas,ptr_current_itn_cas);

  dm_array_multiply_complex_byte(ptr_opt_itn_cas,ptr_spt_bas);


  if (ptr_recon_params->positivity == 1) {
      dm_positivity(ptr_opt_itn_cas, 0,
                    ptr_recon_params->imaginary_only);
  }

  dm_array_multiply_real_scalar(ptr_opt_itn_cas, 
				(dm_array_real)ptr_recon_params->beta);

  /* undo modification to current_itn_cas */
  dm_array_multiply_real_scalar(ptr_current_itn_cas,
				(dm_array_real)(1/ptr_recon_params->gamma_mod));

  /* modulus estimate */
  dm_array_copy_complex(ptr_final_itn_cas, ptr_current_itn_cas);

  dm_array_multiply_complex_byte(ptr_final_itn_cas,ptr_spt_bas);
  
  if (ptr_recon_params->positivity == 1) {
      dm_positivity(ptr_final_itn_cas, 0,
                    ptr_recon_params->imaginary_only);
  }

  dm_array_multiply_real_scalar(ptr_final_itn_cas,
				(dm_array_real)(ptr_recon_params->gamma_spt+1));

  dm_array_multiply_real_scalar(ptr_current_itn_cas,
				(dm_array_real)(ptr_recon_params->gamma_spt));

  dm_array_subtract_complex(ptr_final_itn_cas,ptr_current_itn_cas);

  dm_array_fft(ptr_final_itn_cas,p,DM_ARRAY_FORWARD_FFT,
	       my_rank);

  /* Pass with error array if desired */
  if ((ptr_recon_params->use_error_array) && (ptr_mag_error_ras->npix != 0)) {
    dm_array_transfer_magnitudes(ptr_final_itn_cas,ptr_magnitude_ras,
				 ptr_mag_error_ras,0);
  } else {
    dm_array_transfer_magnitudes(ptr_final_itn_cas,ptr_magnitude_ras,NULL,0);
  }
  
  dm_array_fft(ptr_final_itn_cas,p,DM_ARRAY_INVERSE_FFT,
	       my_rank);
  
  dm_array_multiply_real_scalar(ptr_final_itn_cas, 
				(dm_array_real)ptr_recon_params->beta);

  /* undo modification to my_current_itn_cas */                                 
  dm_array_multiply_real_scalar(ptr_current_itn_cas,
				(dm_array_real)(1/ptr_recon_params->gamma_spt));
  
  /* Calculate new object */
  dm_array_add_complex(ptr_current_itn_cas,ptr_final_itn_cas);

  dm_array_subtract_complex(ptr_current_itn_cas,ptr_opt_itn_cas);

  /* Undo modification on Fourier estimate  */
  dm_array_multiply_real_scalar(ptr_final_itn_cas,
				(dm_array_real)(1./ptr_recon_params->beta));

  /* Calculate error */
  dm_recon_error(ptr_final_itn_cas,ptr_spt_bas,ptr_recon_errors,my_rank,it);

  /* Average here if necessary */
  if (((it+1) > ptr_recon_params->avg_after) && 
      ((it+1) % ptr_recon_params->avg_every == 0)) {
    dm_recon_average(ptr_avg_final_itn_cas,ptr_final_itn_cas,
                     *ptr_avg_count,ptr_recon_params->avg_save_iterates,
                     it,my_rank,p,filename);
    (*ptr_avg_count)++;
    printf("Calculating average with iteration %i. \n", it);
  }

  /* Shrinkwrap here if necessary */
  if ((it >= ptr_recon_params->shrink_first) && 
      (it <= ptr_recon_params->shrink_last) &&
      ((it-ptr_recon_params->shrink_first)
       % ptr_recon_params->shrink_every == 0)) {      
      dm_recon_shrinkwrap(ptr_recon_params,ptr_final_itn_cas, 
                          ptr_opt_itn_cas,ptr_spt_bas,my_rank,p);
      printf("Performed SHRINKWRAP after iteration %i.\n", it);
  }
}

/*------------------------------------------------------------*/
void dm_recon_do_hio(reconParams *ptr_recon_params,
		     dm_array_real_struct *ptr_magnitude_ras,
		     dm_array_real_struct *ptr_mag_error_ras,
		     dm_array_complex_struct *ptr_current_itn_cas,
		     dm_array_complex_struct *ptr_final_itn_cas,
		     dm_array_complex_struct *ptr_opt_itn_cas,
		     dm_array_complex_struct *ptr_avg_final_itn_cas,
		     dm_array_byte_struct *ptr_spt_bas,
		     dm_array_real_struct *ptr_recon_errors,
		     int p, int my_rank, int it, int *ptr_avg_count,
                     char *filename)
  
{
  dm_array_real total_power;
  dm_array_index_t ipix;

  /* NOTE: the array that will be averaged and saved is 
   *     ptr_final_itn_cas.
   */
  dm_array_copy_complex(ptr_final_itn_cas, ptr_current_itn_cas);

  dm_array_fft(ptr_final_itn_cas,p,DM_ARRAY_FORWARD_FFT,
	       my_rank);

  /* Pass with error array if desired */
  if ((ptr_recon_params->use_error_array) && (ptr_mag_error_ras->npix != 0)) {
    dm_array_transfer_magnitudes(ptr_final_itn_cas,ptr_magnitude_ras,
				 ptr_mag_error_ras,0);
  } else {
    dm_array_transfer_magnitudes(ptr_final_itn_cas,ptr_magnitude_ras,NULL,0);
  }

  dm_array_fft(ptr_final_itn_cas,p,DM_ARRAY_INVERSE_FFT,
	       my_rank);

  /* Calculate error here! */
  dm_recon_error(ptr_final_itn_cas,ptr_spt_bas,ptr_recon_errors,my_rank,it);

  /* Average here if necessary */
  if (((it+1) > ptr_recon_params->avg_after) && 
      ((it+1) % ptr_recon_params->avg_every == 0)) {
      dm_recon_average(ptr_avg_final_itn_cas,ptr_final_itn_cas,
                       *ptr_avg_count,ptr_recon_params->avg_save_iterates,
                       it,my_rank,p,filename);
      (*ptr_avg_count)++;
    printf("Calculating average with iteration %i. \n", it);
  }
 
  /* Shrinkwrap here if necessary */
  if ((it >= ptr_recon_params->shrink_first) && 
      (it <= ptr_recon_params->shrink_last) &&
      ((it-ptr_recon_params->shrink_first)
       % ptr_recon_params->shrink_every == 0)) {      
      dm_recon_shrinkwrap(ptr_recon_params,ptr_final_itn_cas, 
                          ptr_opt_itn_cas,ptr_spt_bas,my_rank,p);
      printf("Performed SHRINKWRAP at iteration %i.\n", it);
  }

  /* Next iterate - This is not necessary for last iteration */
  if (it < ptr_recon_params->iterate_count -1) {

    dm_array_multiply_real_scalar(ptr_final_itn_cas,
				  (dm_array_real)(ptr_recon_params->beta));
        
    dm_array_subtract_complex(ptr_current_itn_cas,ptr_final_itn_cas);

    dm_array_multiply_complex_byte(ptr_final_itn_cas,ptr_spt_bas);
    
    if (ptr_recon_params->positivity == 1) {
        dm_positivity(ptr_final_itn_cas, 1,
                      ptr_recon_params->imaginary_only);
    }

    /* Need to loop through both arrays to set pixels that are non zero
     * in the first array to zero in the second. 
     */
    for (ipix=0; ipix<(ptr_final_itn_cas->local_npix); ipix++) {
        if ((c_re(ptr_final_itn_cas->complex_array,ipix) != 0.) &&
            (c_im(ptr_final_itn_cas->complex_array,ipix) != 0.)) {
            c_re(ptr_current_itn_cas->complex_array,ipix) = 0.0;
            c_im(ptr_current_itn_cas->complex_array,ipix) = 0.0;
        }
    }
    
    dm_array_multiply_real_scalar(ptr_final_itn_cas, 
				  (dm_array_real)(1/ptr_recon_params->beta));   

    dm_array_add_complex(ptr_current_itn_cas,ptr_final_itn_cas); 
    
  }
}


/*------------------------------------------------------------*/
void dm_recon_error(dm_array_complex_struct *ptr_final_itn_cas,
		    dm_array_byte_struct *ptr_spt_bas,
		    dm_array_real_struct *ptr_recon_errors, 
		    int my_rank, int it)
{
  dm_array_index_t ipix;
  dm_array_real total_power_outside, total_power_inside;

  total_power_inside = 
    dm_array_total_power_complex(ptr_final_itn_cas,ptr_spt_bas,0);

  total_power_outside = 
    dm_array_total_power_complex(ptr_final_itn_cas,ptr_spt_bas,1);
  
  if (my_rank == 0) {
    *(ptr_recon_errors->real_array + it) = 
      total_power_outside/total_power_inside;
    printf("Error for iteration %d: %f\n",it,
	   (float)*(ptr_recon_errors->real_array + it));
  }
}

/*------------------------------------------------------------*/
void dm_recon_average(dm_array_complex_struct *ptr_avg_final_itn_cas,
		      dm_array_complex_struct *ptr_final_itn_cas,
		      int avg_count,
                      int save_iterate,
                      int iteration,
                      int my_rank,
                      int p,
                      char *filename)
{
  dm_array_complex *this_value;
  dm_array_real this_ratio;
  dm_array_real global_phase;
  hid_t h5_file_id;
  char error_string[STRLEN], this_name[STRLEN], base_name[STRLEN];
  dm_itn_struct this_itn_struct;
  dm_array_real_struct recon_errors;
  size_t num;

   DM_ARRAY_COMPLEX_MALLOC(this_value,1);
   /* Adjust the global phase */
  if (avg_count == 0) {
    dm_array_square_sum_complex(ptr_final_itn_cas,NULL,this_value);

    /* This should give us the right quadrant */
    this_ratio = (dm_array_real)
      atan2((double)c_im(this_value,0),(double)c_re(this_value,0));

    /* New angle according to Chapman */
    this_ratio *= (dm_array_real)(-1./2.);
  } else {
    dm_array_square_sum_complex(ptr_final_itn_cas,ptr_avg_final_itn_cas,
				this_value);

    /* This should give us the right quadrant */
    this_ratio = (dm_array_real)
      atan2((double)c_im(this_value,0),(double)c_re(this_value,0));

    this_ratio *= (dm_array_real)-1.;
    
  }
  printf("Adjusting phase by %f radians.\n",this_ratio);

  /* Convert to {Re,Im} preserving magnitude */  
  c_re(this_value,0) = (dm_array_real)cos((double)this_ratio);
  c_im(this_value,0) = (dm_array_real)sin((double)this_ratio);

  /* multiply by new phase */
  dm_array_multiply_complex_scalar(ptr_final_itn_cas,this_value);

  if (save_iterate) {
      /* this is how we can forego recon_errors */
      recon_errors.nx=0;
      recon_errors.ny=0;
      recon_errors.nz=0;
      recon_errors.npix=0;

      /* save iterate_count */
      this_itn_struct.iterate_count = iteration;

      /* create output filename */
      num = strlen(filename)-3;
      strncpy(base_name,filename,num);
      base_name[num]='\0'; /* make sure it's null terminated */
      sprintf(this_name,"%s_%d.h5",base_name,iteration);

      if (dm_h5_create(this_name,&h5_file_id,error_string,my_rank) 
          != DM_FILEIO_SUCCESS) {
          printf(error_string);
      } else {
      
          /* Save final iterate in object space */
          if (dm_h5_write_itn(h5_file_id,&this_itn_struct, ptr_final_itn_cas,
                              &recon_errors, error_string,my_rank,p) 
              != DM_FILEIO_SUCCESS) {
              dm_h5_close(h5_file_id,my_rank);
              printf(error_string);
          }

          dm_h5_close(h5_file_id,my_rank);
      }
  }
  
  /* Add the new array to the old ones */
  dm_array_add_complex(ptr_avg_final_itn_cas, ptr_final_itn_cas);

  /* Undo phase adjustment */
  c_re(this_value,0) = (dm_array_real)cos((double)(-this_ratio));
  c_im(this_value,0) = (dm_array_real)sin((double)(-this_ratio));
  dm_array_multiply_complex_scalar(ptr_final_itn_cas,this_value);
  
  DM_ARRAY_COMPLEX_FREE(this_value);
}


/*------------------------------------------------------------*/
void dm_modify_spt(dm_array_real_struct *ptr_ras,
		   dm_array_byte_struct *ptr_bas,
		   float threshold)
{ 
  /* modifies byte array such that 
     if ras element is greater than threshold set that bas element to 1, 
     otherwise set it to 0                            
   */
  dm_array_index_t ipix;
  
  if (ptr_ras->npix != ptr_bas->npix) return;
  
  for (ipix=0; ipix<(ptr_ras->local_npix); ipix++) {
    if (*(ptr_ras->real_array+ipix) > threshold) {
      *(ptr_bas->byte_array+ipix) = 1;
    } else {
      *(ptr_bas->byte_array+ipix) = 0;
    }
  }
  
#if USE_MPI 
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE MPI */
}

/*------------------------------------------------------------*/
void dm_positivity(dm_array_complex_struct *ptr_cas,
                   int using_hio,
		   int imaginary_only)
{
  /* 
     Applies positivity constraint to both real and imaginary parts
     unless imaginary_only = 1
  */
  dm_array_index_t ipix;

  /* We have to distinguish between HIO and others, others first. */
  if (using_hio == 0) {
      for (ipix=0; ipix<(ptr_cas->local_npix); ipix++) {
          if ( imaginary_only != 1) {
              if (c_re(ptr_cas->complex_array,ipix) < 0.) {
                  c_re(ptr_cas->complex_array,ipix) = 0.0;
              }
          }
          if (c_im(ptr_cas->complex_array,ipix) < 0.) {
              c_im(ptr_cas->complex_array,ipix) = 0.0;
          }
      }
  } else {
      /* If we are using HIO, then set the pixels violating the
       * positivity constraint in the first array to 0 (both re and im),
       */
      for (ipix=0; ipix<(ptr_cas->local_npix); ipix++) {
          if ( imaginary_only != 1) {
              if ((c_re(ptr_cas->complex_array,ipix) < 0.) ||
                  (c_im(ptr_cas->complex_array,ipix) < 0.)) {
                  c_re(ptr_cas->complex_array,ipix) = 0.0;
                  c_im(ptr_cas->complex_array,ipix) = 0.0;
              }
          } else {
              if (c_im(ptr_cas->complex_array,ipix) < 0.) {
                  c_re(ptr_cas->complex_array,ipix) = 0.0;
                  c_im(ptr_cas->complex_array,ipix) = 0.0;
              }
          }
      }
  }
   
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE MPI */
}

/*------------------------------------------------------------*/
int dm_recon_shrinkwrap(reconParams *ptr_recon_params,
			dm_array_complex_struct *ptr_final_itn_cas,
			dm_array_complex_struct *ptr_opt_itn_cas,
			dm_array_byte_struct *ptr_spt_bas,
			int my_rank,
			int p)
{
  float sigma,threshold;
  dm_array_real_struct real_array_struct;
  dm_array_real max_magnitude, inverse_norm;

  sigma = ptr_recon_params->sigma;
  threshold = ptr_recon_params->threshold;

  printf("Sigma is %f. \n", ptr_recon_params->sigma);
 
  /* Allocate memory for real part */
  real_array_struct.nx = ptr_recon_params->nx;
  real_array_struct.ny = ptr_recon_params->ny;
  real_array_struct.nz = ptr_recon_params->nz;
  real_array_struct.npix = (dm_array_index_t)real_array_struct.nx*
    (dm_array_index_t)real_array_struct.ny*
    (dm_array_index_t)real_array_struct.nz;
  DM_ARRAY_REAL_STRUCT_INIT((&real_array_struct),real_array_struct.npix,p);

  /* store guassian in opt_itn_cas_array_struct to save memory. Calculate 
   * Gaussian in inverse space. Need to take unit conversion into account:
   * delta_f = f_max/nx/2 with f_max=1/(2*delta_x) 
   *    ==> delta_f = 1/(nx*delta_x)
   */
  sigma = ptr_recon_params->nx/(2*PI*sigma);

  dm_array_load_gaussian(ptr_opt_itn_cas,sigma,sigma,sigma, 
			 0,1, p, my_rank);

  /* No need to FFT Gaussian since we initialized it in inverse space. */
  dm_array_fft(ptr_final_itn_cas,p,DM_ARRAY_FORWARD_FFT,
	       my_rank);

  dm_array_multiply_complex(ptr_opt_itn_cas, ptr_final_itn_cas);

  dm_array_fft(ptr_opt_itn_cas,p,DM_ARRAY_INVERSE_FFT,
	       my_rank);

  /* Bring current iterate back to object space */
  dm_array_fft(ptr_final_itn_cas,p,DM_ARRAY_INVERSE_FFT,
	       my_rank);


  /* Look at magnitudes */
  dm_array_magnitude_complex(&real_array_struct, ptr_opt_itn_cas);
 
  /* Look at magnitudes */
  dm_array_magnitude_complex(&real_array_struct, ptr_opt_itn_cas);

  max_magnitude = dm_array_max_real(&real_array_struct, p);

  /* Absolute threshold */
  max_magnitude = max_magnitude*threshold;

  dm_modify_spt(&real_array_struct, ptr_spt_bas, max_magnitude);

  /* Free the array for the real part */
  free(real_array_struct.real_array);

  /* modify sigma each time using shrink_count */
  if (sigma > 1.5) {
    /* reduce sigma by 1% */
    sigma = ptr_recon_params->sigma;
    ptr_recon_params->sigma = sigma - (sigma/100);  
    printf("Sigma reduced to %f. \n", ptr_recon_params->sigma);
  }
 
}


/*------------------------------------------------------------*/
int dm_recon_save_recon(char *filename, 
			reconParams *ptr_recon_params,
			dm_itn_struct *ptr_itn_struct, 
			dm_array_complex_struct *ptr_final_itn_cas, 
			dm_spt_struct *ptr_spt_struct,
			dm_array_byte_struct *ptr_spt_bas,
			dm_array_real_struct *ptr_recon_errors,
			char *error_string,
			int my_rank, 
			int p)
{
  hid_t h5_file_id;
  dm_comment_struct this_comment_struct;
  char this_string[STRLEN];
  int free_all_comment_struct;
  time_t rawtime;

  if (dm_h5_openwrite(filename,&h5_file_id,error_string,my_rank) 
      != DM_FILEIO_SUCCESS) {
    return(-1);
  }

  /* Save final iterate in object space */
  if (dm_h5_write_itn(h5_file_id,ptr_itn_struct, ptr_final_itn_cas,
		      ptr_recon_errors, error_string,my_rank,p) 
      != DM_FILEIO_SUCCESS) {
    dm_h5_close(h5_file_id,my_rank);
    return(-1);
  }

  /* Save the current support */
  if (dm_h5_write_spt(h5_file_id,ptr_spt_struct, ptr_spt_bas,
		      error_string,my_rank,p) != DM_FILEIO_SUCCESS) {
    dm_h5_close(h5_file_id,my_rank);
    return(-1);
  } 


  /* Currently up to 4 strings possible */
  this_comment_struct.string_array = 
    (char *)malloc(sizeof(char)*4*STRLEN);
  this_comment_struct.n_strings_max = 5;
  this_comment_struct.n_strings = 0;
  this_comment_struct.string_length = STRLEN;
  
  time(&rawtime);
  sprintf(this_string,"DM_RECON, %s",ctime(&rawtime));
  
  dm_add_string_to_comments(this_string, &this_comment_struct);

  if (ptr_recon_params->hio) {
    sprintf(this_string,"   Performed %d iterations of HIO.\n",
	    ptr_recon_params->iterate_count);
  } else if (ptr_recon_params->diffmap) {
    sprintf(this_string,"   Performed %d iterations of DM.\n",
	    ptr_recon_params->iterate_count);
  }

  dm_add_string_to_comments(this_string, &this_comment_struct);

  /* Averaging? */
  if (ptr_recon_params->avg_after < ptr_itn_struct->iterate_count){
    sprintf(this_string,"   Averaged every %d after %d iterations\n", 
	    ptr_recon_params->avg_every, ptr_recon_params->avg_after);
    dm_add_string_to_comments(this_string,&this_comment_struct);
  }

  /* Shrinkwrapping? */
  if (ptr_recon_params->shrink_first < ptr_itn_struct->iterate_count) {
    sprintf(this_string,
	   "   Performed shrinkwrap every %d starting with %d up to %d iterations\n",
	   ptr_recon_params->shrink_every, ptr_recon_params->shrink_first,
	   ptr_recon_params->shrink_last);
    dm_add_string_to_comments(this_string, &this_comment_struct);
  }
  
  /* are we adding or creating new? */
  if (dm_h5_comments_group_exists(h5_file_id,my_rank) == 1) {
      free_all_comment_struct = 0;
      if (dm_h5_add_comments(h5_file_id,&this_comment_struct,
                             error_string,my_rank) != DM_FILEIO_SUCCESS) {
          printf("%s\n",error_string);
          dm_h5_close(h5_file_id,my_rank);
          exit(1);
      }
  } else {
    printf("[WARNING] dm_recon: no comments found, creating comments!\n");
    free_all_comment_struct = 1;

    this_comment_struct.specimen_name =
        (char *)malloc(this_comment_struct.string_length);

    this_comment_struct.collection_date =
        (char *)malloc(this_comment_struct.string_length);
    
    dm_add_specimen_name_to_comments("Unkown",&this_comment_struct);
    time(&rawtime);

    dm_add_collection_date_to_comments(ctime(&rawtime),&this_comment_struct);
    
    if (dm_h5_create_comments(h5_file_id,&this_comment_struct,
			      error_string,my_rank) != DM_FILEIO_SUCCESS) {
      printf("%s\n",error_string);
      dm_h5_close(h5_file_id,my_rank);
      exit(1);
    }
  }
 
  free(this_comment_struct.string_array);
  if (free_all_comment_struct) {
      free(this_comment_struct.specimen_name);
      free(this_comment_struct.collection_date);
  }
      

  dm_h5_close(h5_file_id,my_rank);

  printf("Updated file \"%s\"\n",filename);

  return(0);
}

/*------------------------------------------------------------*/
void dm_recon_free_mem(allocatedArrays *ptr_allocated_arrays,
		       dm_array_real_struct *ptr_mag_error_ras,
		       dm_array_real_struct *ptr_magnitude_ras,
		       dm_array_byte_struct *ptr_spt_bas, 
		       dm_array_complex_struct *ptr_current_itn_cas,
		       dm_array_complex_struct *ptr_avg_final_itn_cas, 
		       dm_array_complex_struct *ptr_final_itn_cas,
		       dm_array_complex_struct *ptr_opt_itn_cas,
		       dm_array_real_struct *ptr_recon_errors,
                       int my_rank)
{
  if (ptr_allocated_arrays->adi_error_allocated == 1)
    free(ptr_mag_error_ras->real_array);

  if (ptr_allocated_arrays->magnitude_allocated == 1)
    free(ptr_magnitude_ras->real_array);
  
  if (ptr_allocated_arrays->spt_allocated == 1)
    free(ptr_spt_bas->byte_array);
  
  if (ptr_allocated_arrays->current_itn_allocated == 1)
    DM_ARRAY_COMPLEX_FREE(ptr_current_itn_cas->complex_array);
  
  if (ptr_allocated_arrays->avg_final_itn_allocated == 1)
    DM_ARRAY_COMPLEX_FREE(ptr_avg_final_itn_cas->complex_array);

  if (ptr_allocated_arrays->final_itn_allocated == 1)
    DM_ARRAY_COMPLEX_FREE(ptr_final_itn_cas->complex_array);
  
  if (ptr_allocated_arrays->opt_itn_allocated == 1)
    DM_ARRAY_COMPLEX_FREE(ptr_opt_itn_cas->complex_array);

  if (my_rank == 0) {
      free(ptr_recon_errors->real_array);
  }
}

