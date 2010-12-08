/* This is the file dm_recon.h */

#ifndef DM_RECON_H
#define DM_RECON_H

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "../util/dm.h"
#include "../util/dm_fileio.h"
#include "../util/dm_array.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

    typedef struct {
        int hio;
        int diffmap;
        float beta;
        float gamma_mod;
        float gamma_spt;
        int iterate_count;
        int use_saved_iterate;
        int use_error_array;
        int positivity;
        int imaginary_only;
        int avg_after;
        int avg_every;
        int avg_save_iterates;
        int shrink_first;
        int shrink_every;
        int shrink_last;
        float sigma;
        float threshold;
        int nx;
        int ny;
        int nz;
    } reconParams;


  typedef struct {
    int adi_error_allocated;
    int magnitude_allocated;
    int spt_allocated;
    int current_itn_allocated;
    int avg_final_itn_allocated;
    int final_itn_allocated;
    int opt_itn_allocated;
  } allocatedArrays;

  /*------------------------------------------------------------*/
  int dm_recon_read_script(char *scriptname, 
			   char *filename,
			   reconParams *these_params, 
			   char *error_string);
  /*------------------------------------------------------------*/
  int dm_recon_get_adm(char *filename, 
		       allocatedArrays *ptr_allocated_arrays, 
		       dm_array_real_struct *ptr_adi_error_ras,
		       dm_array_real_struct *ptr_magnitude_ras,
		       reconParams *ptr_recon_params,
		       char *error_string,
		       int my_rank,
		       int p);

  /*------------------------------------------------------------*/
  void dm_recon_crop_2d_adi(dm_adi_struct *ptr_adi_struct,
			    dm_array_real_struct *ptr_adi_ras,
			    dm_array_real_struct *ptr_adi_error_ras,
			    int my_rank,int p);

  /*------------------------------------------------------------*/
  int dm_recon_get_spt(char *filename,
		       reconParams *ptr_recon_params,
		       dm_spt_struct *ptr_spt_struct,
		       dm_array_byte_struct *ptr_spt_bas,
		       allocatedArrays *ptr_allocated_arrays,
		       char *error_string, int my_rank, int p);
    
  /*------------------------------------------------------------*/
  void dm_recon_init_complex_arrays(dm_array_complex_struct *ptr_current_itn_cas, 
				    dm_array_complex_struct *ptr_final_itn_cas,
				    dm_array_complex_struct *ptr_opt_itn_cas,
				    dm_array_complex_struct *ptr_avg_final_itn_cas,
				    allocatedArrays *ptr_allocated_arrays,
				    reconParams *ptr_recon_params,
				    char *error_string, 
				    int p);

    /*------------------------------------------------------------*/
    void dm_recon_init_itn(char *filename,
                           char *error_string,
                           reconParams *ptr_recon_params,
                           dm_array_complex_struct *ptr_current_itn_cas, 
                           dm_array_byte_struct *ptr_spt_bas, 
                           dm_array_real_struct *ptr_magnitude_ras,
                           int p, int my_rank);

  /*------------------------------------------------------------*/
  void dm_recon_do_diffmap(reconParams *ptr_recon_params,
			   dm_array_real_struct *ptr_magnitude_ras,
			   dm_array_real_struct *ptr_adi_error_ras,
			   dm_array_complex_struct *ptr_current_itn_cas,
			   dm_array_complex_struct *ptr_final_itn_cas,
			   dm_array_complex_struct *ptr_opt_itn_cas,
			   dm_array_complex_struct *ptr_avg_final_itn_cas,
			   dm_array_byte_struct *ptr_spt_bas,
			   dm_array_real_struct *ptr_recon_errors,
			   int p, int my_rank, int it, int *ptr_avg_count,
                           char *filename);

 /*------------------------------------------------------------*/
  void dm_recon_do_hio(reconParams *ptr_recon_params,
		       dm_array_real_struct *ptr_magnitude_ras,
		       dm_array_real_struct *ptr_adi_error_ras,
		       dm_array_complex_struct *ptr_current_itn_cas,
		       dm_array_complex_struct *ptr_final_itn_cas,
		       dm_array_complex_struct *ptr_opt_itn_cas,
		       dm_array_complex_struct *ptr_avg_final_itn_cas,
		       dm_array_byte_struct *ptr_spt_bas,
		       dm_array_real_struct *ptr_recon_errors,
		       int p,int my_rank, int it, int *ptr_avg_count,
                       char *filename);

  /*------------------------------------------------------------*/
  void dm_recon_error(dm_array_complex_struct *ptr_final_itn_cas,
		      dm_array_byte_struct *ptr_spt_bas,
		      dm_array_real_struct *ptr_recon_errors, 
		      int my_rank, int it);

  /*------------------------------------------------------------*/
  void dm_recon_average(dm_array_complex_struct *ptr_avg_final_itn_cas,
			dm_array_complex_struct *ptr_final_itn_cas,
			int avg_count,
                        int save_iterate,
                        int iteration,
                        int my_rank,
                        int p,
                        char *filename);

  /** This routine modifies the support for Shrinkwrap */
  void dm_modify_spt(dm_array_real_struct *ptr_ras,
		     dm_array_byte_struct *ptr_bas,
		     float threshold);

  /** This routine applies the positivity constraint */
  void dm_positivity(dm_array_complex_struct *ptr_cas,
                     int using_hio,
                     int imaginary_only);

  /*------------------------------------------------------------*/
  int dm_recon_shrinkwrap(reconParams *ptr_recon_params,
			  dm_array_complex_struct *ptr_final_itn_cas,
			  dm_array_complex_struct *ptr_opt_itn_cas,
			  dm_array_byte_struct *ptr_spt_bas,
			  int my_rank,
			  int p);


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
			  int p);

  /*------------------------------------------------------------*/
  void dm_recon_free_mem(allocatedArrays *ptr_allocated_arrays,
			 dm_array_real_struct *ptr_adi_error_ras,
			 dm_array_real_struct *ptr_magnitude_ras,
			 dm_array_byte_struct *ptr_spt_bas, 
			 dm_array_complex_struct *ptr_current_itn_cas,
			 dm_array_complex_struct *ptr_avg_final_itn_cas, 
			 dm_array_complex_struct *ptr_final_itn_cas,
			 dm_array_complex_struct *ptr_opt_itn_cas,
			 dm_array_real_struct *ptr_recon_errors,
                         int my_rank);


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* #ifndef DM_RECON_H */
