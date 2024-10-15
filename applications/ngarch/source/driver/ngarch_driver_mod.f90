!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the ngarch miniapp.
!>
!> @details This runs only a single physics scheme per timestep, allowing for
!>          faster prototyping of changes within a scheme.
!>
module ngarch_driver_mod

  use constants_mod,              only : i_def, str_def, &
                                         r_def, r_second
  use gungho_modeldb_mod,         only : modeldb_type
  use field_collection_mod,       only : field_collection_type
  use field_mod,                  only : field_type
  use field_array_mod,            only : field_array_type
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_INFO,    &
                                         LOG_LEVEL_ERROR,   &
                                         LOG_LEVEL_TRACE
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use gungho_init_fields_mod,     only : create_model_data, initialise_model_data
  use sci_field_bundle_builtins_mod, only : clone_bundle, &
                                            set_bundle_scalar
  use mr_indices_mod,             only : nummr

  use casim_alg_mod, only : casim_alg


  implicit none

  private

  public step

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Executes a single science scheme for a single timestep
  !> @details In this example the CASIM microphysics scheme is used. This is
  !>          expected to be swapped out for another scheme while prototyping
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine step( program_name, modeldb )

    implicit none

    character(*),       intent(in)    :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    type( field_collection_type ), pointer :: collection
    type( field_type ),            pointer :: theta, rho, dtheta_mphys
    type( field_type )                     :: dmr_mphys(nummr), dcfl, dcff, dbcf
    type( field_array_type ),      pointer :: mr

    type(mesh_type),               pointer :: mesh

    collection => modeldb%fields%get_field_collection( "moisture_fields" )
    call collection%get_field( "mr", mr )
    call clone_bundle( mr%bundle, dmr_mphys, nummr )
    call set_bundle_scalar( 0.0_r_def, dmr_mphys, nummr )

    call modeldb%model_data%microphysics_fields%get_field( 'dtheta_mphys', dtheta_mphys )

    collection => modeldb%fields%get_field_collection( "depository" )
    call collection%get_field( "theta", theta )
    call collection%get_field( "rho", rho )
    mesh => theta%get_mesh()

    ! Call an algorithm
    call log_event( program_name//": Running CASIM", LOG_LEVEL_INFO )
    call casim_alg( mr%bundle, theta, rho, modeldb%model_data%derived_fields, &
                    modeldb%model_data%microphysics_fields,                   &
                    modeldb%model_data%cloud_fields,                          &
                    modeldb%model_data%aerosol_fields,                        &
                    modeldb%model_data%turbulence_fields, mesh,               &
                    dmr_mphys, dtheta_mphys, dcfl, dcff, dbcf )
    call log_event( program_name//": CASIM completed", LOG_LEVEL_INFO )

  end subroutine step

end module ngarch_driver_mod
