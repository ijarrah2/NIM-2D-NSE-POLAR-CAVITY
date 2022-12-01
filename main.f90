!! *********************************************************************************
!! *****  _   _  ___  __  __        ____   ____          _   _  ____   _____    ****
!! ***** | \ | ||_ _||  \/  |      |___ \ |  _ \        | \ | |/ ___| | ____|   ****
!! ***** |  \| | | | | |\/| | _____  __) || | | | _____ |  \| |\___ \ |  _|     ****
!! ***** | |\  | | | | |  | ||_____|/ __/ | |_| ||_____|| |\  | ___) || |___    ****
!! ***** |_| \_||___||_|  |_|      |_____||____/        |_| \_||____/ |_____|   ****
!! *****                                                                        ****
!! *********************************************************************************
!! Ibrahim Jarrah, 2022
!!
include 'include.f90'
program main
  use general
  use user_parameters
  use parameters
  use mesh
  use geometry
  use curvilinear
  use nim
  use bc
  use postprocess
  use nse
  use restart
  use omp_lib
  implicit none

  timing_start=omp_get_wtime()

  call create_geometry
  call fill_global_arrays
  timing_initialize_1=omp_get_wtime()

  call solve_nse

  timing_all=omp_get_wtime()
  initialization_1_time=timing_initialize_1-timing_start
  initialization_2_time=timing_initialize_2-timing_initialize_1
  integration_time=timing_integration-timing_initialize_2
  simulation_time=timing_all-timing_start
  print '("============================================================================")'
  print '("========================== TIMING SUMMARY ==================================")'
  print '("============================================================================")'
  print '("Preparing mesh time = ",f10.5," seconds.")',initialization_1_time
  print '("Initialization time = ",f10.5," seconds.")',initialization_2_time
  print '("Integration    time = ",f10.5," seconds.")',integration_time
  print '("Simulation     time = ",f10.5," seconds.")',simulation_time
  print '("============================================================================")'
  print '("============================================================================")'


end program main
