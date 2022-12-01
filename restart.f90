module restart
use mesh
use general
use user_parameters
use parameters
implicit none
contains

subroutine write_restart
    implicit none
	character (35) :: restart_file
	character (20) :: rest_id

	write (rest_id, *) STEP
	!Create a dynamic file name
	restart_file = "restart/rest.step"//trim(adjustl(rest_id))
	print*,"***************************************************************************"
	print*,"Writing restart file at time = ", TI, "|     Time step = ", STEP
	print*,"***************************************************************************"
	OPEN(111,FILE=restart_file,FORM='UNFORMATTED')
	REWIND(111)
	write(111) TI
	write(111) STEP
	write(111) edges_u
	write(111) edges_v
	write(111) edges_p
	write(111) meshes_uxy
	write(111) meshes_vxy
	write(111) meshes_uxy_old
	write(111) meshes_vxy_old
	write(111) meshes_dila_old
	close(111)
end subroutine write_restart


subroutine read_restart
    implicit none
	character (35) :: restart_file
	character (20) :: rest_id

	write (rest_id, *) read_restart_step
	!Create a dynamic file name
	restart_file = "restart/rest.step"//trim(adjustl(rest_id))
	OPEN(111,FILE=restart_file,FORM='UNFORMATTED')
	REWIND(111)
	read(111) TI
	read(111) STEP
	read(111) edges_u
	read(111) edges_v
	read(111) edges_p
	read(111) meshes_uxy
	read(111) meshes_vxy
	read(111) meshes_uxy_old
	read(111) meshes_vxy_old
	read(111) meshes_dila_old
	close(111)
end subroutine read_restart



end module restart