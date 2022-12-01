module postprocess
use mesh
use general
use geometry
contains
subroutine write_data
  implicit none
  call mesh_export()
end subroutine write_data

subroutine mesh_export()
	include "silo_f9x.inc"
	character(len = 40) :: name,name2
	character(len=6)    :: fmt 
	character(len=6)    :: x1 
	logical             :: exist
	integer err, ierr, dims(2), ndims, NX, NY,i,j,ZX,ZY
	parameter (NX = (NMESH_i+1))
	parameter (NY = NMESH_j+1)
	parameter (ZX = NX-1)
	parameter (ZY = NY-1)

	real x(NX,NY), y(NX,NY)
	real velx_a(2*NX*NY)
	real vely_a(2*NX*NY)
	real pressure_a(2*NX*NY)
	integer e_n
    real dtime


	dims=(/NX,NY/)

	do j=1, NMESH_j
		do i=1,NMESH_i
			x(i,j)=MESHES_ARRAY(i,j)%nodes(1)%x
			y(i,j)=MESHES_ARRAY(i,j)%nodes(1)%y
		enddo
		x(NMESH_i+1,j)=MESHES_ARRAY(NMESH_i,j)%nodes(2)%x
		y(NMESH_i+1,j)=MESHES_ARRAY(NMESH_i,j)%nodes(2)%y
	enddo
	do i=1,NMESH_i
		x(i,NMESH_j+1)=MESHES_ARRAY(i,NMESH_j)%nodes(4)%x
		y(i,NMESH_j+1)=MESHES_ARRAY(i,NMESH_j)%nodes(4)%y
	enddo
	x(NMESH_i+1,NMESH_j+1)=MESHES_ARRAY(NMESH_i,NMESH_j)%nodes(3)%x
	y(NMESH_i+1,NMESH_j+1)=MESHES_ARRAY(NMESH_i,NMESH_j)%nodes(3)%y

	do j=1, NMESH_j
		do i=1,NMESH_i
			e_n=(j-1)*(NMESH_i+1)+i
			velx_a(e_n)=edges_u(meshes_edges(1,i,j))
		enddo
		e_n=(j-1)*(NMESH_i+1)+NMESH_i+1
		velx_a(e_n)=edges_u(meshes_edges(2,NMESH_i,j))
	enddo
	do j=1, NMESH_j
		do i=1,NMESH_i
			e_n=(j-1)*(NMESH_i+1)+i+NX*NY
			velx_a(e_n)=edges_u(meshes_edges(3,i,j))
		enddo
	enddo
	do i=1,NMESH_i
		e_n=(NMESH_j-1)*(NMESH_i+1)+i+NX*NY
		velx_a(e_n)=edges_u(meshes_edges(4,i,NMESH_j))
	enddo

	name='test2.silo'
	fmt = '(I6.6)' ! an integer of width 5 with zeros at the left
	write (x1,fmt) STEP ! converting integer to string using a 'internal file'
	name='plots//step'//trim(x1)
	name2='plots/step'//trim(x1)
	name= trim(name)
  dtime=TI
  if (STEP == 1) dtime=0.0
  print*,"==========================================================================="
  print*,"Printing data at time = ", real(dtime,4), "|     Time step = ", STEP
  print*,"==========================================================================="




	ierr = dbcreate(trim(name), LEN_TRIM(name), DB_CLOBBER, DB_LOCAL, &
				  "Comment about the data", 22, DB_PDB, dbfile)
	if(dbfile.eq.-1) then
	  write (6,*) 'Could not create Silo file!\n'
	endif
    dtime=real(dtime,4)
    err  = dbmkoptlist(2, optlistid)
    err = dbaddiopt(optlistid, DBOPT_CYCLE, STEP)
    if(kind(dtime) == 4) then
       err = dbadddopt(optlistid, DBOPT_TIME, dtime)
    else
       err = dbadddopt(optlistid, DBOPT_DTIME, dtime)
    endif
	ndims = 2
	err = dbputqm (dbfile, "quadmesh", 8, "xc", 2, &
	      "yc", 2, "zc", 2, real(x,4), real(y,4), DB_F77NULL, dims, ndims, &
	      DB_FLOAT, DB_NONCOLLINEAR, optlistid, ierr)

	dims=(/NMESH_i,NMESH_j/)

	ndims = 2
	err = dbputqv1(dbfile, "u", 1, "quadmesh", 8, real(meshes_uxy,4), dims, &
	      ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
	err = dbputqv1(dbfile, "v", 1, "quadmesh", 8, real(meshes_vxy,4), dims, &
	      ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)

	err = dbputqv1(dbfile, "p", 1, "quadmesh", 8, real(meshes_p,4), dims, &
	      ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)


	err = dbputqv1(dbfile, "D", 1, "quadmesh", 8, real(meshes_dila,4), dims, &
	ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
	if(is_analytical) then
		err = dbputqv1(dbfile, "u_T", 3, "quadmesh", 8, real(meshes_uT,4), dims, &
		ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
		err = dbputqv1(dbfile, "v_T", 3, "quadmesh", 8, real(meshes_vT,4), dims, &
		ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
		err = dbputqv1(dbfile, "p_T", 3, "quadmesh", 8, real(meshes_pT,4), dims, &
		ndims, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
		err = dbputdefvars(dbfile, "defvars", 7, 1, &
		"vel_T", 5, DB_VARTYPE_VECTOR, "{u_T,v_T}", 9, optlistid, ierr)
	endif

	err = dbputdefvars(dbfile, "defvars", 7, 1, &
	"vel", 3, DB_VARTYPE_VECTOR, "{u,v}", 5, optlistid, ierr)

	err = dbclose(dbfile)

	inquire(file="plots.visit", exist=exist)
	if (exist) then
	open(12, file="plots.visit", status="old", position="append", action="write")
	else
	open(12, file="plots.visit", status="new", action="write")
	end if
	write(12, '(a)') trim(name2)
	close(12)
end subroutine mesh_export


end module postprocess

