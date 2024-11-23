open (1, file='grid.ascii')

!-----
! Nodes
!-----
nindx = 2
nd = 2                                      !  Number of dimensions
write (1, '(i1,1x,i1)') nindx, nd
nindx = 10                                  !  Index, =10 for nodes, =12 for cells, =13 for faces
write (1, '(i2)') nindx
icol = 6                                    !  Number of columns
izoneid = 0                                 !  Number of zone, defines if the face is interior or boundary
NNT_start = 1                               !  Starting node
NNT_end = imax*jmax                         !  Ending node
itype = 0
node_nd = 2
write (1, 100) icol, nindx, izoneid, NNT_start, NNT_end, itype, node_nd
izoneid = 8
itype = 1
write (1, '(i2)') nindx
write (1, 100) icol, nindx, izoneid, NNT_start, NNT_end, itype, node_nd

do j = 1, jmax
do i = 1, imax
write (1, '(2(f23.16,1x))') xgrid(i, j), ygrid(i, j)
end do
end do

!-----
! Cells
!-----
nindx = 12
izoneid = 0
itype = 0
nte_type = 0
NTE_start = 1
NTE_end = (imax - 1)*(jmax - 1)
write (1, '(i2)') nindx
write (1, 100) icol, nindx, izoneid, NTE_start, NTE_end, itype, nte_type
izoneid = 9
itype = 1
nte_type = 3
write (1, '(i2)') nindx
write (1, 100) icol, nindx, izoneid, NTE_start, NTE_end, itype, nte_type

!-----
! Faces
!-----
nindx = 13
izoneid = 0
itype = 0
NBC_type = 0
NFace_type = 0
nte_type = 0

!-- Calculate number of faces
!-- Interior faces
nfaces_y = (imax - 2)*(jmax - 1)  ! y normal faces
nfaces_x = (imax - 1)*(jmax - 2)  ! x normal faces
nfaces_int = nfaces_x + nfaces_y

!-- Boundary faces
nfaces_ilo = (jmax - 1)
nfaces_ihi = (jmax - 1)
nfaces_jlo = (imax - 1)
nfaces_jhi = (imax - 1)                      !  Number of faces on jhi

nfaces_tot = nfaces_int + nfaces_ilo + nfaces_ihi + nfaces_jlo + nfaces_jhi

NFace_start = 1
NFace_end = nfaces_tot
icol = 6                                      !  Number of columns
write (1, '(i2)') nindx
write (1, 100) icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type

!----
! interior faces
!----
izoneid = 10
NBC_type = 2
NFace_type = 2
NFace_start = 1
NFace_end = nfaces_int
icol = 6                                      !  Number of columns
write (1, '(i2)') nindx
write (1, 200) icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type

icol = 4                                      !  Number of columns
do j = 1, jmax - 2
do i = 2, imax - 1
! vertical face
nod1 = (j - 1)*(imax) + i
nod2 = nod1 + imax
nel1 = (j - 1)*(imax - 1) + i - 1
nel2 = nel1 + 1
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
!horizontal face
nod1 = nod2
nod2 = nod1 - 1
nel1 = nel1
nel2 = nel1 + (imax - 1)
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
end do
i = imax
!horizontal face
nod1 = (j)*(imax) + i
nod2 = nod1 - 1
nel1 = (j - 1)*(imax - 1) + i - 1
nel2 = nel1 + (imax - 1)
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
end do
j = jmax - 1
do i = 2, imax - 1
! vertical face
nod1 = (j - 1)*(imax) + i
nod2 = nod1 + imax
nel1 = (j - 1)*(imax - 1) + i - 1
nel2 = nel1 + 1
write (1, '(i1,1x,4(i7,1x))')


!----
!  j=jmax boundary - wall            
!----
izoneid = 4
NBC_type = 3 
NFace_type = 2
NFace_start = Nface_end + 1
NFace_end = Nface_start + nfaces_jhi - 1
icol = 6                                      !  Number of columns
write (1, '(i2)') nindx
write (1, 300) icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type

icol = 4                                      !  Number of columns
do i = 1, imax - 1
nod1 = (jmax - 1)*imax + i + 1
nod2 = nod1 - 1
nel1 = (jmax - 2)*(imax - 1) + i
nel2 = 0
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
end do

!----
!  j=1 boundary - wall
!----

izoneid = 4
NBC_type = 3
NFace_type = 2
NFace_start = Nface_end + 1
NFace_end = Nface_start + nfaces_jlo - 1
icol = 6                                      !  Number of columns
write (1, '(i2)') nindx
write (1, 300) icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type

icol = 4                                      !  Number of columns
do i = 1, imax - 1
nod1 = i
nod2 = i + 1
nel1 = i
nel2 = 0
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
end do

!----
!  i=imax boundary - Outflow
!----
izoneid = 5
NBC_type = 0
NFace_type = 2
NFace_start = Nface_end + 1
NFace_end = NFace_start + nfaces_ihi - 1
icol = 6                                      !  Number of columns
write (1, '(i2)') nindx
write (1, 300) icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type

icol = 4                                      !  Number of columns
do j = 1, jmax - 1
nod1 = (j)*imax
nod2 = (j + 1)*imax
nel1 = j*(imax - 1)
nel2 = 0
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
end do

!----
!  i=1 boundary - Inflow
!----
izoneid = 4
NBC_type = 1
NFace_type = 2
NFace_start = Nface_end + 1
NFace_end = NFace_start + nfaces_ilo - 1
icol = 6                                      !  Number of columns
write (1, '(i2)') nindx
write (1, 300) icol, nindx, izoneid, NFace_start, NFace_end, NBC_type, NFace_type

icol = 4                                      !  Number of columns
do j = 1, jmax - 1
nod1 = (j)*imax + 1
nod2 = (j - 1)*imax + 1
nel1 = (j - 1)*(imax - 1) + 1
nel2 = 0
write (1, '(i1,1x,4(i7,1x))') icol, nod1, nod2, nel1, nel2
end do


100     format(i1, 1x, i2, 1x, i1, 1x, i1, 1x, i8, 1x, i1, 1x, i1)
200     format(i1, 1x, i2, 1x, i2, 1x, i1, 1x, i8, 1x, i1, 1x, i1)
300     format(i1, 1x, i2, 1x, i2, 1x, i8, 1x, i8, 1x, i2, 1x, i1)

