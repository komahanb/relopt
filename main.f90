program reliability
  use dimrel

  implicit none

  include 'mpif.h'

  !*******************************************
  !          VARIABLE DECLARATIONS           !
  !*******************************************


  integer::i,j,k,ii,jj,kk,iii,jjj,kkk
  integer::seed
  integer*8::NMCS
  integer::ndim
  double precision, allocatable, dimension(:)   :: MNCf
  double precision, allocatable, dimension(:,:) :: MNCx
  double precision :: dinvnorm
  double precision :: pfail
  integer :: is,ie,idec

  integer::beta

  integer:: idx,idxglb,id

  !*******************************************
  !          MAIN PROGRAM                    !
  !*******************************************

  call MPI_START

  ndim=1  
  NMCS=100000000

  beta= 1 !reliability index

  xavg(:)=2.0d0
  xstd(:)=0.2d0
  
  ! Allocate memory for MC samples and corresponding f(x)
  
  allocate(MNCf(NMCS))
  allocate(MNCx(ndim,NMCS))

  ! Master produces NMCS samples using LHS

  if (id_proc.eq.0) then
     call get_seed(seed)
     call latin_random(ndim,NMCS,seed,MNCx)
  end if

  ! Broadcast MNCS to all threads

  call MPI_BCAST(MNCx,NMCS*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! Each thread now normally distributes the LHS sample

  idec = dble(NMCS)/dble(num_proc)
  is   = idec*id_proc + 1
  ie   = idec*(id_proc+1)
  if(id_proc.eq.num_proc-1)ie = NMCS 

  ! Normally distribute with given mean and standard deviation

  do j = is,ie !1, NMCS
     do k=1,ndim 
        MNCx(k,j)=xavg(k)+dinvnorm(MNCx(k,j))*xstd(k)
     end do
  end do

  ! Information Sharing
  do id=0,num_proc-1
     is   = idec*id + 1
     ie   = idec*(id+1)
     if(id.eq.num_proc-1)ie = NMCS
     call MPI_BCAST(MNCx(:,is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
  end do


  idec = dble(NMCS)/dble(num_proc)
  is   = idec*id_proc + 1
  ie   = idec*(id_proc+1)
  if(id_proc.eq.num_proc-1)ie = NMCS 

  if (id_proc.eq.0)print*,"NMCS = ",NMCS
  if (id_proc.eq.0)print*,'====================================================================='
  if (id_proc.eq.0) write(*,'(4a)') '   Rel. Index','  Violated   ','   P_fail  ','                P_k=1.0-P_fail'
if (id_proc.eq.0)print*,'======================================================================'
  do beta=0,6 

     idx=0

     do ii=is,ie

        if ((MNCx(1,ii)).gt. (xavg(1)+dble(beta)*xstd(1)) ) idx=idx+1
        ! we are having a linear limit state function, that says failure is when the samples fall below beta*sigma 's

     end do

     call MPI_ALLREDUCE(idx,idxglb,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

     pfail=dble(idxglb)/dble(nmcs)

     if (id_proc.eq.0) print*,beta,idxglb,pfail,1.0-pfail

  end do
if (id_proc.eq.0)print*,'======================================================================'

  if (id_proc.eq.0) then

     if (NMCS.lt.100000) then

        open(120,file='MCsamp.dat',form='formatted',status='unknown')
        !     write(120,'(2i8)') NMCS,ndim
        !     write(120,*) (xavg(i),i=1,ndim)
        !     write(120,*) (xstd(i),i=1,ndim)
        do j=1,NMCS
           write(120,*) (MNCx(k,j),k=1,ndim)
        end do
        close(120)
     end if

  end if


  call MPI_Barrier(MPI_COMM_WORLD,ierr) ! All wait until master gets here
  !*******************************************
  !          END OF MAIN PROGRAM             !
  !*******************************************

!!$  if (id_proc.eq.0) then
!!$     write(filenum,*)
!!$     write(filenum,*)'>> Program terminated successfully'
!!$     write(filenum,*) 
!!$  end if

  deallocate(MNCf,MNCx) 
  call stop_all

end program reliability

!+++++++++++++++++++++++++++++++++

!!$
!!$
!!$subroutine get_f()
!!$  implicit none
!!$  
!!$  
!!$  
!!$  return
!!$end subroutine get_f

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine CDF(xin,xc,st,vout)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: vout
  double precision :: vtmp
  !       vout = 0.5d0 * (1.d0 + erf( (xin-xc)/(st*dsqrt(2.d0)) ))
  call ERF_MINE1( (xin-xc)/(st*dsqrt(2.d0)), vtmp )
  vout = 0.5d0 * (1.d0 + vtmp)
end subroutine CDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine DCDF(xin,xc,st,dvout)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: dvout
  double precision :: dvtmp
  call DERF_MINE( (xin-xc)/(st*dsqrt(2.d0)), dvtmp )
  dvout = 0.5d0*dvtmp/(st*dsqrt(2.d0))
end subroutine DCDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine CDF_Numerical(xin,xc,st,cdf)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: cdf
  double precision :: vtmp
  integer :: i,num
  double precision :: xs,xe,dx,x1,x2,pdf1,pdf2
  if(xin.lt.xc)then
     cdf = 0.d0
     xs  = xin -2.d0
     xe  = xin
  else if(xin.ge.xc)then
     cdf = 0.5d0
     xs  = xc
     xe  = xin
  end if
  num = 1001
  dx  = (xe-xs)/dble(num-1)
  do i=1,num-1
     x1 = xs + dble(i-1)*dx
     x2 = xs + dble(i  )*dx
     call normal_dist(x1,xc,st,pdf1)
     call normal_dist(x2,xc,st,pdf2)
     cdf = cdf + (pdf1+pdf2)*dx*0.5d0
  end do
  return
end subroutine CDF_Numerical
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine normal_dist(xin,xc,st,y)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: y
  double precision :: pi

  pi = 4.d0*datan(1.d0)
  y = exp(-1.d0*(xin-xc)**2/2.d0/st/st)/dsqrt(2.d0*pi*st**2)
  return
end subroutine normal_dist

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ERF_MINE1(xin,yout)
  implicit none
  double precision, intent(in)  :: xin
  double precision, intent(out) :: yout
  integer :: i,k,n
  double precision :: vsum,kai
  ! n is the order of Taylor
  ! Maybe accurate within the range of [-4:4] with n=100
  n = 100
  vsum = 0.d0
  do i=0,n
     kai = 1.d0
     do k=1,i
        kai = kai * (-1.d0) * xin**2 / dble(k)
     end do
     vsum = vsum + kai*xin/dble(2*i+1)
  end do
  yout = vsum*2.d0/(dsqrt(3.141592653589793238d0))

  if(yout.gt.1.d0)write(*,'(a,2e15.5)')'*ERF>1 ',xin,yout-1.d0
end subroutine ERF_MINE1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine DERF_MINE(xin,dyout)
  implicit none
  double precision, intent(in)  :: xin
  double precision, intent(out) :: dyout
  double precision :: vsum

  vsum  = exp(-1.d0*xin**2)
  dyout = vsum*2.d0/(dsqrt(3.141592653589793238d0))
  
end subroutine DERF_MINE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine find_st(pst,dx,st)
  implicit none
  ! find st which ensures the probability within [-dx:dx] is pst
  double precision, intent(in)  :: pst,dx
  double precision, intent(out) :: st
  integer :: iter
  double precision :: pout,s,ds,vv,vp,dv,sini

  if(pst.le.0.d0.or.pst.ge.1.d0)stop'pst in find_st'
  pout = (1.d0-pst)/2.d0
  sini = 1.d0
190 continue
  s    = sini
  ds   = sini*1.d-3
  iter = 0
200 continue
  iter = iter + 1
  call CDF(-1.d0*dx,0.d0,s,   vv)
  if(dabs(vv-pout).le.1.d-10)go to 210
  call CDF(-1.d0*dx,0.d0,s+ds,vp)
  dv = (vp-vv)/ds
  !       write(filenum,'(5e15.5)')s,vv,pout,vv-pout,dv
  if(dv.eq.0.d0)stop'dv = 0.d0 in find_st'
  if(iter.ge.100)stop'iter>100 in find_st'
  s = s - (vv-pout)/dv
  if(s.le.0.d0)then
     sini = sini * 0.1d0
     go to 190
  end if
  go to 200
210 continue
  !       write(filenum,'(4e15.5)')s,vv,pout,vv-pout
  st = s

end subroutine find_st
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine find_x(mode,ran,xc,st,xout)
  implicit none
  ! find xout which ensures the CDF at xout is ran
  ! mode=0 : analytical CDF (fast but less robust?)
  ! mode=1 : numerical  CDF (time comsuming, but robust)
  integer, intent(in) :: mode
  double precision, intent(in)  :: ran,xc,st
  double precision, intent(out) :: xout
  integer :: iter 
  double precision :: x,vv,dv,vp,dx

  x    = xc
  dx   = 1.d-4
  iter = 0
200 continue
  iter = iter + 1
  if(mode.eq.0)then
     call CDF(x,xc,st,vv)
  else
     call CDF_Numerical(x,xc,st,vv)
  end if
  if(dabs(vv-ran).le.1.d-10)go to 210
  if(mode.eq.0)then
     call CDF(x+dx,xc,st,vp)
  else
     call CDF_Numerical(x+dx,xc,st,vp)
  end if
  dv = (vp-vv)/dx
  !       write(filenum,'(4e15.5)')x,vv,vv-ran,dv
  if(dv.eq.0.d0)stop'dv=0 in find_x'
  if(iter.ge.100)stop'iter>100 in find_x'
  x = x - (vv-ran)/dv
  go to 200
210 continue
  !       write(filenum,'(3e15.5)')x,vv,vv-ran
  xout = x

end subroutine find_x
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function DINVNORM(p)

  implicit none

  real*8 :: dinvnorm,p,p_low,p_high
  real*8 :: a1,a2,a3,a4,a5,a6
  real*8 :: b1,b2,b3,b4,b5
  real*8 :: c1,c2,c3,c4,c5,c6
  real*8 :: d1,d2,d3,d4
  real*8 :: z,q,r
  a1=-39.6968302866538
  a2=220.946098424521
  a3=-275.928510446969
  a4=138.357751867269
  a5=-30.6647980661472
  a6=2.50662827745924
  b1=-54.4760987982241
  b2=161.585836858041
  b3=-155.698979859887
  b4=66.8013118877197
  b5=-13.2806815528857
  c1=-0.00778489400243029
  c2=-0.322396458041136
  c3=-2.40075827716184
  c4=-2.54973253934373
  c5=4.37466414146497
  c6=2.93816398269878
  d1=0.00778469570904146
  d2=0.32246712907004
  d3=2.445134137143
  d4=3.75440866190742
  p_low=0.02425
  p_high=1-p_low
  if(p.lt.p_low) goto 201
  if(p.ge.p_low) goto 301
201 q=dsqrt(-2*dlog(p))
  z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
  goto 204
301 if((p.ge.p_low).and.(p.le.p_high)) goto 202
  if(p.gt.p_high) goto 302
202 q=p-0.5
  r=q*q
  z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
  goto 204
302 if((p.gt.p_high).and.(p.lt.1)) goto 203
203 q=dsqrt(-2*dlog(1-p))
  z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204 dinvnorm=z

  return

end function dinvnorm

