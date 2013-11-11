subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fctindx)
!  use dimpce,only:fctindx
  implicit none

  integer  :: ndvart,fctindx
  double precision :: X(ndvart),fobj,dfdD(ndvart),x3
  double precision ::  rho, L, sigmay, pi, p, E, Fs  
  double precision :: DAT(20)

!  call get_f(ndvart,12,x,fobj)
!  call get_df(ndvart,12,x,dfDD)

  !Problem data and other constants
  dat(1)=10.0 !height ref
  dat(2)=1.0e7 !E
  dat(3)=0.1 !gamma
  dat(4)=50.0*pi/180.0
  dat(5)=30000.0

  ! Max constraint values

  !Tensile
  dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
  dat(7)=10000.0    ! psi tensile_sigma2_max=dat(7)
  dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
  !Comessive
  dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
  dat(10)=10000.0   ! psi comp_sigma2_max=dat(10)
  dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
  !Displacement
  dat(12)=0.005    ! in  max_u_disp=dat(12)
  dat(13)=0.005    ! in  max_v_disp=dat(12)
  dat(14)=1.0      ! Factor of safety
  dat(20)=77      ! filenum for PC

  call threebarf(fctindx,dat,x,fobj)
 ! call threebardf(fctindx,dat,x,dfdD)

  return
end subroutine CalcstuffBFGS
