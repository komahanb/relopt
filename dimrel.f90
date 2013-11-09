Module dimrel
      implicit none

      double precision, parameter :: pi=4.0d0*datan(1.0d0)

      integer :: id_proc, num_proc, ierr
      integer,parameter::MAXNMCS=1000000

      integer::casemode
      integer::filenum
      integer::fctindx
      integer::dyncyccnt
      integer::OS,evlfnc,dynamics
      integer::probtype(20) ! 1=fixed 2=variable 
      integer::ndimt 

      double precision fmean,fvar,fstd
      double precision, dimension(20)   :: xavg,xstd,fmeanprime,fvarprime,xavgt,xstdt
      double precision,dimension(20,20)::fmeandbleprime,fvardbleprime
      
      double precision,dimension(20)::dat
    
      logical::mainprog
      integer::OUUflag

    End Module dimrel
