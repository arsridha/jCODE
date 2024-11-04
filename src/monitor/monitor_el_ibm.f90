module monitor_el_ibm

   ! External modules
   use monitor

   implicit none

   ! Global variables
   integer :: i1, i2
   real(WP) :: threshold
  
end module monitor_el_ibm
  
  
  ! ========================================= !
  ! Setup the routine to monitor a Shock Tube !
  ! ========================================= !
subroutine monitor_el_ibm_setup

   ! Internal modules
   use monitor_el_ibm

   ! External modules
   use parser
   use parallel
   use simulation_flags
   use solver_options
   use geometry
   use grid_functions
   use particle

   implicit none

   ! Local variables
   integer :: i, j, k
   real(WP) :: x1, x2
   real(WP), dimension(:), allocatable :: x

   !if (.not. useIBM .or. trim(simulationName).ne.'el ibm') return

   ! Read from input
   call parser_read('particle threshold', threshold, 0.1_WP)

   ! Set the monitor names
   call monitor_create('el_ibm', 14)
   call monitor_set_header(1, 'N_ibm', 'r')
   call monitor_set_header(2, 'N_lpt', 'r')
   call monitor_set_header(3, 'xp_ibm_1%', 'r')
   call monitor_set_header(4, 'xp_ibm_99%', 'r')
   call monitor_set_header(5, 'xp_lpt_1%', 'r')
   call monitor_set_header(6, 'xp_lpt_99%', 'r')
   call monitor_set_header(7, 'xp_total_1%', 'r')
   call monitor_set_header(8, 'xp_total_99%', 'r')
   call monitor_set_header(9, '<x>_ibm','r')
   call monitor_set_header(10, '<x>_lpt','r')
   call monitor_set_header(11, 's','r')
   call monitor_set_header(12,'ncol_partpart','r')
   call monitor_set_header(13,'ncol_partibm','r')
   
   return
end subroutine monitor_el_ibm_setup
  
! ================================== !
! Compute Shock Tube statistics      !
! ================================== !
subroutine monitor_el_ibm_timestep
    ! Internal modules
    use monitor_el_ibm
  
    ! External modules
    use parallel
    use solver_options
    use grid_functions
    use state
    use particle
    use ibm
    use particle_solver
  
    implicit none
  
    ! Local variables
    integer, parameter :: nbin = 100
    integer :: i, j, k
    real(WP) :: minValue, minPart, maxPart, buf, xleft_part, xright_part,     &
      pdf_part(nbin), minIBM, maxIBM,xleft_ibm,xright_ibm,pdf_ibm(nbin),        &
      minAll, maxAll, xleft_all, xright_all, pdf_all(nbin), meanX_part, meanX_ibm, seg
  
   !if (.not. useIBM .or. trim(simulationName).ne.'el ibm') return

   pdf_ibm = 0.0_WP
   pdf_part = 0.0_WP   
   pdf_all = 0.0_WP
   meanX_ibm=0.0_WP
   meanX_part=0.0_WP

   xleft_ibm = -huge(1.0_WP)
   xright_ibm = -huge(1.0_WP)
   xleft_part = -huge(1.0_WP)
   xright_part = -huge(1.0_WP)
   xleft_all = -huge(1.0_WP)
   xright_all = -huge(1.0_WP) 
      
   ! Determine minimum and maximum IBM particle x-position 
   minIBM =  huge(1.0_WP)
   maxIBM = -huge(1.0_WP)
   minAll = huge(1.0_WP)
   maxAll = -huge(1.0_WP)
   do i = 1, nObjects
      minIBM = min(minIBM, object(i)%position(1))
      maxIBM = max(maxIBM, object(i)%position(1))
   end do

   if (useParticles) then
      ! Determine minimum and maximum Lagrangian particle x-position
      minPart =  huge(1.0_WP)
      maxPart = -huge(1.0_WP)
      do i = 1, nParticles
         minPart = min(minPart, particles(i)%position(1))
         maxPart = max(maxPart, particles(i)%position(1))
      end do
      call parallel_min(minPart)
      call parallel_max(maxPart)

      ! Determine minimum and maximum position of all particles
      minAll=min(minPart, minIBM)
      maxAll=max(maxPart, maxIBM)
   else
      minAll = minIBM
      maxAll = maxIBM
   end if

   ! Generate PDF of IBM particle position
   do i = 1, nObjects
      j = min(floor((object(i)%position(1)-minIBM) / (maxIBM-minIBM) * nbin) + 1, nbin)
      pdf_ibm(j) = pdf_ibm(j) + 1.0_WP
      k = min(floor((object(i)%position(1)-minAll) / (maxAll-minAll) * nbin) + 1, nbin)
      pdf_all(k) = pdf_all(k)+1.0_WP
      meanX_ibm = meanX_ibm + object(i)%position(1) 
   end do
   meanX_ibm = meanX_ibm / real(nObjects, WP)
   pdf_ibm = pdf_ibm / real(nObjects, WP)

   ! Find position of 1st percentile
   i=1; buf=pdf_ibm(1)
   do while (buf.lt.0.01_WP .and. nObjects.ne.0)
      i=i+1
      buf=buf+pdf_ibm(i)
   end do
   xleft_ibm = real(i,WP) / real(nbin,WP) * (maxIBM - minIBM) + minIBM

   ! Find position of 99th percentile
   i=1; buf=pdf_ibm(1)
   do while (buf.lt.0.99_WP .and. nObjects.ne.0)
      i=i+1
      buf=buf+pdf_ibm(i)
   end do
   xright_ibm = real(i,WP) / real(nbin,WP) * (maxIBM - minIBM) + minIBM 
   
   ! Generate PDF of particle position
   if (useParticles) then
      do i = 1, nParticles
         j = min(floor((particles(i)%position(1)-minPart) / (maxPart-minPart) * nbin) + 1, nbin)
         pdf_part(j) = pdf_part(j) + 1.0_WP
         k = min(floor((particles(i)%position(1)-minAll) / (maxAll-minAll) * nbin) + 1, nbin)
         pdf_all(k) = pdf_all(k) + 1.0_WP
         meanX_part = meanX_part + particles(i)%position(1)
      end do
      call parallel_sum(pdf_part)
      call parallel_sum(pdf_all)
      call parallel_sum(meanX_part)
      meanX_part = meanX_part / real(nParticlesGlobal, WP)
      if (nParticlesGlobal.gt.0) pdf_part = pdf_part / real(nParticlesGlobal, WP)

      ! Find position of 1st percentile
      i=1; buf=pdf_part(1)
      do while (buf.lt.0.01_WP .and. nParticlesGlobal.ne.0)
         i=i+1
         buf=buf+pdf_part(i)
      end do
      xleft_part = real(i,WP) / real(nbin,WP) * (maxPart - minPart) + minPart

      ! Find position of 99th percentile
      i=1; buf=pdf_part(1)
      do while (buf.lt.0.99_WP .and. nParticlesGlobal.ne.0)
         i=i+1
         buf=buf+pdf_part(i)
      end do
      xright_part = real(i,WP) / real(nbin,WP) * (maxPart - minPart) + minPart

      if (nParticlesGlobal.gt.0 .and. nObjects.gt.0) then
         pdf_all = pdf_all / real(nParticlesGlobal + nObjects, WP)
      end if

      ! Find position of 1st percentile
      i=1; buf=pdf_all(1)
      do while (buf.lt.0.01_WP .and. nParticlesGlobal.ne.0 .and. nObjects.ne.0)
         i=i+1
         buf=buf+pdf_all(i)
      end do
      xleft_all = real(i,WP) / real(nbin,WP) * (maxAll - minAll) + minAll

      ! Find position of 99th percentile
      i=1; buf=pdf_all(1)
      do while (buf.lt.0.99_WP .and. nParticlesGlobal.ne.0 .and. nObjects.ne.0)
         i=i+1
         buf=buf+pdf_all(i)
      end do
      xright_all = real(i,WP) / real(nbin,WP) * (maxAll - minAll) + minAll

      seg= meanX_part/meanX_ibm - 1.0_WP
   else
      seg = 0.0_WP
      nParticlesGlobal = 0
      nParticleParticleCollisions = 0
      nParticleIBMCollisions = 0
      meanX_part = 0.0_WP
   end if

   !print *, 'Past particles pdf gen'
   ! Set the shock tube parameters
   call monitor_select('el_ibm')
   call monitor_set_single_value(1, real(nObjects, WP))
   call monitor_set_single_value(2, real(nParticlesGlobal, WP))
   call monitor_set_single_value(3, xleft_ibm)
   call monitor_set_single_value(4, xright_ibm)
   call monitor_set_single_value(5, xleft_part)
   call monitor_set_single_value(6, xright_part)
   call monitor_set_single_value(7, xleft_all)
   call monitor_set_single_value(8, xright_all)
   call monitor_set_single_value(9, meanX_ibm)
   call monitor_set_single_value(10, meanX_part)
   call monitor_set_single_value(11, seg)
   call monitor_set_single_value(12, real(nParticleParticleCollisions,WP))
   call monitor_set_single_value(13, real(nParticleIBMCollisions,WP))
   
   return
  end subroutine monitor_el_ibm_timestep
  
