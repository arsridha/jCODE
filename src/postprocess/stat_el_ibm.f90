module stat_el_ibm

  ! External modules
  use stat
  use string
  
  implicit none

  ! Global variables
  integer :: nStatPoints
  real(WP) :: xmin, xmax, dx, ymin, ymax, dy
  real(WP), dimension(:), allocatable :: grid1D
 
end module stat_el_ibm


! =========================!
! Initialize 1D statistics !
! =========================!
subroutine stat_el_ibm_setup

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use parser
  use parallel
  use simulation_flags
  use geometry
  use grid

  implicit none

  ! Local variables
  integer :: i
  
  ! Return if IBM is not used
  if (.not. useIBM .and. .not. useParticles) return

  ! Get stat size
  call parser_read('stat ibm bin size', nStatPoints, globalGridSize(1))

  ! Allocate arrays
  allocate(grid1D(nStatPoints))

  ! Generate the 1D grid
  xmin = minval(coordinates(:,1))
  call parallel_min(xmin)
  xmax = maxval(coordinates(:,1))
  call parallel_max(xmax)
  dx = (xmax - xmin) / real(nStatPoints-1, WP)
  do i = 1, nStatPoints
     grid1D(i) = xmin + real(i-1,WP) * dx
  end do
  
  return
end subroutine stat_el_ibm_setup

! ======================!
! Cleanup 1D statistics !
! ======================!
subroutine stat_el_ibm_cleanup

  ! Internal modules
  use stat_el_ibm
  
  implicit none

  if (allocated(grid1D)) deallocate(grid1D)

  return
end subroutine stat_el_ibm_cleanup

! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_el_ibm_compute

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use ibm
  use particle
  use grid
  use state
  use grid_functions
  use grid_levelset
  use math

  implicit none

  ! Local variables
  integer :: i, j
  real(WP), dimension(nStatPoints, nDimensions) :: favreVelocity
  real(WP), dimension(nStatPoints) :: meanRho, volume, totalVolume

  ! Return if IBM is not used
  if (.not. useIBM .and. .not. useParticles) return      
  
  return
end subroutine stat_el_ibm_compute

! =================== !
! Read ibm statistics !
! =================== !
subroutine stat_el_ibm_read

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use parallel
  use ibm

  implicit none

  ! Nothing to do
  
  return
end subroutine stat_el_ibm_read


! ================== !
! Read 1D statistics !
! ================== !
subroutine stat_el_ibm_write

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use time_info
  use parallel
  use fileio
  use math
  use geometry
  use ibm
  use ibm_ghost_point
  use particle
  use particle_exchange
  use grid
  use state
  use grid_functions
  use grid_levelset

  implicit none

  ! Local variables
  character(len=str_medium) :: filename
  integer  :: i, j, n, iunit, ierror
  real(WP), dimension(nDimensions) :: uDoublePrime
  real(WP), dimension(nStatPoints, nDimensions) :: meanVelocity, favreVelocity,       &
       meanVelocity_ibm, meanVelocity_part      
  real(WP), dimension(nStatPoints) :: meanRho, meanUU, meanVV, meanWW, meanUV, meanUW,&
       meanVW, alpha_ibm, alpha_lpt, volume_ibm, volume_lpt, totalVolume,             &
       mean_granTemp, colfreqPartPart, colfreqPartIBM, colMomPartIBM
  
  ! Return if IBM is not used
  if (.not. useIBM .and. .not. useParticles) return
    
  ! Compute mean stats
  meanRho = 0.0_WP
  meanVelocity = 0.0_WP
  favreVelocity = 0.0_WP
  totalVolume = 0.0_WP
  volume_ibm = 0.0_WP
  volume_lpt = 0.0_WP
  mean_granTemp = 0.0_WP
  colfreqPartPart = 0.0_WP
  colfreqPartIBM = 0.0_WP
  colMomPartIBM = 0.0_WP
  meanVelocity_part = 0.0_WP
  meanVelocity_ibm = 0.0_WP
  
  do i = 1, nGridPoints
    j = 1 + nint((coordinates(i,1) - xmin) / dx)
     ! Sum up grid volumes
     volume_ibm(j) = volume_ibm(j) + gridNorm(i,1)
     volume_lpt(j) = volume_lpt(j) + primitiveGridNorm(i,1) * (1.0_WP - volumeFraction(i,1))
     totalVolume(j) = totalVolume(j) + primitiveGridNorm(i,1)
     ! Mean stats
     meanRho(j) = meanRho(j) + gridNorm(i,1) * conservedVariables(i,1)                     
     meanVelocity(j,1:nDimensions) = meanVelocity(j,1:nDimensions) +                         &
          gridNorm(i,1) * velocity(i,1:nDimensions)
     favreVelocity(j,1:nDimensions) = favreVelocity(j,1:nDimensions) +                       &
          gridNorm(i,1) * conservedVariables(i,2:nDimensions+1)
     meanVelocity_part(j,1:nDimensions) = meanVelocity_part(j,1:nDimensions) +               &
          primitiveGridNorm(i,1) * particleVelocity(i,1:nDimensions) *                       &
          (1.0_WP - volumeFraction(i,1))
     meanVelocity_ibm(j,1:nDimensions) =  meanVelocity_ibm(j,1:nDimensions) +                &
          primitiveGridNorm(i,1) * (1.0_WP - indicatorFunction(i,1)) *                       &
          ibmVelocity(i,1:nDimensions)
     ! Collision frequency and granular temperature
     mean_granTemp(j) = mean_granTemp(j) + granularTemperature(i,1) *                        &
          primitiveGridNorm(i,1) * (1.0_WP - volumeFraction(i,1) )
     colfreqPartPart(j) = colfreqPartPart(j) + collisionFrequency(i,1) *                     &
          primitiveGridNorm(i,1) * (1.0_WP - volumeFraction(i,1) )
     colfreqPartIBM(j) = colfreqPartIBM(j) + collisionFrequencyPartIBM(i,1) *                &
          primitiveGridNorm(i,1) * (1.0_WP - volumeFraction(i,1) )
     colMomPartIBM(j) = colMomPartIBM(j) + collisionMomentum(i,1) *                          &
          primitiveGridNorm(i,1) * (1.0_WP - volumeFraction(i,1) )
  end do
  
  ! Sum them over procs
  call parallel_sum(volume_ibm)
  call parallel_sum(volume_lpt)
  call parallel_sum(totalVolume)
  call parallel_sum(meanRho)    
  call parallel_sum(mean_granTemp)
  call parallel_sum(colfreqPartPart)
  call parallel_sum(colfreqPartIBM)
  call parallel_sum(colMomPartIBM)
  do i = 1, nDimensions
     call parallel_sum(meanVelocity(:,i))
     call parallel_sum(favreVelocity(:,i))
     call parallel_sum(meanVelocity_part(:,i))
     call parallel_sum(meanVelocity_ibm(:,i))     
  end do
  
  ! Normalize
  do i = 1, nStatPoints
     alpha_ibm(i) = 1.0_WP - (volume_ibm(i) / totalVolume(i))            
     alpha_lpt(i) = volume_lpt(i) / totalVolume(i)            
     meanRho(i) = meanRho(i) / volume_ibm(i)
     meanVelocity(i,:) = meanVelocity(i,:) / volume_ibm(i)
     favreVelocity(i,:) = favreVelocity(i,:) / volume_ibm(i) / meanRho(i)
     meanVelocity_ibm(i,:) = meanVelocity_ibm(i,:) / volume_ibm(i) 
     if (volume_lpt(i) .gt. 0.0_WP) then 
          mean_granTemp(i) = mean_granTemp(i)/volume_lpt(i)
          colfreqPartPart(i) = colfreqPartPart(i)/volume_lpt(i) 
          colfreqPartIBM(i) = colfreqPartIBM(i)/volume_lpt(i)
          colMomPartIBM(i) = colMomPartIBM(i)/volume_lpt(i)
          meanVelocity_part(i,:) = meanVelocity_part(i,:) / volume_lpt(i) 
     end if
  end do

  ! Compute variance
  meanUU = 0.0_WP; meanVV = 0.0_WP; meanWW = 0.0_WP
  meanUV = 0.0_WP; meanUW = 0.0_WP; meanVW = 0.0_WP  
   
  ! Compute higher order statistics
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     uDoublePrime = velocity(i,1:nDimensions) - favreVelocity(j,1:nDimensions)
     
     ! Reynolds stresses
     meanUU(j) = meanUU(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1)**2
     if (nDimensions .gt. 1) then
        meanVV(j) = meanVV(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(2)**2
        meanUV(j) = meanUV(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1) *  &
             uDoublePrime(2)
     end if
     if (nDimensions .gt. 2) then
        meanWW(j) = meanWW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(3)**2
        meanUW(j) = meanUW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1) *  &
             uDoublePrime(3)
        meanVW(j) = meanVW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(2) *  &
             uDoublePrime(3)
     end if
  end do
  call parallel_sum(meanUU); meanUU = meanUU / meanRho / volume_ibm
  call parallel_sum(meanVV); meanVV = meanVV / meanRho / volume_ibm
  call parallel_sum(meanWW); meanWW = meanWW / meanRho / volume_ibm
  call parallel_sum(meanUV); meanUV = meanUV / meanRho / volume_ibm
  call parallel_sum(meanUW); meanUW = meanUW / meanRho / volume_ibm
  call parallel_sum(meanVW); meanVW = meanVW / meanRho / volume_ibm
  
  ! Root process writes
  if (iRank .eq. iRoot) then

     ! ------- WRITE THE ASCI FILE ------
     
     filename= "stat/stat-el_ibm"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES21.12)') "t=",time
     write(iunit,'(27a21)') "X", "alphaibm", "alphalpt", "Volume", "<U>", "<V>", "<W>",     &
          "<U'U'>", "<V'V'>", "<W'W'>", "<U'V'>", "<V'W'>", "<U'W'>", "<rho>",              &
          "Utilde", "Vtilde", "Wtilde", "Theta" , "colfreqPartPart", "colfreqPartIBM",      &
          "colMomPartIBM", "<Uibm>","<Vibm>","<Wibm>", "<Ulpt>","Vlpt>","<Wlpt>"
     do i = 1, nStatPoints
        write(iunit,'(10000ES21.12E3)') grid1D(i), alpha_ibm(i), alpha_lpt(i), totalVolume(i), &
             meanVelocity(i,1), meanVelocity(i,2), meanVelocity(i,3), meanUU(i), meanVV(i),  &
             meanWW(i), meanUV(i), meanVW(i), meanUW(i), meanRho(i),                         &
             favreVelocity(i,1), favreVelocity(i,2), favreVelocity(i,3), mean_granTemp(i),   &
             colfreqPartPart(i), colfreqPartIBM(i), colMomPartIBM(i), meanVelocity_ibm(i,1), &
             meanVelocity_ibm(i,2), meanVelocity_ibm(i,3), meanVelocity_part(i,1),           &
             meanVelocity_part(i,2), meanVelocity_part(i,3) 
     end do
     close(iclose(iunit))
     
  end if

  ! Log
  call monitor_log("IBM-LPT STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_el_ibm_write

