! =======================================================================
      program Heisenberg_3D
      implicit none
      
      
      
!-----------MPI-------------
      include 'mpif.h'
      integer myid, ierror, num_of_procs, tag, source 
      integer destination, count
      integer status(MPI_STATUS_SIZE)
      double precision Energy_current
      integer master

! ------------------------------
      integer N_Spins, N_Couples
      
      double precision, allocatable :: S_Value(:), Exchange(:), 
     &                                 Conf(:, :)
      integer, allocatable :: Couple(:, :)
! ------------------------------
      integer N_Frozen_Spins
      integer, allocatable ::          Frozen_Spin(:)
      double precision, allocatable :: Frozen_Spin_Conf(:, :)
      
      integer, allocatable :: Spin_Type(:)
     
      integer N_Fitting_Spins
      integer, allocatable :: Fitting_Spin(:),
     &                        Fitting_Spin_NumByID(:)     
      integer Fit
! ------------------------------
      integer Hmag_Dir
      double precision Hmag
! ------------------------------
      integer N_Steps
! ------------------------------
      double precision Pi
      
      integer N_Adjusting_Parameters
      
      double precision Energy
      double precision Energy_Coefficient
      integer N_States
      double precision, allocatable :: Spectrum(:)
! ------------------------------
      integer Lessened_output
      integer Write_State_or_not
! ------------------------------
      character*1 Nsp1
      character*2 Nsp2
      character*3 Nsp3
      character*4 Nsp4
      character*5 Nsp5
! ------------------------------
      integer debug
      integer i
      integer num
      
      integer step
! ++++++++++++++++++++++++++++++
      
! --- Welcoming message ---
      write(*,*)
      write(*,*) '  Hello! It is Heisenberg 3D written by Kashin I. V.'
      write(*,*) '    Current version is 19 apr 2018 - MaxFn 10^6  '
      write(*,*)
! -------------------------

! --- Lets prepare the Output folder ---
      call system('mkdir Output')
      call system('rm ./Output/*')
      write(*,*)
! --------------------------------------
      
! ==== Lets read the parameters of the 3D Heisenberg model ====
! --- The number of spins and its absolute values ---
      open(unit = 1, file = './Input/Basic', status = 'old')
      read(1, *) ! The number of spins
      read(1, *) N_Spins
      read(1, *) ! Magnetic field
      read(1, *) Hmag_Dir, Hmag
      read(1, *) ! Energy Coefficient, turning H to [meV]
      read(1, *) Energy_Coefficient
      read(1, *) ! Number of steps
      read(1, *) N_Steps
      read(1, *) ! Lessened output mode
      read(1, *) Lessened_output
      read(1, *) ! Debug mode
      read(1, *) debug
      
! --- Allocate the Spin Value array ----
      allocate(S_Value(N_Spins))
! --------------------------------------

      read(1, *) ! Absolute values
      do i = 1, N_Spins
       read(1, *) num, S_Value(i)
      enddo !i
      close(1)
      
      write(*, '("   Number of spins: ", i4)') N_Spins
      write(*,*)
      if (Hmag_Dir.eq.1) then
       write(*, '("   Magnetic field: ", f6.3, " eV, X-directed")')
     &                                   Hmag
      endif !Magnetic Field is directed along X axis
      if (Hmag_Dir.eq.2) then
       write(*, '("   Magnetic field: ", f6.3, " eV, Y-directed")')
     &                                   Hmag
      endif !Magnetic Field is directed along Y axis
      if (Hmag_Dir.eq.3) then
       write(*, '("   Magnetic field: ", f6.3, " eV, Z-directed")')
     &                                   Hmag
      endif !Magnetic Field is directed along Z axis
      write(*,*)
      write(*, '("   Energy coefficient: ", f10.3)') 
     &               Energy_Coefficient
      write(*,*)
      write(*, '("   Number of steps: ", i8)') N_Steps
      write(*,*)
      if (Lessened_output.eq.0) then 
       write(*,*) '  Full output stream   '
      else
       write(*,*) '  Lessened output stream   '
      endif !debug mode
      write(*,*)
      if (debug.eq.0) then 
       write(*,*) '  Debug mode disabled   '
      else
       write(*,*) '  Debug mode enabled   '
      endif !debug mode
      write(*,*)
      write(*,*) '  Absolute values: '
      do i = 1, N_Spins
       write(*, '(i6, f8.3)') i, S_Value(i)
      enddo !i

      write(*,*)
! -----------------------------------------------

! --- Heisenberg exchange interactions -----------------------
      open(unit = 2, file = './Input/Exchange', status = 'old')
      read(2, *) ! The number of couples
      read(2, *) N_Couples
      
! --- Allocate the Exchange couplings arrays ----
      allocate(Exchange(N_Couples))
      allocate(Couple(N_Couples, 2))
! -----------------------------------------------
      
      read(2, *)
      do i = 1, N_Couples
       read(2, *) num, Couple(i, 1), Couple(i, 2), Exchange(i)
      enddo !i
      close(2)
      
      write(*,*) '  Exchange interactions (eV):   '
      do i = 1, N_Couples
       write(*, '(2i6, i4, f14.6)') 
     &       i, Couple(i, 1), Couple(i, 2), Exchange(i)
      enddo !i
      write(*,*)
! ------------------------------------------------------------

! --- Frozen Spins configuration -----------------------------
      open(unit = 3, file = './Input/Frozen', status = 'old')
      read(3, *) ! The number of Frozen Spins
      read(3, *) N_Frozen_Spins
      
      if (N_Frozen_Spins.gt.0) then 
! --- Allocate the Frozen Spins arrays ----
       allocate(Frozen_Spin(N_Frozen_Spins))
       allocate(Frozen_Spin_Conf(N_Frozen_Spins, 2))
! -----------------------------------------
      
       read(3, *) ! No., ID, Configuration (Long, Trans)
       do i = 1, N_Frozen_Spins
        read(3, *) num, Frozen_Spin(i),
     &                  Frozen_Spin_Conf(i, 1),
     &                  Frozen_Spin_Conf(i, 2)
       enddo !i
       close(3)
      
       write(*,*) '  Frozen Spins configuration
     & (No., ID, Long, Trans): '
       do i = 1, N_Frozen_Spins
        write(*, '(2i6, 2f14.6)') i, Frozen_Spin(i),
     &                               Frozen_Spin_Conf(i, 1),
     &                               Frozen_Spin_Conf(i, 2)
       enddo !i
       write(*,*)
      else
       close(3)
       write(*,*) '  No Frozen Spins  '
       write(*,*)
      endif !N_Frozen_Spins>0
! ------------------------------------------------------------
! =============================================================

! --- Prepare the Spin_Type array ---------
      allocate(Spin_Type(N_Spins))
      Spin_Type = 1
      do i = 1, N_Frozen_Spins
       Spin_Type(Frozen_Spin(i)) = 0
      enddo !i
! -----------------------------------------

! --- Prepare the Fitting Spins arrays ----
      N_Fitting_Spins = N_Spins - N_Frozen_Spins
      allocate(Fitting_Spin(N_Fitting_Spins))
      allocate(Fitting_Spin_NumByID(N_Spins))
      Fitting_Spin_NumByID = 0
      
      Fit = 0
      do i = 1, N_Spins
       if (Spin_Type(i).eq.1) then
        Fit = Fit + 1
        Fitting_Spin(Fit) = i
        Fitting_Spin_NumByID(i) = Fit
       endif !If Fitting
      enddo !i
      if (Fit.ne.N_Fitting_Spins) stop 'Fit Initiation Error'
      
      if (N_Frozen_Spins.gt.0 .or. debug.gt.0) then 
       write(*,*) '  Fitting Spins (No., ID): '
       do i = 1, N_Fitting_Spins
        write(*, '(2i6)') i, Fitting_Spin(i)
       enddo !i
       write(*,*)
      endif !Only if we have frozen spins (or debug mode)
! -----------------------------------------


! --- Allocate the Configuration array ----
      allocate(Conf(N_Spins, 2))
! -----------------------------------------

! --- Allocate the Spectrum array ----
      allocate(Spectrum(N_Steps))
      Spectrum = 0.d0
      N_States = 0
! ------------------------------------

! --- Create output files ----
      open( unit = 4,
     &      file = './Output/Longitudinal',
     &      status = 'replace' )
     
      open( unit = 5,
     &      file = './Output/Transverse',
     &      status = 'replace' )
     
      open( unit = 8,
     &      file = './Output/States_Guide',
     &      status = 'replace' )
     
      if (debug.gt.0) then
       open( unit = 101,
     &       file = './Output/Initial_Configurations',
     &       status = 'replace' )
     
       write(101, *) '-----------------------------------------'
       write(101, *) '     No. | Longitudinal | Transverse  '
       write(101, *) '-----------------------------------------'
      endif !debug
! -----------------------------

! --- Interpret N_Spins as a character ---
      if (N_Spins.lt.10) then
       write(Nsp1, '(i1)') N_Spins
       if (debug.gt.0) 
     &     write(*,*) '  Output format: (' // Nsp1 // 'f10.1)'
       if (debug.gt.0) write(*,*)
      endif !Np1
      
      if (N_Spins.ge.10 .and. N_Spins.lt.100) then
       write(Nsp2, '(i2)') N_Spins
       if (debug.gt.0) 
     &     write(*,*) '  Output format: (' // Nsp2 // 'f10.1)'
       if (debug.gt.0) write(*,*)
      endif !Np2
      
      if (N_Spins.ge.100 .and. N_Spins.lt.1000) then
       write(Nsp3, '(i3)') N_Spins
       if (debug.gt.0) 
     &     write(*,*) '  Output format: (' // Nsp3 // 'f10.1)'
       if (debug.gt.0) write(*,*)
      endif !Np3
      
      if (N_Spins.ge.1000 .and. N_Spins.lt.10000) then
       write(Nsp4, '(i4)') N_Spins
       if (debug.gt.0) 
     &     write(*,*) '  Output format: (' // Nsp4 // 'f10.1)'
       if (debug.gt.0) write(*,*)
      endif !Np4
      
      if (N_Spins.ge.10000 .and. N_Spins.lt.100000) then
       write(Nsp5, '(i5)') N_Spins
       if (debug.gt.0) 
     &     write(*,*) '  Output format: (' // Nsp5 // 'f10.1)'
       if (debug.gt.0) write(*,*)
      endif !Np4
      
      if (N_Spins.gt.100000) stop 'Too many spins'
! -----------------------------------------

! --- Set Pi ---
      Pi = dacos(-1.d0)
      if (debug.gt.0) then 
       write(*, '("   Pi = ", f8.6)') Pi
       write(*,*)
      endif !debug
! --------------

! --- Convert Frozen Spins Configuration into Radians ----
       Frozen_Spin_Conf = Frozen_Spin_Conf * Pi / 180.d0
! --------------------------------------------------------

! --- Plant a random seed ---
      call plant_random_seed
! ---------------------------

!-----------MPI-------------
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_of_procs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
!-----------MPI-------------

! -------- Main Step cycle ------------------
      do step = 1, N_Steps

! --- Set the initial chaotic configuration ------
       call random_number(Conf)
! ++++ Longitudinal ++++
       Conf(:, 1) = dacos((Conf(:, 1) - 0.5d0) * 2.d0)
! ++++ Transverse ++++
       Conf(:, 2) = (Conf(:, 2) - 0.5d0) * 2.d0 * Pi
! ------------------------------------------------

! --- Set the Frozen Spins ---------------
       do i = 1, N_Frozen_Spins
        Conf(Frozen_Spin(i), 1) = Frozen_Spin_Conf(i, 1)
        Conf(Frozen_Spin(i), 2) = Frozen_Spin_Conf(i, 2)
       enddo !i
! ----------------------------------------

! --- Write the initial configuration into file if requested ---
      if (debug.gt.0) then
       write(101, '("  Configuration No. ", i4)') step
       write(101, *) '=========================='
       do i = 1, N_Spins
        write(101, '(i8, f12.1, f14.1)') 
     &        i, Conf(i, 1) * 180.d0 / Pi,
     &           Conf(i, 2) * 180.d0 / Pi
       enddo !i
       write(101, *) '-----------------------------------------'
      endif !debug
! --------------------------------------------------------------

! --- Looking for the Energy Minimum ---------
       N_Adjusting_Parameters = 2 * N_Fitting_Spins
       call Search_for_the_Energy_minimum
     &    ( N_Adjusting_Parameters,
     &      N_Spins, N_Couples,
     &      Spin_Type,
     &      N_Fitting_Spins,
     &      Fitting_Spin, Fitting_Spin_NumByID,
     &      S_Value, Exchange, Couple, Conf,
     &      Hmag_Dir, Hmag, Energy )
! --------------------------------------------

       



! --- Convert the Result into Degrees ---
       Conf = Conf * 180.d0 / Pi
! ---------------------------------------

! --- Check the results for period ------
       do i = 1, N_Spins
!  +++++ Longitudinal +++++
        do while (Conf(i, 1).lt.-180.d0)
         Conf(i, 1) = Conf(i, 1) + 360.d0
        enddo !if <0
        
        do while (Conf(i, 1).gt.180.d0)
         Conf(i, 1) = Conf(i, 1) - 360.d0
        enddo !if >180
        
!  +++++ Transverse +++++
        do while (Conf(i, 2).lt.-180.d0)
         Conf(i, 2) = Conf(i, 2) + 360.d0
        enddo !if <-180
        
        do while (Conf(i, 2).gt.180.d0)
         Conf(i, 2) = Conf(i, 2) - 360.d0
        enddo !if >180
       enddo !i
! ---------------------------------------

! --- Transform negative Longitudinal angle ---
       do i = 1, N_Spins
        if (Conf(i, 1).lt.0.d0) then
         Conf(i, 1) = -Conf(i, 1)
         
         if (Conf(i, 2).gt.0.d0) then
          Conf(i, 2) = Conf(i, 2) - 180.d0
         else
          Conf(i, 2) = Conf(i, 2) + 180.d0
         endif !Transverse correction
        endif !Client
       enddo !i
! ---------------------------------------------

!-----------MPI-------------

      integer current_id
      if (myid.ne.master) then
         call MPI_SEND(Energy, count, MPI_DOUBLE_PRECISION, master, tag,
     & MPI_COMM_WORLD, ierror)
         call MPI_RECV(New_or_not, count, MPI_INTEGER, master, tag,
     & MPI_COMM_WORLD, status, ierror)
      endif
      
      if (myid.eq.master) then
        do current_id = num_of_procs-1, 0, -1
            if (current_id.gt.master) then
                call MPI_RECV(Energy_current, count, 
     & MPI_DOUBLE_PRECISION, master, tag, MPI_COMM_WORLD, 
     & status, ierror)
           else
            Energy_current = Energy
           endif
           
! --- !!!!!!!!!!!!!!!!! Lets find out novelty of the Energy value ---
           New_or_not = 1
           do st = 1, N_States
            if (dabs(Energy - Spectrum(st)).lt.1.d-5) then
             New_or_not = 0
             exit
            endif !Not new
           enddo !st
! -------------------------------------------------  
! --- Perform the gained state analysis --------         
           if (New_or_not.eq.1) then call State_Analysis
     &    ( Energy, N_Spins, Conf, S_Value, Hmag_Dir,
     &      Lessened_output, Write_State_or_not,
     &      N_States, N_Steps, Spectrum )
! ----------------------------------------------
              if (current_id.gt.master) then 
              call MPI_SEND(Energy, count, 
          & MPI_DOUBLE_PRECISION, master, tag,
          & MPI_COMM_WORLD, ierror)
             enddo !current_id
            endif !myid
    
           

!-----------STILL MPI-------------



! --- Perform the gained state analysis --------
!       call State_Analysis
!     &    ( Energy, N_Spins, Conf, S_Value, Hmag_Dir,
!     &      Lessened_output, Write_State_or_not,
!     &      N_States, N_Steps, Spectrum )
! ----------------------------------------------
       if (New_or_not.eq.1) then
! --- Write the output line ---------------
       if (Write_State_or_not.eq.1) then
       
        if (N_Spins.lt.10) then
         write(4, '(i6, f18.6, "    ", ' // Nsp1 // 'f10.1)') 
     &         step, Energy, 
     &        (Conf(i, 1), i = 1, N_Spins)
     
         write(5, '(i6, f18.6, "    ", ' // Nsp1 // 'f10.1)')  
     &         step, Energy, 
     &        (Conf(i, 2), i = 1, N_Spins)
        endif !Nsp1
       
        if (N_Spins.ge.10 .and. N_Spins.lt.100) then
         write(4, '(i6, f18.6, "    ", ' // Nsp2 // 'f10.1)') 
     &         step, Energy, 
     &        (Conf(i, 1), i = 1, N_Spins)
     
         write(5, '(i6, f18.6, "    ", ' // Nsp2 // 'f10.1)')  
     &         step, Energy, 
     &        (Conf(i, 2), i = 1, N_Spins)
        endif !Nsp2
       
        if (N_Spins.ge.100 .and. N_Spins.lt.1000) then
         write(4, '(i6, f18.6, "    ", ' // Nsp3 // 'f10.1)') 
     &         step, Energy, 
     &        (Conf(i, 1), i = 1, N_Spins)
     
         write(5, '(i6, f18.6, "    ", ' // Nsp3 // 'f10.1)')  
     &         step, Energy, 
     &        (Conf(i, 2), i = 1, N_Spins)
        endif !Nsp3
       
        if (N_Spins.ge.1000 .and. N_Spins.lt.10000) then
         write(4, '(i6, f18.6, "    ", ' // Nsp4 // 'f10.1)') 
     &         step, Energy, 
     &        (Conf(i, 1), i = 1, N_Spins)
     
         write(5, '(i6, f18.6, "    ", ' // Nsp4 // 'f10.1)')  
     &         step, Energy, 
     &        (Conf(i, 2), i = 1, N_Spins)
        endif !Nsp4
       
        if (N_Spins.ge.10000 .and. N_Spins.lt.100000) then
         write(4, '(i6, f18.6, "    ", ' // Nsp5 // 'f10.1)') 
     &         step, Energy, 
     &        (Conf(i, 1), i = 1, N_Spins)
     
         write(5, '(i6, f18.6, "    ", ' // Nsp5 // 'f10.1)')  
     &         step, Energy, 
     &        (Conf(i, 2), i = 1, N_Spins)
        endif !Nsp5
        
       endif !If it is Lessened mode - write only the lowest states
! -----------------------------------------
      enddo !step
! ----------------------------------------------------

! --- Write out the Spectrum into file ------
      open( unit = 9,
     &      file = './Output/Spectrum',
     &      status = 'replace' )
     
      write(9, *) '  State (rel. units) | State(meV) |
     & Excitation(meV) | Excitation(K)'
      write(9, *)
      do i = 1, N_States
       write(9, '(3f20.6, f25.3)') 
     &    Spectrum(i),
     &    Spectrum(i) * Energy_Coefficient,
     &  ( Spectrum(i) - Spectrum(1) ) * Energy_Coefficient,
     &  ( Spectrum(i) - Spectrum(1) ) * Energy_Coefficient *
     &                                  11.604522167676d0
      enddo !i
      
      close(9)
! ---------------------------------

! --- Close the output files ---
      close(4)
      close(5)
      close(8)
      if (debug.gt.0) close(101)
! ------------------------------
! MPI
       endif
       
! --------------ENDMPI-------

! --- Deallocation block ---
      deallocate(S_Value)
      deallocate(Exchange)
      deallocate(Couple)
      if (N_Frozen_Spins.gt.0) then
       deallocate(Frozen_Spin)
       deallocate(Frozen_Spin_Conf)
      endif !Frozen
      deallocate(Spin_Type)
      deallocate(Fitting_Spin)
      deallocate(Fitting_Spin_NumByID)
      deallocate(Conf)
      deallocate(Spectrum)
! --------------------------

! --- Final message --------
      write(*,*) '  All done. Have a nice day!'
      write(*,*)
! --------------------------
      end
! =======================================================================
      subroutine State_Analysis
     &         ( Energy, N_Spins, Conf, S_Value, Hmag_Dir,
     &           Lessened_output, Write_State_or_not,
     &           N_States, N_Steps, Spectrum )
      implicit none
! ------------------------------
      double precision Energy
      integer N_Spins
      double precision Conf(N_Spins, 2)
      double precision S_Value(N_Spins)
      integer Hmag_Dir
      
      integer Lessened_output
      integer Write_State_or_not
      
      integer N_States
      integer N_Steps
      double precision Spectrum(N_Steps)
! ------------------------------
      double precision Mag_Moment
      double precision Total_Mag_Moment
! ------------------------------
      integer New_or_not
      integer Place
      
      integer st
      
      character*1 Nst1
      character*2 Nst2
      character*3 Nst3
      character*4 Nst4
      character*5 Nst5
      integer i
      
      double precision Rad
! ------------------------------

! --- Degrees into Radians ---
      Rad = dacos(-1.d0) / 180.d0
! ----------------------------

! --- Depending on the mode, we write the state or not ---
      if (Lessened_output.eq.0) then
       Write_State_or_not = 1
      else
       Write_State_or_not = 0
      endif !Full output - write, Lessened - not to (at this step)
! --------------------------------------------------------

! --- Lets find out novelty of the Energy value ---
!      New_or_not = 1
!      do st = 1, N_States
!       if (dabs(Energy - Spectrum(st)).lt.1.d-5) then
!        New_or_not = 0
!        exit
!       endif !Not new
!      enddo !st
! -------------------------------------------------

! --- Treat the new one ----------------------
      if (New_or_not.eq.1) then
       do Place = 1, N_States
        if (Energy.lt.Spectrum(Place)) exit
       enddo !Place
       
       if (Place.lt.N_States + 1) then
        do st = N_States, Place, -1
         Spectrum(st + 1) = Spectrum(st)
        enddo !st
       endif !Shift all higher states
       
       Spectrum(Place) = Energy !Establish the new state
       N_States = N_States + 1
       
       if (Place.eq.1) Write_State_or_not = 1
       
! +++++ Write the configuration into files +++++++
       if (Lessened_output.eq.0 .or. Place.eq.1) then
       
        if (N_States.lt.10) then
         write(Nst1, '(i1)') N_States
         open( unit = 6,
     &         file = './Output/State_' // Nst1,
     &         status = 'replace' )
         open( unit = 7,
     &         file = './Output/State_' // Nst1 // '_Plot',
     &         status = 'replace' )
        endif ! <10
       
        if (N_States.ge.10 .and. N_States.lt.100) then
         write(Nst2, '(i2)') N_States
         open( unit = 6,
     &         file = './Output/State_' // Nst2,
     &         status = 'replace' )
         open( unit = 7,
     &         file = './Output/State_' // Nst2 // '_Plot',
     &         status = 'replace' )
        endif ! >10 and <100
       
        if (N_States.ge.100 .and. N_States.lt.1000) then
         write(Nst3, '(i3)') N_States
         open( unit = 6,
     &         file = './Output/State_' // Nst3,
     &         status = 'replace' )
         open( unit = 7,
     &         file = './Output/State_' // Nst3 // '_Plot',
     &         status = 'replace' )
        endif ! >100 and <1000
       
        if (N_States.ge.1000 .and. N_States.lt.10000) then
         write(Nst4, '(i4)') N_States
         open( unit = 6,
     &         file = './Output/State_' // Nst4,
     &         status = 'replace' )
         open( unit = 7,
     &         file = './Output/State_' // Nst4 // '_Plot',
     &         status = 'replace' )
        endif ! >1000 and <10000
       
        if (N_States.ge.10000 .and. N_States.lt.100000) then
         write(Nst5, '(i5)') N_States
         open( unit = 6,
     &         file = './Output/State_' // Nst5,
     &         status = 'replace' )
         open( unit = 7,
     &         file = './Output/State_' // Nst5 // '_Plot',
     &         status = 'replace' )
        endif ! >10000 and <100000
       
        if (N_States.ge.100000) stop 'Too many states'
       
        write(6, *) '-----------------------------
     &------------------------------'
        write(6, '("   Energy = ", f30.6)') Energy
        write(6, *) '-----------------------------
     &------------------------------'
        write(6, *) '     No. | Longitudinal | Transverse |
     & Magnetic Moment   '
        write(6, *)
       
        Total_Mag_Moment = 0.d0
        do i = 1, N_Spins
         if (Hmag_Dir.eq.1)
     &    Mag_Moment = S_Value(i) * dsin(Rad * Conf(i, 1)) *
     &                              dcos(Rad * Conf(i, 2))
         if (Hmag_Dir.eq.2)
     &    Mag_Moment = S_Value(i) * dsin(Rad * Conf(i, 1)) *
     &                              dsin(Rad * Conf(i, 2))
         if (Hmag_Dir.eq.3)
     &    Mag_Moment = S_Value(i) * dcos(Rad * Conf(i, 1))
     
         Total_Mag_Moment = Total_Mag_Moment + Mag_Moment
     
         write(6, '(i8, f12.1, f14.1, f17.6)') i,
     &                                         Conf(i, 1), Conf(i, 2),
     &                                         Mag_Moment
         write(7, '(i8, f12.1, f14.1, f17.6)') i,
     &                                         Conf(i, 1), Conf(i, 2),
     &                                         Mag_Moment
        enddo !i
        write(6, *) '-----------------------------
     &------------------------------'
        write(6, '("                  Total magnetic moment:  ", f9.6)')
     &                                Total_Mag_Moment
        write(6, *) '-----------------------------
     &------------------------------'
        close(6)
        close(7)
       
        write(8, '("State No.", i8, "  |  Energy = ", f20.6)')
     &              N_States, Energy
     
       endif !Only if it is the lowest one or Lessened mode is off
! ++++++++++++++++++++++++++++++++++++++++++++++++
      endif !Only new one is to be treated
! --------------------------------------------
      end
! =======================================================================
      subroutine Search_for_the_Energy_minimum
     &         ( n,
     &           N_Spins, N_Couples,
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin, Fitting_Spin_NumByID,
     &           S_Value, Exchange, Couple, Conf,
     &           Hmag_Dir, Hmag, Energy )
      implicit none
! ------------------------------
      integer n
! ------------------------------
      integer N_Spins, N_Couples
      integer Spin_Type(N_Spins)
      integer N_Fitting_Spins
      integer Fitting_Spin(N_Fitting_Spins),
     &        Fitting_Spin_NumByID(N_Spins)
      double precision S_Value(N_Spins)
      double precision Exchange(N_Couples)
      integer Couple(N_Couples, 2)
      double precision Conf(N_Spins, 2)
      integer Hmag_Dir
      double precision Hmag
      double precision Energy
! ------------------------------
      external funct
      integer mode, maxfn, iprint, iexit
      double precision x(n), f, g(n), h(n*(n+1)/2), w(3*n), dfn,
     &                 xm(n), hh, eps
! ------------------------------
      integer i
! ------------------------------


! ----- Set the parameters array ------
      do i = 1, N_Fitting_Spins
!    -- Longitudinal --
       x(i) = Conf(Fitting_Spin(i), 1)
!    --  Transverse  --
       x(i + N_Fitting_Spins) = Conf(Fitting_Spin(i), 2)
      enddo !i
! -------------------------------------
      
! ---- Set the all other VA10AD initial parameters ---
      f = 0.d0
      g = 0.d0
      h = 0.d0
      w = 0.d0
      dfn = 0.d0
      xm = 0.001d0
      hh = 1.d-6
      eps = 0.001d0
      mode = 1
      maxfn = 1000000
      iprint = 0
      iexit = 0
! ---- End of setting ----------------------
      
! ----- Perform the minimization ------      
      call VA10AD
     &   ( funct, n, x, f, g, h, w, dfn, xm,
     &     hh, eps, mode, maxfn, iprint, iexit,
     &     N_Spins, N_Couples, 
     &     Spin_Type,
     &     N_Fitting_Spins,
     &     Fitting_Spin_NumByID,
     &     Conf, S_Value, Exchange, Couple,
     &     Hmag_Dir, Hmag )
! -------------------------------------

! ----- Slice and dice ------     
      Energy = f
      
      do i = 1, N_Fitting_Spins
       Conf(Fitting_Spin(i), 1) = x(i)
       Conf(Fitting_Spin(i), 2) = x(i + N_Fitting_Spins)
      enddo !i
! ---------------------------
      end
! =======================================================================
      subroutine funct
     &         ( n, x, f,
     &           N_Spins, N_Couples, 
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin_NumByID,
     &           Conf, S_Value, Exchange, Couple,
     &           Hmag_Dir, Hmag )
      implicit none
! ------------------------------------
      integer n
      double precision x(n), f
! ------------------------------------
      integer N_Spins, N_Couples
      integer Spin_Type(N_Spins)
      integer N_Fitting_Spins
      integer Fitting_Spin_NumByID(N_Spins)
      double precision Conf(N_Spins, 2)
      double precision S_Value(N_Spins)
      double precision Exchange(N_Couples)
      integer Couple(N_Couples, 2)
      integer Hmag_Dir
      double precision Hmag
! ------------------------------------
      integer Sp1, Sp2
      double precision Theta1, Phi1,
     &                 Theta2, Phi2
      double precision S1_x, S1_y, S1_z,
     &                 S2_x, S2_y, S2_z
     
      double precision Theta, Phi 
! ------------------------------------
      integer i
! ------------------------------------

      f = 0.d0

! +++ Exchange interactions +++++++   
      do i = 1, N_Couples
       Sp1 = Couple(i, 1)
       Sp2 = Couple(i, 2)
       
! ---- Angle definition for the first spin -----
       if (Spin_Type(Sp1).eq.1) then !Fitting case
        Theta1 = x(Fitting_Spin_NumByID(Sp1))
        Phi1   = x(Fitting_Spin_NumByID(Sp1) + N_Fitting_Spins)
       else !Frozen case
        Theta1 = Conf(Sp1, 1) ! Longitudinal
        Phi1   = Conf(Sp1, 2) ! Transverse
       endif !Kind of spin
! ---- Angle definition for the second spin ----
       if (Spin_Type(Sp2).eq.1) then !Fitting case
        Theta2 = x(Fitting_Spin_NumByID(Sp2))
        Phi2   = x(Fitting_Spin_NumByID(Sp2) + N_Fitting_Spins)
       else !Frozen case
        Theta2 = Conf(Sp2, 1) ! Longitudinal
        Phi2   = Conf(Sp2, 2) ! Transverse
       endif !Kind of spin
       
       S1_x = S_Value(Sp1) * dsin(Theta1) * dcos(Phi1)
       S1_y = S_Value(Sp1) * dsin(Theta1) * dsin(Phi1)
       S1_z = S_Value(Sp1) * dcos(Theta1)
       
       S2_x = S_Value(Sp2) * dsin(Theta2) * dcos(Phi2)
       S2_y = S_Value(Sp2) * dsin(Theta2) * dsin(Phi2)
       S2_z = S_Value(Sp2) * dcos(Theta2)
      
       f = f + Exchange(i) * ( S1_x * S2_x +
     &                         S1_y * S2_y +
     &                         S1_z * S2_z )
      enddo !i
! ++++++++++++++++++++++++++++++++

! +++ Magnetic field +++++++++++++
      do i = 1, N_Spins
       if (Spin_Type(i).eq.1) then !Fitting case
        Theta = x(Fitting_Spin_NumByID(i))
        Phi   = x(Fitting_Spin_NumByID(i) + N_Fitting_Spins)
       else !Frozen case
        Theta = Conf(i, 1) ! Longitudinal
        Phi   = Conf(i, 2) ! Transverse
       endif !Kind of spin
      
       if (Hmag_Dir.eq.1) then
        f = f + Hmag * S_Value(i) * dsin(Theta) * dcos(Phi)
       endif !X
       
       if (Hmag_Dir.eq.2) then
        f = f + Hmag * S_Value(i) * dsin(Theta) * dsin(Phi)
       endif !Y
       
       if (Hmag_Dir.eq.3) then
        f = f + Hmag * S_Value(i) * dcos(Theta)
       endif !Z
      enddo !i
! ++++++++++++++++++++++++++++++++

      f = -f ! Main model minus
      end
! =======================================================================
      subroutine VA10AD
     &         ( funct, n, x, f, g, h, w, dfn, xm,
     &           hh, eps, mode, maxfn, iprint, iexit,
     &           N_Spins, N_Couples, 
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin_NumByID,
     &           Conf, S_Value, Exchange, Couple,
     &           Hmag_Dir, Hmag )
c     SUBROUTINE VA10AD(FUNCT,N,X,F,G,H,W,DFN,XM,HH,EPS,MODE,MAXFN,
c    +                  IPRINT,IEXIT)

! ------------------------------------
      integer N_Spins, N_Couples
      integer Spin_Type(N_Spins)
      integer N_Fitting_Spins
      integer Fitting_Spin_NumByID(N_Spins)
      double precision Conf(N_Spins, 2)
      double precision S_Value(N_Spins)
      double precision Exchange(N_Couples)
      integer Couple(N_Couples, 2)
      integer Hmag_Dir
      double precision Hmag
! ------------------------------------

      DOUBLE PRECISION DFN,EPS,F,HH
      INTEGER IEXIT,IPRINT,MAXFN,MODE,N
      DOUBLE PRECISION G(*),H(*),W(*),X(*),XM(*)
      EXTERNAL FUNCT
      DOUBLE PRECISION AEPS,ALPHA,DF,DGS,EPSMCH,F1,F2,FF,GS0,GYS,SIG,
     +                 TOT,Z,ZZ
      INTEGER I,IDIFF,IFN,IG,IGG,IJ,INT,IR,IS,ITN,J,LINK,NN
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
      EXTERNAL MC11AD,MC11BD,MC11ED
      INTRINSIC DABS,MOD
      EPSMCH = FD05AD(1)*10.0D0
      IF (IPRINT.NE.0) WRITE (6,FMT=1000)
 1000 FORMAT ('1ENTRY TO VA10AD',/)
      NN = N* (N+1)/2
      IG = N
      IGG = N + N
      IS = IGG
      IDIFF = 1
      IEXIT = 0
      IR = N
      IF (MODE.EQ.3) GO TO 15
      IF (MODE.EQ.2) GO TO 10
      IJ = NN + 1
      DO 5 I = 1,N
        DO 6 J = 1,I
          IJ = IJ - 1
    6   H(IJ) = 0.D0
    5 H(IJ) = 1.D0
      GO TO 15
   10 CONTINUE
      CALL MC11BD(H,N,IR)
      IF (IR.LT.N) RETURN
   15 CONTINUE
      Z = F
      ITN = 0
      CALL FUNCT(N,X,F,
     &           N_Spins, N_Couples, 
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin_NumByID,
     &           Conf, S_Value, Exchange, Couple,
     &           Hmag_Dir, Hmag)
      IFN = 1
      DF = DFN
      IF (DFN.EQ.0D0) DF = F - Z
      IF (DFN.LT.0D0) DF = DABS(DF*F)
      IF (DF.LE.0D0) DF = 1.D0
   17 CONTINUE
      LINK = 1
      IF (IDIFF-1) 100,100,110
   18 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
   20 CONTINUE
      IF (IPRINT.EQ.0) GO TO 21
      IF (MOD(ITN,IPRINT).NE.0) GO TO 21
      WRITE (6,FMT=1001) ITN,IFN
 1001 FORMAT (24I5)
      WRITE (6,FMT=1002) F
 1002 FORMAT (5D24.16)
      IF (IPRINT.LT.0) GO TO 21
      WRITE (6,FMT=1002) (X(I),I=1,N)
      WRITE (6,FMT=1002) (W(IG+I),I=1,N)
   21 CONTINUE
      ITN = ITN + 1
      DO 22 I = 1,N
   22 W(I) = -W(IG+I)
      CALL MC11ED(H,N,W,G,IR)
      Z = 0.D0
      GS0 = 0.D0
      DO 29 I = 1,N
        W(IS+I) = W(I)
        IF (Z*XM(I).GE.DABS(W(I))) GO TO 29
        Z = DABS(W(I))/XM(I)
   29 GS0 = GS0 + W(IG+I)*W(I)
      AEPS = EPS/Z
      IEXIT = 2
      IF (GS0.GE.0D0) GO TO 92
      ALPHA = -2D0*DF/GS0
      IF (ALPHA.GT.1D0) ALPHA = 1D0
      FF = F
      TOT = 0.D0
      INT = 0
      IEXIT = 1
   30 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
      DO 31 I = 1,N
   31 W(I) = X(I) + ALPHA*W(IS+I)
      CALL FUNCT(N,W,F1,
     &           N_Spins, N_Couples, 
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin_NumByID,
     &           Conf, S_Value, Exchange, Couple,
     &           Hmag_Dir, Hmag)
      IFN = IFN + 1
      IF (F1.GE.F) GO TO 40
      F2 = F
      TOT = TOT + ALPHA
   32 CONTINUE
      DO 33 I = 1,N
   33 X(I) = W(I)
      F = F1
      IF (INT-1) 35,49,50
   35 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
      DO 34 I = 1,N
   34 W(I) = X(I) + ALPHA*W(IS+I)
      CALL FUNCT(N,W,F1,
     &           N_Spins, N_Couples, 
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin_NumByID,
     &           Conf, S_Value, Exchange, Couple,
     &           Hmag_Dir, Hmag)
      IFN = IFN + 1
      IF (F1.GE.F) GO TO 50
      IF (F1+F2.GE.F+F .AND. 7D0*F1+5D0*F2.GT.12D0*F) INT = 2
      TOT = TOT + ALPHA
      ALPHA = 2D0*ALPHA
      GO TO 32
   40 CONTINUE
      IF (ALPHA.LT.AEPS) GO TO 92
      IF (IFN.GE.MAXFN) GO TO 90
      ALPHA = .5D0*ALPHA
      DO 41 I = 1,N
   41 W(I) = X(I) + ALPHA*W(IS+I)
      CALL FUNCT(N,W,F2,
     &           N_Spins, N_Couples, 
     &           Spin_Type,
     &           N_Fitting_Spins,
     &           Fitting_Spin_NumByID,
     &           Conf, S_Value, Exchange, Couple,
     &           Hmag_Dir, Hmag)
      IFN = IFN + 1
      IF (F2.GE.F) GO TO 45
      TOT = TOT + ALPHA
      F = F2
      DO 42 I = 1,N
   42 X(I) = W(I)
      GO TO 49
   45 CONTINUE
      Z = .1D0
      IF (F1+F.GT.F2+F2) Z = 1D0 + .5D0* (F-F1)/ (F+F1-F2-F2)
      IF (Z.LT..1D0) Z = .1D0
      ALPHA = Z*ALPHA
      INT = 1
      GO TO 30
   49 CONTINUE
      IF (TOT.LT.AEPS) GO TO 92
   50 CONTINUE
      ALPHA = TOT
      DO 56 I = 1,N
   56 W(I) = W(IG+I)
      LINK = 2
      IF (IDIFF-1) 100,100,110
   54 CONTINUE
      IF (IFN.GE.MAXFN) GO TO 90
      GYS = 0.D0
      DO 55 I = 1,N
        GYS = GYS + W(IG+I)*W(IS+I)
   55 W(IGG+I) = W(I)
      DF = FF - F
      DGS = GYS - GS0
      IF (DGS.LE.0D0) GO TO 20
      IF (DGS+ALPHA*GS0.GT.0D0) GO TO 70
      SIG = 1D0/GS0
      IR = -IR
      CALL MC11AD(H,N,W,SIG,G,IR,1,0D0)
      DO 60 I = 1,N
   60 G(I) = W(IG+I) - W(IGG+I)
      SIG = 1D0/ (ALPHA*DGS)
      IR = -IR
      CALL MC11AD(H,N,G,SIG,W,IR,0,0D0)
      GO TO 20
   70 CONTINUE
      ZZ = ALPHA/ (DGS-ALPHA*GS0)
      SIG = -ZZ
      CALL MC11AD(H,N,W,SIG,G,IR,1,EPSMCH)
      Z = DGS*ZZ - 1.D0
      DO 71 I = 1,N
   71 G(I) = W(IG+I) + Z*W(IGG+I)
      SIG = 1.D0/ (ZZ*DGS**2)
      CALL MC11AD(H,N,G,SIG,W,IR,0,0D0)
      GO TO 20
   90 CONTINUE
      IEXIT = 3
      GO TO 94
   92 CONTINUE
      IF (IDIFF.EQ.2) GO TO 94
      IDIFF = 2
      GO TO 17
   94 CONTINUE
      DO 95 I = 1,N
   95 G(I) = W(IG+I)
      IF (IPRINT.EQ.0) RETURN
      WRITE (6,FMT=1001) ITN,IFN,IEXIT
      WRITE (6,FMT=1002) F
      WRITE (6,FMT=1002) (X(I),I=1,N)
      WRITE (6,FMT=1002) (G(I),I=1,N)
      RETURN
  100 CONTINUE
      DO 101 I = 1,N
        Z = HH*XM(I)
        ZZ = X(I)
        X(I) = ZZ + Z
        CALL FUNCT(N,X,F1,
     &             N_Spins, N_Couples, 
     &             Spin_Type,
     &             N_Fitting_Spins,
     &             Fitting_Spin_NumByID,
     &             Conf, S_Value, Exchange, Couple,
     &             Hmag_Dir, Hmag)
        W(IG+I) = (F1-F)/Z
  101 X(I) = ZZ
      IFN = IFN + N
  102 GO TO (18,54) LINK
  110 CONTINUE
      DO 111 I = 1,N
        Z = HH*XM(I)
        ZZ = X(I)
        X(I) = ZZ + Z
        CALL FUNCT(N,X,F1,
     &             N_Spins, N_Couples, 
     &             Spin_Type,
     &             N_Fitting_Spins,
     &             Fitting_Spin_NumByID,
     &             Conf, S_Value, Exchange, Couple,
     &             Hmag_Dir, Hmag)
        X(I) = ZZ - Z
        CALL FUNCT(N,X,F2,
     &             N_Spins, N_Couples, 
     &             Spin_Type,
     &             N_Fitting_Spins,
     &             Fitting_Spin_NumByID,
     &             Conf, S_Value, Exchange, Couple,
     &             Hmag_Dir, Hmag)
        W(IG+I) = (F1-F2)/ (2D0*Z)
  111 X(I) = ZZ
      IFN = IFN + N + N
      GO TO 102
      END
*######DATE 1 Feb 1993 COPYRIGHT AEA Technology
C       Toolpack tool decs employed.
C       Arg dimensions set to *.
C
      SUBROUTINE MC11AD(A,N,Z,SIG,W,IR,MK,EPS)
      DOUBLE PRECISION EPS,SIG
      INTEGER IR,MK,N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION AL,B,GM,R,TI,TIM,V,Y
      INTEGER I,IJ,IP,J,MM,NP
      IF (N.GT.1) GO TO 1
      A(1) = A(1) + SIG*Z(1)**2
      IR = 1
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
    1 CONTINUE
      NP = N + 1
      IF (SIG.GT.0.0D0) GO TO 40
      IF (SIG.EQ.0.0D0 .OR. IR.EQ.0) RETURN
      TI = 1.0D0/SIG
      IJ = 1
      IF (MK.EQ.0) GO TO 10
      DO 7 I = 1,N
        IF (A(IJ).NE.0.0D0) TI = TI + W(I)**2/A(IJ)
    7 IJ = IJ + NP - I
      GO TO 20
   10 CONTINUE
      DO 11 I = 1,N
   11 W(I) = Z(I)
      DO 15 I = 1,N
        IP = I + 1
        V = W(I)
        IF (A(IJ).GT.0.0D0) GO TO 12
        W(I) = 0.D0
        IJ = IJ + NP - I
        GO TO 15
   12   CONTINUE
        TI = TI + V**2/A(IJ)
        IF (I.EQ.N) GO TO 14
        DO 13 J = IP,N
          IJ = IJ + 1
   13   W(J) = W(J) - V*A(IJ)
   14   IJ = IJ + 1
   15 CONTINUE
   20 CONTINUE
      IF (IR.LE.0) GO TO 21
      IF (TI.GT.0.0D0) GO TO 22
      IF (MK-1) 40,40,23
   21 TI = 0.D0
      IR = -IR - 1
      GO TO 23
   22 TI = EPS/SIG
      IF (EPS.EQ.0.0D0) IR = IR - 1
   23 CONTINUE
      MM = 1
      TIM = TI
      DO 30 I = 1,N
        J = NP - I
        IJ = IJ - I
        IF (A(IJ).NE.0.0D0) TIM = TI - W(J)**2/A(IJ)
        W(J) = TI
   30 TI = TIM
      GO TO 41
   40 CONTINUE
      MM = 0
      TIM = 1.0D0/SIG
   41 CONTINUE
      IJ = 1
      DO 66 I = 1,N
        IP = I + 1
        V = Z(I)
        IF (A(IJ).GT.0.0D0) GO TO 53
        IF (IR.GT.0 .OR. SIG.LT.0.0D0 .OR. V.EQ.0.0D0) GO TO 52
        IR = 1 - IR
        A(IJ) = V**2/TIM
        IF (I.EQ.N) RETURN
        DO 51 J = IP,N
          IJ = IJ + 1
   51   A(IJ) = Z(J)/V
        RETURN
   52   CONTINUE
        TI = TIM
        IJ = IJ + NP - I
        GO TO 66
   53   CONTINUE
        AL = V/A(IJ)
        IF (MM) 54,54,55
   54   TI = TIM + V*AL
        GO TO 56
   55   TI = W(I)
   56   CONTINUE
        R = TI/TIM
        A(IJ) = A(IJ)*R
        IF (R.EQ.0.0D0) GO TO 70
        IF (I.EQ.N) GO TO 70
        B = AL/TI
        IF (R.GT.4.0D0) GO TO 62
        DO 61 J = IP,N
          IJ = IJ + 1
          Z(J) = Z(J) - V*A(IJ)
   61   A(IJ) = A(IJ) + B*Z(J)
        GO TO 64
   62   GM = TIM/TI
        DO 63 J = IP,N
          IJ = IJ + 1
          Y = A(IJ)
          A(IJ) = B*Z(J) + Y*GM
   63   Z(J) = Z(J) - V*Y
   64   CONTINUE
        TIM = TI
        IJ = IJ + 1
   66 CONTINUE
   70 CONTINUE
      IF (IR.LT.0) IR = -IR
      RETURN
      END
      SUBROUTINE MC11BD(A,N,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER I,II,IJ,IK,IP,JK,NI,NP
      IR = N
      IF (N.GT.1) GO TO 100
      IF (A(1).GT.0.0D0) RETURN
      A(1) = 0.D0
      IR = 0
      RETURN
  100 CONTINUE
      NP = N + 1
      II = 1
      DO 104 I = 2,N
        AA = A(II)
        NI = II + NP - I
        IF (AA.GT.0.0D0) GO TO 101
        A(II) = 0.D0
        IR = IR - 1
        II = NI + 1
        GO TO 104
  101   CONTINUE
        IP = II + 1
        II = NI + 1
        JK = II
        DO 103 IJ = IP,NI
          V = A(IJ)/AA
          DO 102 IK = IJ,NI
            A(JK) = A(JK) - A(IK)*V
  102     JK = JK + 1
  103   A(IJ) = V
  104 CONTINUE
      IF (A(II).GT.0.0D0) RETURN
      A(II) = 0.D0
      IR = IR - 1
      RETURN
      END
      SUBROUTINE MC11CD(A,N)
      INTEGER N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER II,IJ,IK,IP,JK,NI,NIP,NP
      IF (N.EQ.1) RETURN
      NP = N + 1
      II = N*NP/2
      DO 202 NIP = 2,N
        JK = II
        NI = II - 1
        II = II - NIP
        AA = A(II)
        IP = II + 1
        IF (AA.GT.0.0D0) GO TO 203
        DO 204 IJ = IP,NI
  204   A(IJ) = 0.D0
        GO TO 202
  203   CONTINUE
        DO 201 IJ = IP,NI
          V = A(IJ)*AA
          DO 200 IK = IJ,NI
            A(JK) = A(JK) + A(IK)*V
  200     JK = JK + 1
  201   A(IJ) = V
  202 CONTINUE
      RETURN
      END
      SUBROUTINE MC11DD(A,N,Z,W)
      INTEGER N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION Y
      INTEGER I,II,IJ,IP,J,K,N1,NP
      IF (N.GT.1) GO TO 300
      Z(1) = Z(1)*A(1)
      W(1) = Z(1)
      RETURN
  300 CONTINUE
      NP = N + 1
      II = 1
      N1 = N - 1
      DO 303 I = 1,N1
        Y = Z(I)
        IF (A(II).EQ.0.0D0) GO TO 302
        IJ = II
        IP = I + 1
        DO 301 J = IP,N
          IJ = IJ + 1
  301   Y = Y + Z(J)*A(IJ)
  302   Z(I) = Y*A(II)
        W(I) = Z(I)
  303 II = II + NP - I
      Z(N) = Z(N)*A(II)
      W(N) = Z(N)
      DO 311 K = 1,N1
        I = N - K
        II = II - NP + I
        IF (Z(I).EQ.0.0D0) GO TO 311
        IP = I + 1
        IJ = II
        Y = Z(I)
        DO 310 J = IP,N
          IJ = IJ + 1
  310   Z(J) = Z(J) + A(IJ)*Z(I)
  311 CONTINUE
      RETURN
      END
      SUBROUTINE MC11ED(A,N,Z,W,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*),W(*),Z(*)
      DOUBLE PRECISION V
      INTEGER I,I1,II,IJ,IP,J,NIP,NP
      IF (IR.LT.N) RETURN
      W(1) = Z(1)
      IF (N.GT.1) GO TO 400
      Z(1) = Z(1)/A(1)
      RETURN
  400 CONTINUE
      DO 402 I = 2,N
        IJ = I
        I1 = I - 1
        V = Z(I)
        DO 401 J = 1,I1
          V = V - A(IJ)*Z(J)
  401   IJ = IJ + N - J
        W(I) = V
  402 Z(I) = V
      Z(N) = Z(N)/A(IJ)
      NP = N + 1
      DO 411 NIP = 2,N
        I = NP - NIP
        II = IJ - NIP
        V = Z(I)/A(II)
        IP = I + 1
        IJ = II
        DO 410 J = IP,N
          II = II + 1
  410   V = V - A(II)*Z(J)
  411 Z(I) = V
      RETURN
      END
      SUBROUTINE MC11FD(A,N,IR)
      INTEGER IR,N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION AA,V
      INTEGER I,I1,II,IJ,IK,IP,J,JK,K,N1,NI,NIP,NP
      IF (IR.LT.N) RETURN
      A(1) = 1.0D0/A(1)
      IF (N.EQ.1) RETURN
      NP = N + 1
      N1 = N - 1
      II = 2
      DO 511 I = 2,N
        A(II) = -A(II)
        IJ = II + 1
        IF (I.EQ.N) GO TO 502
        DO 501 J = I,N1
          IK = II
          JK = IJ
          V = A(IJ)
          DO 500 K = I,J
            JK = JK + NP - K
            V = V + A(IK)*A(JK)
  500     IK = IK + 1
          A(IJ) = -V
  501   IJ = IJ + 1
  502   CONTINUE
        A(IJ) = 1.0D0/A(IJ)
        II = IJ + 1
        AA = A(IJ)
        IJ = I
        IP = I + 1
        NI = N - I
        DO 511 J = 2,I
          V = A(IJ)*AA
          IK = IJ
          K = IJ - IP + J
          I1 = IJ - 1
          NIP = NI + IJ
          DO 510 JK = K,I1
            A(JK) = A(JK) + V*A(IK)
  510     IK = IK + NIP - JK
          A(IJ) = V
  511 IJ = IJ + NP - J
      RETURN
      END
C#####################################################
      DOUBLE PRECISION FUNCTION FD05AD( INUM )
      INTEGER INUM
      DOUBLE PRECISION DC( 5 )
C
C  REAL CONSTANTS (DOUBLE PRECISION ARITHMETIC).
C
C  OBTAINED FROM H.S.L. SUBROUTINE ZE02AM.
C  NICK GOULD AND SID MARLOW, HARWELL, JULY 1988.
C
C  DC(1) THE 'SMALLEST' POSITIVE NUMBER: 1 + DC(1) > 1.
C  DC(2) THE 'SMALLEST' POSITIVE NUMBER: 1 - DC(2) < 1.
C  DC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
C  DC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
C  DC(5) THE LARGEST FINITE +VE REAL NUMBER.
C
      SAVE DC
      DATA DC( 1 ) /    2.220446049250314D-016 /
      DATA DC( 2 ) /    1.110223024625158D-016 /
      DATA DC( 3 ) /    2.225073858507202D-308 /
      DATA DC( 4 ) /    2.225073858507202D-308 /
      DATA DC( 5 ) /    1.797693134862315D+308 /
      IF ( INUM .LE. 0 .OR. INUM .GE. 6 ) THEN
         WRITE(6, 2000) INUM
         STOP
      ELSE
         FD05AD = DC( INUM )
      ENDIF
      RETURN
 2000 FORMAT( ' INUM =', I3, ' OUT OF RANGE IN FD05AD.',
     *        ' EXECUTION TERMINATED.' )
      END
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: tql2
C TYPE   : subroutine
C PURPOSE: Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix
C I/O    : 
C VERSION:
C COMMENT: Minimally modified SLATEC routine from netlib.
C          To learn about netlib, send an otherwise empty
C          message to netlib@research.att.com
C          containing 'send index' in the subject header)
C          The WWW address of netlib is
C          http://netlib.att.com/netlib/search.html
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
*DECK TQL2
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
C***BEGIN PROLOGUE  TQL2
C***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (TQL2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional DOUBLE PRECISION array, 
C          dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional DOUBLE PRECISION array, 
C          dimensioned E(N).
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.  Z is a two-dimensional DOUBLE PRECISION array,
C          dimensioned Z(NM,N).
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1, 2, ..., IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  TQL2
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(*),E(*),Z(NM,*)
      DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      DOUBLE PRECISION PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: pythag 
C TYPE   : function
C PURPOSE: compute dsqrt(a**2+b**2) without overflow or
C          destructive underflow
C I/O    : 
C VERSION:
C COMMENT: 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      double precision function pythag(a,b)
      double precision a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
! =======================================================================
      subroutine plant_random_seed
      implicit none
! ------------------------
      integer n
      integer, allocatable :: seed(:)
! ------------------------
      integer ios
      integer(4) buf
! ------------------------
      integer i
! ------------------------

      call random_seed(size = n)
      allocate(seed(n))
      
! ----- Tricky seed planter -----
      open( unit = 777,
     &      file = '/dev/urandom',
     &      access = 'stream',
     &      form = 'UNFORMATTED',
     &      iostat = ios,
     &      status = 'old',
     &      action = 'read' )
      if (ios.ne.0) stop 'Error opening random device file'
      do i = 1, n
       read(777) buf
       seed(i) = buf
      enddo !i
      close(777)
! -------------------------------
      
      call random_seed(put = seed)
      deallocate(seed)
      
      end
! =======================================================================
