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
