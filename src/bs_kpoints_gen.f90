  program bs_kpoints_gen 
     use highsymmetrykpoints
     implicit none 

!     call gen_kpoints4bs()
      call readkpoints()
      call Dealloc_high_sym_kpoints()
       
    end program
