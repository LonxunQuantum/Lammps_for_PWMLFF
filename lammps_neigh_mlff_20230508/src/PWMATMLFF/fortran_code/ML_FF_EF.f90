subroutine ML_FF_EF(num_neigh,list_neigh,dR_neigh,Etot,fatom,virial)
        use IFPORT
        use mod_data, only : natoms, nall, ntypes, iflag_model    
        use calc_deepMD, only : m_neigh, cal_energy_force_deepMD        

        implicit none

        integer, dimension(ntypes,natoms), intent(in) :: num_neigh
        integer, dimension(m_neigh,ntypes,natoms), intent(in) :: list_neigh
        real*8, dimension(3,m_neigh,ntypes,natoms), intent(in) :: dR_neigh
        real*8, intent(out) :: Etot
        real*8, dimension(3, nall), intent(out) :: fatom 
        real*8, dimension(6), intent(out) :: virial
        integer :: i, j, t

        !write(*,*) "address 2 ", LOC(num_neigh)    
        !write(*,*) "address 2 ", LOC(list_neigh)    

        !write(*,*) " --- ML CHECK --- "
        !write(*,*) "iflag ", iflag_model, m_neigh
        !write(*,*) "ntypes ", ntypes, natoms
        !do i = 1, natoms
        !  write(*,*) i
        !  do t = 1, ntypes
        !    write(*,*) " ", t, num_neigh(t,i)
        !    do j = 1, 6
        !      write(*,*) " l ", list_neigh(j,t,i), dR_neigh(1,j,t,i)
        !    enddo
        !  enddo
        !enddo

        ! ***************************************
        !              inference
        ! ***************************************
        if(iflag_model.eq.5) then
          ! DP model  
          call cal_energy_force_deepMD(num_neigh,list_neigh,dR_neigh,Etot,fatom,virial)
        endif
               
        return

end subroutine ML_FF_EF
        

