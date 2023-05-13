subroutine f2c_calc_energy_force(i_model_lvn, &
                                     inatoms, &
                                       inall, &
                                     intypes, &
                                     icatype, &
                                      iatype, &
                                  inum_neigh, &
                                 ilist_neigh, &
                                   idR_neigh, &
                                  e_atom_out, &
                                   e_tot_out, &
                                  f_atom_out, &
                                  virial_out, &
                                      ff_idx ) bind(c,name="f2c_calc_energy_force")

    !**************************************************
    !  Inference routine run on a SINGLE thread
    !  with neighbor list AND network paras passed-in 
    !  L. Wang 2023.1 
    !**************************************************
    use iso_c_binding
    use mod_data  ! natoms, nall, ntypes, catype, atype, iflag_model
    use calc_deepMD, only: load_model_deepMD, set_image_info_deepMD
    use calc_deepMD_f,only: load_model_deepMD_f,set_image_info_deepMD_f
    
    implicit none

    integer, intent(in) :: i_model_lvn, inatoms, inall, intypes
    integer, dimension(inall), intent(in) :: icatype
    integer, dimension(inall), intent(in) :: iatype
    integer, dimension(intypes, inatoms), intent(in) :: inum_neigh
    integer, dimension(100, intypes, inatoms), intent(in) :: ilist_neigh
    real*8, dimension(3, 100, intypes, inatoms), intent(in) :: idR_neigh

    !! $ spend almost two days to change inall to inatoms.
    !! $ inum_neigh value reset to 0.
    real*8, dimension(inatoms), intent(out) :: e_atom_out
    real*8, intent(out) :: e_tot_out
    real*8, dimension(3,inall), intent(out) :: f_atom_out
    real*8, dimension(6), intent(out) :: virial_out
    integer, intent(in) :: ff_idx

    logical :: atype_check
    integer :: i, j, t

    !write(*,*) " --- f2c CHECK --- "
    !write(*,*) "imodel ", i_model_lvn
    !write(*,*) "nlocal ", inatoms
    !write(*,*) "nall   ", inall
    !write(*,*) "ntypes ", intypes
    !write(*,*) "ff_idx ", ff_idx

    !!!do i = 1, inall
    !!!  write(*,*) i, icatype(i), iatype(i)
    !!!enddo

    !do i = 1, inatoms
    !  write(*,*) i
    !  do t = 1, intypes
    !    write(*,*) "  ", t, inum_neigh(t,i)
    !    !do j = 1, inum_neigh(t,i)
    !    !  write(*,"('     ', A, 2(I5, 1X), 3(F12.6, 3X))") &
    !    !            "     ", j, ilist_neigh(j,t,i), &
    !    !    idR_neigh(1,j,t,i), idR_neigh(2,j,t,i), idR_neigh(3,j,t,i)
    !    !enddo
    !  enddo
    !enddo

    natoms = inatoms
    nall = inall
    ntypes = intypes
    catype(1:nall) = icatype(1:nall)   ! lammps atom type
    atype(1:nall) = iatype(1:nall)     ! periodic table index
    iflag_model = i_model_lvn
        
    ! ********************************************
    !   Initialization: load control & net paras 
    ! ********************************************

    if(iflag_model.eq.5) then

        ! load feature cutoff, shift and norm
        call load_model_deepMD_f(ff_idx)
        ! allocate memroy for features
        call set_image_info_deepMD_f(ff_idx)

        ! load the network parameters
        call load_model_deepMD(ff_idx)
        ! check atom type
        atype_check = .True.
        call set_image_info_deepMD(atype_check)

    endif

    ! ***************************************
    !           FF calculation
    ! *************************************** 

    e_atom_out = 0.0
    e_tot_out = 0.0
    f_atom_out = 0.0
    virial_out = 0.0

    call ML_FF_EF(inum_neigh, ilist_neigh, idR_neigh, &
           e_tot_out, f_atom_out, virial_out)

    ! a) the force calculated by ML_FF_EF is \partial E / \partial x,
    !    lacking a minus sign.
    ! b) the same as virial.
    f_atom_out(1:3,1:nall) = -1.0*f_atom_out(1:3,1:nall)
    virial_out(1:6) = -1.0*virial_out(1:6)

end


