module li12_ff_mod

   ! **************************************************************
   !           Overarching data module for forece field
   !   Will be called from LAMMPS pair constructor and destructor
   ! **************************************************************
   use iso_c_binding
   implicit none

   ! wrap up a ff with an struct
   type single_ff
      ! local temporary variables
      !character(200) ff_name

      ! force field paras
      integer(4) :: li_ff_model_idx

      ! feat info part
      integer(4) :: li_ff_iflag_PCA
      integer(4) :: li_ff_nfeat_type     ! number of feat type
      integer(4) :: li_ff_ifeat_type(10)        ! feat type list
      integer(4) :: li_ff_num_type
      ! integer(4) :: li_ff_itype_atom_sumfe(10)        ! a list with atomtype and sum of feats
      integer(4),allocatable,dimension(:,:) :: li_ff_itype_atom_sumfe
      integer(4) :: li_ff_itype_atom_perfe(10)        ! a list with per feats of each atomtype

      ! gen_xxx_feature.in part
      real*8 li_ff_Rc_M                  ! cutoff
      integer(4) li_ff_max_neigh         ! max nieghbor num
      integer li_ff_iat_type(100)              ! iat_type
      real*8 li_ff_Rc_type(100),li_ff_Rm_type(100),li_ff_fact_grid_type(100),li_ff_dR_grid1_type(100)
      integer li_ff_iflag_grid_type(100)
      integer li_ff_n2b_type(100)                 ! numOf2bfeat
      real*8 li_ff_E_tolerance
      integer li_ff_iflag_ftype,li_ff_recalc_grid
      real*8 li_ff_Rc2_type(100),li_ff_dR_grid2_type(100)
      integer li_ff_n3b1_type(100),li_ff_n3b2_type(100)

      ! grid part
      ! real*8,allocatable,dimension (:,:) :: li_ff_grid2
      real*8,allocatable,dimension (:,:,:) :: li_ff_grid2_2,li_ff_grid31_2,li_ff_grid32_2
      integer li_ff_iat_feat,temp1,temp2,temp3,n2b_tmp,n3b1_tmp,n3b2_tmp
      integer n2b_tmp_idx,n3b1_tmp_idx,n3b2_tmp_idx
      ! integer li_ff_2bgrid_flag,li_ff_3bgrid_flag_b1b2,li_ff_3bgrid_flag_cb12

      ! shift and scale
      real*8,allocatable,dimension (:,:) :: li_ff_feat2_shift,li_ff_feat2_scale
      ! weight part
      real*8,allocatable,dimension(:,:) :: w_feat

      ! fit coeff part
      real*8,allocatable,dimension(:) :: BB
      integer :: nfeat2tot              ! nfeature1 + nfeature2
      integer :: ifeat2tot              ! idx for nfeat2tot

      ! vdw part
      integer(4) :: nterm
      real*8,allocatable,dimension(:) :: rad_atom,E_ave_vdw
      real*8,allocatable,dimension(:,:,:) :: wp_atom

   end type single_ff

   ! structs for ff data
   type(single_ff) ff(4)

   ! loop index
   integer(4) i,j,jj,k,kkk
   integer(4) itype,itype1,itype2,num,k1,k2,k12,ii_f

   ! others tmps
   logical(2) alive_li
   integer n2bm      ! numOf2bfeat
   integer n3b1m,n3b2m
   ! integer nfeat0m_1, nfeat0m_2   ! _1 for feature1, _2 for feature2
   integer iflag_grid,iflag_ftype
   integer :: tmp_nfeat2
   integer(4) :: ntype_vdw,itype_vdw


contains

   subroutine li_ff_load(name_ptr, ff_idx, slen) bind(c,name="li_ff_load")
      character, dimension(300), intent(in) :: name_ptr    ! name string pointer
      integer, intent(in) :: slen            ! length of name string
      integer, intent(in) :: ff_idx          ! index of ff to be loaded
      ! real*8, intent(out) :: ocut
      character(300) temp

      !write(*,*) "name_ptr(1) ", name_ptr(1)
      temp = trim(name_ptr(1))
      !write(*,*) "start", temp
      do i=2,slen
         temp = trim(temp)//trim(name_ptr(i))
      enddo

      inquire(file=temp,exist=alive_li)
      if (alive_li.ne..true.) then
         write(*,*) "force field named ", temp, " not found. Terminate."
         stop
      endif

      ! *********************************************
      !      Allocate the ff pointers and load
      ! *********************************************

      ! reading feat.info part
      open(10, file=temp)
      read(10,*) ff(ff_idx)%li_ff_model_idx
      if (ff(ff_idx)%li_ff_model_idx.ne.1) then
         write(*,*) "Model type error. Should be 1"
         stop
      endif

      read(10,*) ff(ff_idx)%li_ff_iflag_PCA   ! this can be used to turn off degmm part
      read(10,*) ff(ff_idx)%li_ff_nfeat_type  ! number of feat type
      do i=1,ff(ff_idx)%li_ff_nfeat_type
         read(10,*) ff(ff_idx)%li_ff_ifeat_type(i)
      enddo
      read(10,*) ff(ff_idx)%li_ff_num_type
      allocate(ff(ff_idx)%li_ff_itype_atom_sumfe(ff(ff_idx)%li_ff_num_type,3))
      tmp_nfeat2=0
      do i=1,ff(ff_idx)%li_ff_num_type
         ! read(10,*) (ff(ff_idx)%li_ff_itype_atom_sumfe(k),k=1,ff(ff_idx)%li_ff_num_type)
         read(10,*) (ff(ff_idx)%li_ff_itype_atom_sumfe(i,k),k=1,3)
         if (ff(ff_idx)%li_ff_itype_atom_sumfe(i,3).gt.tmp_nfeat2) then
            tmp_nfeat2=ff(ff_idx)%li_ff_itype_atom_sumfe(i,3)
         endif
         ! write(*,*) ff(ff_idx)%li_ff_itype_atom_sumfe(i,1),ff(ff_idx)%li_ff_itype_atom_sumfe(i,2),ff(ff_idx)%li_ff_itype_atom_sumfe(i,3)
         ! read(10,*) iatom_tmp,nfeat1(i),nfeat2(i)
      enddo
      do i=1,ff(ff_idx)%li_ff_num_type
         read(10,*) (ff(ff_idx)%li_ff_itype_atom_perfe(k), k=1,ff(ff_idx)%li_ff_nfeat_type)
      enddo

      ! reading gen_2b_feature.in part
      read(10,*) ff(ff_idx)%li_ff_Rc_M, ff(ff_idx)%li_ff_max_neigh
      read(10,*) ff(ff_idx)%li_ff_num_type

      do i=1,ff(ff_idx)%li_ff_num_type
         read(10,*) ff(ff_idx)%li_ff_iat_type(i)
         read(10,*) ff(ff_idx)%li_ff_Rc_type(i), &
            ff(ff_idx)%li_ff_Rm_type(i), &
            ff(ff_idx)%li_ff_iflag_grid_type(i), &
            ff(ff_idx)%li_ff_fact_grid_type(i), &
            ff(ff_idx)%li_ff_dR_grid1_type(i)
         read(10,*) ff(ff_idx)%li_ff_n2b_type(i)

         if(ff(ff_idx)%li_ff_Rc_type(i).gt.ff(ff_idx)%li_ff_Rc_M) then
            write(6,*) "Rc_type must be smaller than Rc_M, gen_3b_feature.in", &
               i,ff(ff_idx)%li_ff_Rc_type(i),ff(ff_idx)%li_ff_Rc_M
            stop
         endif
      enddo

      read(10,*) ff(ff_idx)%li_ff_E_tolerance
      read(10,*) ff(ff_idx)%li_ff_iflag_ftype
      read(10,*) ff(ff_idx)%li_ff_recalc_grid

      ! reading gen_3b_feature.in part
      read(10,*) ff(ff_idx)%li_ff_Rc_M, ff(ff_idx)%li_ff_max_neigh
      read(10,*) ff(ff_idx)%li_ff_num_type

      do i=1,ff(ff_idx)%li_ff_num_type
         read(10,*) ff(ff_idx)%li_ff_iat_type(i)
         read(10,*) ff(ff_idx)%li_ff_Rc_type(i), &
            ff(ff_idx)%li_ff_Rc2_type(i), &
            ff(ff_idx)%li_ff_Rm_type(i), &
            ff(ff_idx)%li_ff_iflag_grid_type(i), &
            ff(ff_idx)%li_ff_fact_grid_type(i), &
            ff(ff_idx)%li_ff_dR_grid1_type(i), &
            ff(ff_idx)%li_ff_dR_grid2_type(i)
         read(10,*) ff(ff_idx)%li_ff_n3b1_type(i), ff(ff_idx)%li_ff_n3b2_type(i)

         if(ff(ff_idx)%li_ff_Rc_type(i).gt.ff(ff_idx)%li_ff_Rc_M) then
            write(6,*) "Rc_type must be smaller than Rc_M, gen_3b_feature.in", &
               i,ff(ff_idx)%li_ff_Rc_type(i),ff(ff_idx)%li_ff_Rc_M
            stop
         endif

         if(ff(ff_idx)%li_ff_Rc2_type(i).gt.2*ff(ff_idx)%li_ff_Rc_type(i)) then
            write(6,*) "Rc2_type must be smaller than 2*Rc_type, gen_3b_feature.in", &
               i,ff(ff_idx)%li_ff_Rc_type(i),ff(ff_idx)%li_ff_Rc2_type(i)
            stop
         endif
      enddo

      read(10,*) ff(ff_idx)%li_ff_E_tolerance
      read(10,*) ff(ff_idx)%li_ff_iflag_ftype
      read(10,*) ff(ff_idx)%li_ff_recalc_grid

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! calculate features of all types
      n2bm=0
      n3b1m=0
      n3b2m=0

      iflag_ftype=ff(ff_idx)%li_ff_iflag_ftype

      do i=1,ff(ff_idx)%li_ff_num_type
         if(iflag_ftype.eq.3.and.ff(ff_idx)%li_ff_iflag_grid_type(i).ne.3) then
            write(6,*) "if iflag_ftype.eq.3, iflag_grid must equal 3, stop"
            stop
         endif
         ! for 2b
         if(ff(ff_idx)%li_ff_n2b_type(i).gt.n2bm) then
            n2bm=ff(ff_idx)%li_ff_n2b_type(i)
         endif
         ! for 3b
         if(ff(ff_idx)%li_ff_n3b1_type(i).gt.n3b1m) then
            n3b1m=ff(ff_idx)%li_ff_n3b1_type(i)
         endif
         if(ff(ff_idx)%li_ff_n3b2_type(i).gt.n3b2m) then
            n3b2m=ff(ff_idx)%li_ff_n3b2_type(i)
         endif
      enddo

      ! nfeat0m_1=ff(ff_idx)%li_ff_num_type*n2bm

      num=0
      do itype2=1,ff(ff_idx)%li_ff_num_type
         do itype1=1,itype2
            do k1=1,n3b1m
               do k2=1,n3b1m
                  do k12=1,n3b2m
                     ii_f=0
                     if(itype1.ne.itype2) then
                        ii_f=1
                     endif
                     if(itype1.eq.itype2.and.k1.le.k2) then
                        ii_f=1
                     endif
                     if(ii_f.gt.0) then
                        num=num+1
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! nfeat0m_2=num

      ! cccccccccccccccccccccccccccccccc
      ! reading grid part
      ! allocate(li_ff_grid2(0:n2bm+1,ff(ff_idx)%li_ff_num_type))
      allocate(ff(ff_idx)%li_ff_grid2_2(2,n2bm+1,ff(ff_idx)%li_ff_num_type))
      allocate(ff(ff_idx)%li_ff_grid31_2(2,n3b1m,ff(ff_idx)%li_ff_num_type))
      allocate(ff(ff_idx)%li_ff_grid32_2(2,n3b2m,ff(ff_idx)%li_ff_num_type))

      ! for 2b
      do kkk=1,ff(ff_idx)%li_ff_num_type  ! center atom???
         iflag_grid=ff(ff_idx)%li_ff_iflag_grid_type(kkk)

         if(iflag_grid.eq.3) then
            ! read grid2b_type3
            read(10,*) ff(ff_idx)%li_ff_iat_feat,ff(ff_idx)%temp1,ff(ff_idx)%temp2,ff(ff_idx)%temp3
            read(10,*) ff(ff_idx)%n2b_tmp

            if(ff(ff_idx)%n2b_tmp.ne.ff(ff_idx)%li_ff_n2b_type(kkk)) then
               write(6,*) "n2b_tmp.ne.li_ff_n2b_type,in grid2b_type3", ff(ff_idx)%n2b_tmp,ff(ff_idx)%li_ff_n2b_type(kkk)
               stop
            endif

            do i=1,ff(ff_idx)%li_ff_n2b_type(kkk)
               read(10,*) ff(ff_idx)%n2b_tmp_idx, &
                  ff(ff_idx)%li_ff_grid2_2(1,i,kkk), &
                  ff(ff_idx)%li_ff_grid2_2(2,i,kkk)
               if(ff(ff_idx)%li_ff_grid2_2(2,i,kkk).gt.ff(ff_idx)%li_ff_Rc_type(kkk)) then
                  write(6,*) "grid2_2.gt.Rc",ff(ff_idx)%li_ff_grid2_2(2,i,kkk),ff(ff_idx)%li_ff_Rc_type(kkk)
               endif
            enddo
         endif
      enddo

      ! for 3b
      do kkk=1,ff(ff_idx)%li_ff_num_type
         iflag_grid=ff(ff_idx)%li_ff_iflag_grid_type(kkk)

         if(iflag_grid.eq.3) then
            ! read grid3b_b1b2_type3
            read(10,*) ff(ff_idx)%li_ff_iat_feat,ff(ff_idx)%temp1,ff(ff_idx)%temp2,ff(ff_idx)%temp3
            read(10,*) ff(ff_idx)%n3b2_tmp

            if(ff(ff_idx)%n3b2_tmp.ne.ff(ff_idx)%li_ff_n3b2_type(kkk)) then
               write(6,*) "n3b2_tmp.ne.li_ff_n3b2_type,in grid2b_type3", ff(ff_idx)%n3b2_tmp,ff(ff_idx)%li_ff_n3b2_type(kkk)
               stop
            endif

            do i=1,ff(ff_idx)%li_ff_n3b2_type(kkk)
               read(10,*) ff(ff_idx)%n3b2_tmp_idx, &
                  ff(ff_idx)%li_ff_grid32_2(1,i,kkk), &
                  ff(ff_idx)%li_ff_grid32_2(2,i,kkk)
               if(ff(ff_idx)%li_ff_grid32_2(2,i,kkk).gt.ff(ff_idx)%li_ff_Rc_type(kkk)) then
                  write(6,*) "grid32_2.gt.Rc",ff(ff_idx)%li_ff_grid32_2(2,i,kkk),ff(ff_idx)%li_ff_Rc_type(kkk)
               endif
            enddo
         endif
      enddo

      do kkk=1,ff(ff_idx)%li_ff_num_type
         iflag_grid=ff(ff_idx)%li_ff_iflag_grid_type(kkk)

         if(iflag_grid.eq.3) then
            ! read grid3b_cb12_type3
            read(10,*) ff(ff_idx)%li_ff_iat_feat,ff(ff_idx)%temp1,ff(ff_idx)%temp2,ff(ff_idx)%temp3
            read(10,*) ff(ff_idx)%n3b1_tmp

            if(ff(ff_idx)%n3b1_tmp.ne.ff(ff_idx)%li_ff_n3b1_type(kkk)) then
               write(6,*) "n3b1_tmp.ne.li_ff_n3b1_type,in grid2b_type3", ff(ff_idx)%n3b1_tmp,ff(ff_idx)%li_ff_n3b1_type(kkk)
               stop
            endif

            do i=1,ff(ff_idx)%li_ff_n3b1_type(kkk)
               read(10,*) ff(ff_idx)%n3b1_tmp_idx, &
                  ff(ff_idx)%li_ff_grid31_2(1,i,kkk), &
                  ff(ff_idx)%li_ff_grid31_2(2,i,kkk)
               if(ff(ff_idx)%li_ff_grid31_2(2,i,kkk).gt.ff(ff_idx)%li_ff_Rc_type(kkk)) then
                  write(6,*) "grid31_2.gt.Rc",ff(ff_idx)%li_ff_grid31_2(2,i,kkk),ff(ff_idx)%li_ff_Rc_type(kkk)
               endif
            enddo
         endif
      enddo

      ! reading shift and scale
      allocate(ff(ff_idx)%li_ff_feat2_shift(tmp_nfeat2,ff(ff_idx)%li_ff_num_type))
      allocate(ff(ff_idx)%li_ff_feat2_scale(tmp_nfeat2,ff(ff_idx)%li_ff_num_type))
      do itype=1,ff(ff_idx)%li_ff_num_type
         do j=1,tmp_nfeat2
            read(10,*) ff(ff_idx)%li_ff_feat2_shift(j,itype),ff(ff_idx)%li_ff_feat2_scale(j,itype)
         enddo
      enddo

      ! reading weight part         ! be careful for li_ff_itype_atom_sumfe
      allocate(ff(ff_idx)%w_feat(tmp_nfeat2,ff(ff_idx)%li_ff_num_type))  ! Temporarily use
      do itype=1,ff(ff_idx)%li_ff_num_type
         do j=1,tmp_nfeat2
            read(10,*) jj,ff(ff_idx)%w_feat(j,ff(ff_idx)%li_ff_num_type)
            ff(ff_idx)%w_feat(j,ff(ff_idx)%li_ff_num_type)=ff(ff_idx)%w_feat(j,ff(ff_idx)%li_ff_num_type)**2
         enddo
      enddo

      ! reading linear_fib.ntype
      read(10,*) ff(ff_idx)%nfeat2tot
      allocate(ff(ff_idx)%BB(ff(ff_idx)%nfeat2tot))
      do i=1,ff(ff_idx)%nfeat2tot
         read(10,*) ff(ff_idx)%ifeat2tot,ff(ff_idx)%BB(i)
      enddo

      ! reading vdw_fitB.ntype
      read(10,*) ntype_vdw,ff(ff_idx)%nterm
      if(ff(ff_idx)%nterm.gt.2) then
         write(6,*) "nterm.gt.2,stop"
         stop
      endif
      if(ntype_vdw.ne.ff(ff_idx)%li_ff_num_type) then
         write(6,*) "ntypes not same in vdw_fitB.ntype,something wrong"
         stop
      endif
      allocate(ff(ff_idx)%rad_atom(ff(ff_idx)%li_ff_num_type))
      allocate(ff(ff_idx)%E_ave_vdw(ff(ff_idx)%li_ff_num_type))
      allocate(ff(ff_idx)%wp_atom(ff(ff_idx)%li_ff_num_type,ff(ff_idx)%li_ff_num_type,2))
      do itype=1,ff(ff_idx)%li_ff_num_type
         read(10,*) itype_vdw,ff(ff_idx)%rad_atom(itype),ff(ff_idx)%E_ave_vdw(itype),&
            ((ff(ff_idx)%wp_atom(i,itype,j),i=1,ff(ff_idx)%li_ff_num_type),j=1,ff(ff_idx)%nterm)
      enddo

   end subroutine

   subroutine li_ff_deallocate(ff_idx) bind(c,name="li_ff_deallocate")
      ! deallocate the ff pointers
      integer, intent(in) :: ff_idx

      deallocate(ff(ff_idx)%li_ff_grid2_2)
      deallocate(ff(ff_idx)%li_ff_grid31_2)
      deallocate(ff(ff_idx)%li_ff_grid32_2)
      deallocate(ff(ff_idx)%w_feat)
      deallocate(ff(ff_idx)%BB)
      deallocate(ff(ff_idx)%li_ff_feat2_shift)
      deallocate(ff(ff_idx)%li_ff_feat2_scale)
      deallocate(ff(ff_idx)%li_ff_itype_atom_sumfe)
      deallocate(ff(ff_idx)%rad_atom)
      deallocate(ff(ff_idx)%E_ave_vdw)
      deallocate(ff(ff_idx)%wp_atom)

   end subroutine

end module li12_ff_mod