module calc_ftype1
   ! ******************************************
   !      feature calc on a SINGLE core
   ! ******************************************
   use li_ff_mod, only: ff

   use mod_data, only : natoms, ntypes, catype

   IMPLICIT NONE
   integer :: i,itype
   integer :: max_neigh,max_neigh_M
   real(8) :: E_tolerance
   integer :: recalc_grid
   real(8),allocatable,dimension (:,:,:) :: grid2_2
   integer :: n2b_t,it
   integer :: n2b_type(100),n2bm                              ! numOf2bfeat
   real(8) :: Rc_M
   
   integer :: m_neigh,m_neigh1
   integer :: num,num_M,natom1
   integer :: n2b,nfeat0m
   integer :: nfeat0(100)       ! ?
   integer :: iat_type(100)
   integer :: iflag_grid,iflag_ftype,iflag_grid_type(100)
   real(8) :: Rc_type(100),Rm_type(100),fact_grid_type(100),dR_grid1_type(100)

   !ccccccccccccccccccccc  The variables to be used in global feature type
   real(8), allocatable, dimension(:,:) :: feat_M1
   real(8), allocatable, dimension(:,:,:,:) :: dfeat_M1
   integer, allocatable, dimension(:,:) :: list_neigh_alltypeM1
   integer, allocatable, dimension(:) :: num_neigh_alltypeM1
   integer :: nfeat0M1
   !ccccccccccccccccccccc  The variables to be used in global feature type

contains

   subroutine load_model_type1(ff_idx)
      integer, intent(in) :: ff_idx
      integer :: kkk

      ! gen_2b_feature.in
      Rc_M=ff(ff_idx)%li_ff_Rc_M
      ! m_neigh=ff(ff_idx)%li_ff_max_neigh
      ! ntype=ff(ff_idx)%li_ff_num_type
      do i=1,ff(ff_idx)%li_ff_num_type
         iat_type(i)=ff(ff_idx)%li_ff_iat_type(i)
         Rc_type(i)=ff(ff_idx)%li_ff_Rc_type(i)
         Rm_type(i)=ff(ff_idx)%li_ff_Rm_type(i)
         iflag_grid_type(i)=ff(ff_idx)%li_ff_iflag_grid_type(i)
         fact_grid_type(i)=ff(ff_idx)%li_ff_fact_grid_type(i)
         dR_grid1_type(i)=ff(ff_idx)%li_ff_dR_grid1_type(i)
         n2b_type(i)=ff(ff_idx)%li_ff_n2b_type(i)
      enddo

      E_tolerance=ff(ff_idx)%li_ff_E_tolerance
      iflag_ftype=ff(ff_idx)%li_ff_iflag_ftype
      recalc_grid=ff(ff_idx)%li_ff_recalc_grid

      ! m_neigh1=m_neigh
      ! ccccccccccccccccccccccccccccccccccccccccccccccccc
      ! calculate features of all types
      n2bm=0
      do i=1,ntypes
         if(n2b_type(i).gt.n2bm) then
            n2bm=n2b_type(i)
         endif
      enddo
      !cccccccccccccccccccccccccccccccccccccccccccccccc
      nfeat0m=ntypes*n2bm
      do itype=1,ntypes
         nfeat0(itype)=n2b_type(itype)*ntypes
      enddo
      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.allocated(grid2_2)) then
         allocate(grid2_2(2,n2bm+1,ntypes))
      endif
      ! grid2_2=ff(ff_idx)%li_ff_grid2_2

      do kkk=1,ntypes    ! center atom
         iflag_grid=iflag_grid_type(kkk)
         n2b=n2b_type(kkk)

         if(iflag_grid.eq.3) then
            ! read grid2b_type3
            n2b_t=ff(ff_idx)%n2b_tmp
            do i=1,n2b
               it=ff(ff_idx)%n2b_tmp_idx
               grid2_2(1,i,kkk)=ff(ff_idx)%li_ff_grid2_2(1,i,kkk)
               grid2_2(2,i,kkk)=ff(ff_idx)%li_ff_grid2_2(2,i,kkk)
            enddo
         endif
      enddo ! kkk=1,ntypes

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  FInish the initial grid treatment
   end subroutine load_model_type1

   subroutine set_image_info_type1(ff_idx)
      integer, intent(in) :: ff_idx

      m_neigh=ff(ff_idx)%li_ff_max_neigh
      m_neigh1=m_neigh
      natom1=natoms
      
   end subroutine set_image_info_type1

   subroutine gen_feature_type1(num_neigh,list_neigh,dR_neigh)
      integer, dimension(ntypes,natoms), intent(in) :: num_neigh
      integer, dimension(m_neigh,ntypes,natoms), intent(in) :: list_neigh
      real(8), dimension(3,m_neigh,ntypes,natoms), intent(in) :: dR_neigh

      integer(4) :: j
      integer :: ii,jj,jjm,iat,iat1

      integer,allocatable,dimension (:,:) :: map2neigh_alltypeM      ! -->list_tmp
      integer,allocatable,dimension (:,:) :: list_tmp                ! -->num_M
      integer,allocatable,dimension (:,:,:) :: map2neigh_M
      integer,allocatable,dimension (:,:,:) :: list_neigh_M
      integer,allocatable,dimension (:,:) :: num_neigh_M
      integer,allocatable,dimension (:,:) :: list_neigh_alltype
      integer,allocatable,dimension (:) :: num_neigh_alltype
      real*8,allocatable,dimension (:,:) :: feat
      real*8,allocatable,dimension (:,:,:,:) :: dfeat
      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! the dimension of these array, should be changed to natom_n
      ! instead of natoms. Each process only needs to know its own natom_n

      !write(*,*) "m_neigh,ntypes,natoms", m_neigh,ntypes,natoms
      if (allocated(feat_M1)) then
         deallocate(feat_M1)
      endif
      if (allocated(dfeat_M1)) then
         deallocate(dfeat_M1)
      endif
      if (allocated(list_neigh_alltypeM1)) then
         deallocate(list_neigh_alltypeM1)
      endif
      if (allocated(num_neigh_alltypeM1)) then
         deallocate(num_neigh_alltypeM1)
      endif

      allocate(feat_M1(nfeat0m,natoms))
      allocate(dfeat_M1(nfeat0m,natoms,m_neigh,3))
      allocate(list_neigh_alltypeM1(m_neigh,natoms))
      allocate(num_neigh_alltypeM1(natoms))
      ! allocate(list_neigh(m_neigh,ntypes,natoms))
      ! allocate(dR_neigh(3,m_neigh,ntypes,natoms))   ! d(neighbore)-d(center) in xyz
      ! allocate(num_neigh(ntypes,natoms))
      allocate(map2neigh_M(m_neigh,ntypes,natoms)) ! from list_neigh(of this feature) to list_neigh_all (of Rc_M
      allocate(list_neigh_M(m_neigh,ntypes,natoms)) ! the neigh list of Rc_M
      allocate(num_neigh_M(ntypes,natoms))
      allocate(list_neigh_alltype(m_neigh,natoms))
      allocate(num_neigh_alltype(natoms))
      allocate(map2neigh_alltypeM(m_neigh,natoms)) ! from list_neigh(of this feature) to list_neigh_all (of Rc_M
      allocate(list_tmp(m_neigh,ntypes))
      allocate(feat(nfeat0m,natoms))         ! each note, only know its own feat
      allocate(dfeat(nfeat0m,natoms,m_neigh,3))  ! dfeat is the derivative from the neighboring dR,

      feat = 0.d0
      dfeat = 0.d0
      feat_M1 = 0.d0
      dfeat_M1 = 0.d0

      max_neigh=-1
      num_neigh_alltype=0
      max_neigh_M=-1
      num_neigh_alltypeM1=0
      list_neigh_alltypeM1=0

      do iat=1,natoms
         list_neigh_alltype(1,iat)=iat
         list_neigh_alltypeM1(1,iat)=iat

         num_M=1
         do itype=1,ntypes
            ! do j=1,num_neigh_M(itype,iat)     ! need to check
            do j=1,num_neigh(itype,iat)
               num_M=num_M+1

               if(num_M.gt.m_neigh) then
                  write(6,*) "Error: total num_neigh > m_neigh. Stop",m_neigh
                  stop
               endif
               ! list_neigh_alltypeM1(num_M,iat)=list_neigh_M(j,itype,iat)
               list_neigh_alltypeM1(num_M,iat)=list_neigh(j,itype,iat)
               list_tmp(j,itype)=num_M
            enddo
         enddo

         num=1
         map2neigh_alltypeM(1,iat)=1
         do itype=1,ntypes
            do j=1,num_neigh(itype,iat)
               num=num+1
               list_neigh_alltype(num,iat)=list_neigh(j,itype,iat)
               ! map2neigh_alltypeM(num,iat)=list_tmp(map2neigh_M(j,itype,iat),itype)
               map2neigh_alltypeM(num,iat)=num
               ! the neighbore atom, from list_neigh_alltype list to list_neigh_alltypeM list
               ! map2neigh_M(j,itype,iat), maps the jth neigh in list_neigh(Rc) to jth' neigh in list_neigh_M(Rc_M)
            enddo
         enddo

         num_neigh_alltype(iat)=num
         num_neigh_alltypeM1(iat)=num_M

         if(num.gt.max_neigh) max_neigh=num
         if(num_M.gt.max_neigh_M) max_neigh_M=num_M
      enddo  ! iat

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! This num_neigh_alltype(iat) include itself !

      !   if (inode .eq. 1) then
      !       write(*,*) "before type3, dR_neigh:"
      !       do i = 1,4
      !           write(*,*) dR_neigh(1:3,i,1,1)
      !       enddo
      !   endif

      if(iflag_ftype.eq.3) then
         !  iflag_ftype.eq.3, the sin peak span over the two ends specified by grid31_2,grid32_2
         !  So, there could be many overlaps between different sin peaks
         call find_feature_2b_type3(num_neigh,dR_neigh,m_neigh,list_neigh,&
            n2b_type,n2bm,grid2_2,&
            feat,dfeat,nfeat0m)
      endif
      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iat1=1,natoms
         do ii=1,nfeat0m
            feat_M1(ii,iat1)=feat(ii,iat1)
         enddo
      enddo

      iat1=0
      do iat=1,natoms
         iat1=iat1+1
         do jj=1,num_neigh_alltype(iat)
            jjm=map2neigh_alltypeM(jj,iat)
            do ii=1,nfeat0m
               dfeat_M1(ii,iat1,jjm,1)=dfeat(ii,iat1,jj,1)  ! this is the feature stored in neigh list of Rc_M
               dfeat_M1(ii,iat1,jjm,2)=dfeat(ii,iat1,jj,2)
               dfeat_M1(ii,iat1,jjm,3)=dfeat(ii,iat1,jj,3)
            enddo
         enddo
      enddo

      nfeat0M1=nfeat0m    ! the number of features for feature type 1

      ! deallocate(list_neigh)
      ! deallocate(dR_neigh)
      ! deallocate(num_neigh)
      deallocate(list_neigh_alltype)
      deallocate(num_neigh_alltype)

      deallocate(list_neigh_M)
      deallocate(num_neigh_M)
      deallocate(map2neigh_M)
      ! deallocate(feat_M1)  !
      ! deallocate(dfeat_M1) !
      ! deallocate(list_neigh_alltypeM1) !
      ! deallocate(num_neigh_alltypeM1)  !
      deallocate(map2neigh_alltypeM)
      deallocate(list_tmp)
      deallocate(feat)  
      deallocate(dfeat) 
      ! mem leak
      ! deallocate(grid2_2)  !
      !--------------------------------------------------------
   end subroutine gen_feature_type1

end module calc_ftype1
