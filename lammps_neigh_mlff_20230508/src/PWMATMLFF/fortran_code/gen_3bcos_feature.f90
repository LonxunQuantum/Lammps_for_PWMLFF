module calc_3bcos_feature
   use ff_mod, only: ff

   use mod_data, only : natoms, ntypes, catype

   IMPLICIT NONE
   integer :: i,j,itype
   integer :: max_neigh,max_neigh_M
   integer :: iat_type(100)
   integer :: n3b_type(50),n3bm
   real(8) :: Rc_M,Rc_type(50)
   real(8) :: eta_type(200,50),w_type(100,50),alamda_type(100,50)
   real(8) :: E_tolerance

   integer :: m_neigh,m_neigh4
   integer :: num,num_M,natom4
   integer :: nfeat0m,nfeat0(100)

   real(8),allocatable,dimension (:,:) :: feat_M4
   real(8),allocatable,dimension (:,:,:,:) :: dfeat_M4
   integer,allocatable,dimension (:,:) :: list_neigh_alltypeM4
   integer,allocatable,dimension (:) :: num_neigh_alltypeM4
   integer :: nfeat0M4


contains
   subroutine load_model_type4(ff_idx)
      integer, intent(in) :: ff_idx
      integer :: itype1,itype2

      ! gen_3bcos_feature.in
      Rc_M=ff(ff_idx)%ff_Rc_M
      do i=1,ff(ff_idx)%ff_num_type
         iat_type(i)=ff(ff_idx)%ff_iat_type(i)
         Rc_type(i)=ff(ff_idx)%ff_Rc_type(i)
         n3b_type(i)=ff(ff_idx)%ff_n3b_type(i)
         do j=1,n3b_type(i)
            eta_type(j,i)=ff(ff_idx)%ff_eta_type(j,i)
            w_type(j,i)=ff(ff_idx)%ff_w_type(j,i)
            alamda_type(j,i)=ff(ff_idx)%ff_alamda_type(j,i)
         enddo
      enddo
      E_tolerance=ff(ff_idx)%ff_E_tolerance

      !cccccccccccccccccccccccccccccccccccccccc
      n3bm=0
      do i=1,ntypes
         if(n3b_type(i).gt.n3bm) n3bm=n3b_type(i)
      enddo

      nfeat0m=0
      do itype=1,ntypes
         num=0
         do itype2=1,ntypes
            do itype1=1,itype2
               num=num+1
            enddo
         enddo
         nfeat0(itype)=num*n3b_type(itype)
         if(nfeat0(itype).gt.nfeat0m) nfeat0m=nfeat0(itype)
      enddo

      !cccccccccccccccccccccccccccccccccccccccc
      !  FInish the initial grid treatment
   end subroutine load_model_type4

   subroutine set_image_info_type4(ff_idx)
      integer, intent(in) :: ff_idx

      m_neigh=ff(ff_idx)%ff_max_neigh
      m_neigh4=m_neigh
      natom4=natoms

   end subroutine set_image_info_type4

   subroutine gen_3bcos_feature(num_neigh,list_neigh,dR_neigh)
      integer, dimension(ntypes,natoms), intent(in) :: num_neigh
      integer, dimension(m_neigh,ntypes,natoms), intent(in) :: list_neigh
      real(8), dimension(3,m_neigh,ntypes,natoms), intent(in) :: dR_neigh

      integer :: ii,jj,jjm,iat,iat1

      integer,allocatable,dimension (:,:) :: map2neigh_alltypeM
      integer,allocatable,dimension (:,:) :: list_tmp
      integer,allocatable,dimension (:,:,:) :: map2neigh_M
      integer,allocatable,dimension (:,:,:) :: list_neigh_M
      integer,allocatable,dimension (:,:) :: num_neigh_M
      integer,allocatable,dimension (:,:) :: list_neigh_alltype
      integer,allocatable,dimension (:) :: num_neigh_alltype
      real(8), allocatable,dimension (:,:) :: feat
      real(8), allocatable,dimension (:,:,:,:) :: dfeat
      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (allocated(dfeat_M4)) then
         deallocate(feat_M4)
         deallocate(dfeat_M4)
         deallocate(list_neigh_alltypeM4)
         deallocate(num_neigh_alltypeM4)
      endif

      allocate(map2neigh_M(m_neigh,ntypes,natoms)) ! from list_neigh(of this feature) to list_neigh_all (of Rc_M
      allocate(list_neigh_M(m_neigh,ntypes,natoms)) ! the neigh list of Rc_M
      allocate(num_neigh_M(ntypes,natoms))
      !   allocate(list_neigh(m_neigh,ntypes,natoms))
      !   allocate(num_neigh(ntypes,natoms))
      !   allocate(dR_neigh(3,m_neigh,ntypes,natoms))   ! d(neighbore)-d(center) in xyz
      allocate(list_neigh_alltype(m_neigh,natoms))
      allocate(num_neigh_alltype(natoms))

      allocate(list_neigh_alltypeM4(m_neigh,natoms))
      allocate(num_neigh_alltypeM4(natoms))
      allocate(map2neigh_alltypeM(m_neigh,natoms)) ! from list_neigh(of this feature) to list_neigh_all (of Rc_M
      allocate(list_tmp(m_neigh,ntypes))

      allocate(feat(nfeat0m,natoms))         ! each note, only know its own feat
      allocate(dfeat(nfeat0m,natoms,m_neigh,3))  ! dfeat is the derivative from the neighboring dR,
      allocate(feat_M4(nfeat0m,natoms))
      allocate(dfeat_M4(nfeat0m,natoms,m_neigh,3))

      feat = 0.d0
      dfeat = 0.d0
      feat_M4 = 0.d0
      dfeat_M4 = 0.d0

      max_neigh=-1
      num_neigh_alltype=0
      max_neigh_M=-1
      num_neigh_alltypeM4=0
      list_neigh_alltypeM4=0

      do iat=1,natoms
         list_neigh_alltype(1,iat)=iat
         list_neigh_alltypeM4(1,iat)=iat

         num_M=1
         do itype=1,ntypes
            ! do j=1,num_neigh_M(itype,iat)
            do j=1,num_neigh(itype,iat)
               num_M=num_M+1
               if(num_M.gt.m_neigh) then
                  write(6,*) "total num_neigh.gt.m_neigh,stop",m_neigh
                  stop
               endif
               ! list_neigh_alltypeM4(num_M,iat)=list_neigh_M(j,itype,iat)
               list_neigh_alltypeM4(num_M,iat)=list_neigh(j,itype,iat)
               list_tmp(j,itype)=num_M
            enddo
         enddo

         num=1
         map2neigh_alltypeM(1,iat)=1
         do itype=1,ntypes
            do j=1,num_neigh(itype,iat)
               num=num+1
               list_neigh_alltype(num,iat)=list_neigh(j,itype,iat)
               map2neigh_alltypeM(num,iat)=num
               ! map2neigh_alltypeM(num,iat)=list_tmp(map2neigh_M(j,itype,iat),itype)
               ! map2neigh_M(j,itype,iat), maps the jth neigh in list_neigh(Rc) to jth' neigh in list_neigh_M(Rc_M)
            enddo
         enddo

         num_neigh_alltype(iat)=num
         num_neigh_alltypeM4(iat)=num_M
         if(num.gt.max_neigh) max_neigh=num
         if(num_M.gt.max_neigh_M) max_neigh_M=num_M
      enddo  ! iat

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! This num_neigh_alltype(iat) include itself !
      !    dfeat=0.d0
      !    feat=0.d0

      call find_feature_3bcos(num_neigh,dR_neigh,m_neigh,list_neigh,&
         n3b_type,n3bm,eta_type,w_type,alamda_type,Rc_type,&
         feat,dfeat,nfeat0m)
      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iat1=1,natoms
         do ii=1,nfeat0m
            feat_M4(ii,iat1)=feat(ii,iat1)
         enddo
      enddo

      iat1=0
      do iat=1,natoms
         iat1=iat1+1
         do jj=1,num_neigh_alltype(iat)
            jjm=map2neigh_alltypeM(jj,iat)
            do ii=1,nfeat0m
               dfeat_M4(ii,iat1,jjm,1)=dfeat(ii,iat1,jj,1)  ! this is the feature stored in neigh list of Rc_M
               dfeat_M4(ii,iat1,jjm,2)=dfeat(ii,iat1,jj,2)
               dfeat_M4(ii,iat1,jjm,3)=dfeat(ii,iat1,jj,3)
            enddo
         enddo
      enddo

      nfeat0M4=nfeat0m    ! the number of features for feature type 1

    !   deallocate(list_neigh)
    !   deallocate(dR_neigh)
    !   deallocate(num_neigh)
      deallocate(list_neigh_alltype)
      deallocate(num_neigh_alltype)

      deallocate(list_neigh_M)
      deallocate(num_neigh_M)
      deallocate(map2neigh_M)
      ! deallocate(list_neigh_alltypeM4)
      ! deallocate(num_neigh_alltypeM4)
      deallocate(map2neigh_alltypeM)
      deallocate(list_tmp)
      deallocate(feat)
      deallocate(dfeat)

   end subroutine gen_3bcos_feature

end module calc_3bcos_feature
