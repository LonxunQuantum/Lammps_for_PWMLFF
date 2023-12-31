subroutine find_feature_2b_type2(natom,itype_atom,Rc,n2b, &
        num_neigh,list_neigh, &
        dR_neigh,iat_neigh,ntype,grid2, &
!        feat_all,dfeat_allR,nfeat0m,m_neigh,n2bm,nfeat_atom)
        feat_all,dfeat_all,nfeat0m,m_neigh,n2bm,nfeat_atom)
    ! ******************************************
    !      feature calc on a SINGLE core  
    ! ******************************************
    use mod_mpi
    implicit none
    integer ntype
    integer natom,n2b(ntype)
    integer m_neigh
    integer itype_atom(natom)
    
    real*8 Rc(ntype)
    real*8 dR_neigh(3,m_neigh,ntype,natom)
    real*8 dR_neigh_alltype(3,m_neigh,natom)
    integer iat_neigh(m_neigh,ntype,natom),list_neigh(m_neigh,ntype,natom)
    integer num_neigh(ntype,natom)
    integer num_neigh_alltype(natom)
    integer nperiod(3)
    integer iflag,i,j,num,iat,itype
    integer i1,i2,i3,itype1,itype2,j1,j2,iat1,iat2
    real*8 d,dx1,dx2,dx3,dx,dy,dz,dd
    real*8 grid2(0:n2bm+1,ntype)
    real*8 pi,pi2,x,f1
    integer iflag_grid
    integer itype0,nfeat0m,n2bm

    integer ind_f(2,m_neigh,ntype,natom)
    real*8 f32(2),df32(2,2,3)
    integer inf_f32(2),k,k1,k2,k12,j12,ii_f,jj,jj1,jj2,nneigh,ii
    real*8 y,y2
    integer itype12,ind_f32(2)
    integer ind_all_neigh(m_neigh,ntype,natom),list_neigh_alltype(m_neigh,natom)

    real*8 feat2(n2bm,ntype,natom)
    real*8 dfeat2(n2bm,ntype,natom,m_neigh,3)
    real*8 feat_all(nfeat0m,natom),dfeat_all(nfeat0m,natom,m_neigh,3)


    integer nfeat_atom(natom)
    integer nfeat_atom_tmp(natom)
    integer ierr

    
    num_neigh_alltype=0
    
    do iat=1,natom
        !if(mod(iat-1,nnodes).eq.inode-1) then
            num=1
            list_neigh_alltype(1,iat)=iat   ! the first neighbore is itself
            dR_neigh_alltype(:,1,iat)=0.d0

            do  itype=1,ntype
            do   j=1,num_neigh(itype,iat)
            num=num+1
                if(num.gt.m_neigh) then
                write(6,*) "total num_neigh.gt.m_neigh,stop",m_neigh
                stop
                endif
            ind_all_neigh(j,itype,iat)=num
            list_neigh_alltype(num,iat)=list_neigh(j,itype,iat)
            dR_neigh_alltype(:,num,iat)=dR_neigh(:,j,itype,iat)
            enddo
            enddo
            num_neigh_alltype(iat)=num
        !endif
    enddo

    !ccccccccccccccccccccccccccccccccccccccccc

    pi=4*datan(1.d0)
    pi2=2*pi

    feat2=0.d0
    dfeat2=0.d0

    iat1=0
    
    do iat=1,natom
        !if(mod(iat-1,nnodes).eq.inode-1) then
            iat1=iat1+1

            itype0=itype_atom(iat)

            do 1000 itype=1,ntype
            do 1000 j=1,num_neigh(itype,iat)

            jj=ind_all_neigh(j,itype,iat)

            dd=dR_neigh(1,j,itype,iat)**2+dR_neigh(2,j,itype,iat)**2+dR_neigh(3,j,itype,iat)**2
            d=dsqrt(dd)


            do k=1,n2b(itype0)+1
                if(grid2(k,itype0).ge.d) exit
            enddo
            k=k-1


            if(k.gt.n2b(itype0))  k=n2b(itype0)
                !  This point will be in the intervalL [grid2(k),grid2(k+1)]
                !  It will have feature contribution:
                !  feat2(k,itype,iat),feat2(k+1,itype,iat)
            if(k.le.n2b(itype0)-2) then   ! k+2 feature
                x=(d-grid2(k,itype0))/(grid2(k+3,itype0)-grid2(k,itype0))
                y=(x-0.5d0)*pi2
                f1=0.5d0*(cos(y)+1)
                feat2(k+2,itype,iat1)=feat2(k+2,itype,iat1)+f1
                y2=-pi*sin(y)/(d*(grid2(k+3,itype0)-grid2(k,itype0)))
                dfeat2(k+2,itype,iat1,jj,:)=dfeat2(k+2,itype,iat1,jj,:)+y2*dR_neigh(:,j,itype,iat)
                dfeat2(k+2,itype,iat1,1,:)=dfeat2(k+2,itype,iat1,1,:)-y2*dR_neigh(:,j,itype,iat)
                    ! Note, (k+1,itype) is the feature inde
            endif

            if(k.le.n2b(itype0)-1) then   ! k+1 feature
                x=(d-grid2(k-1,itype0))/(grid2(k+2,itype0)-grid2(k-1,itype0))
                y=(x-0.5d0)*pi2
                f1=0.5d0*(cos(y)+1)
                feat2(k+1,itype,iat1)=feat2(k+1,itype,iat1)+f1
                y2=-pi*sin(y)/(d*(grid2(k+2,itype0)-grid2(k-1,itype0)))
                dfeat2(k+1,itype,iat1,jj,:)=dfeat2(k+1,itype,iat1,jj,:)+y2*dR_neigh(:,j,itype,iat)
                dfeat2(k+1,itype,iat1,1,:)=dfeat2(k+1,itype,iat1,1,:)-y2*dR_neigh(:,j,itype,iat)
                    ! Note, (k+1,itype) is the feature inde
            endif
            
            if(k.gt.0.and.k.le.n2b(itype0)) then   ! k feature
                x=(d-grid2(k-2,itype0))/(grid2(k+1,itype0)-grid2(k-2,itype0))
                y=(x-0.5d0)*pi2
                f1=0.5d0*(cos(y)+1)
                feat2(k,itype,iat1)=feat2(k,itype,iat1)+f1
                y2=-pi*sin(y)/(d*(grid2(k+1,itype0)-grid2(k-2,itype0)))
                dfeat2(k,itype,iat1,jj,:)=dfeat2(k,itype,iat1,jj,:)+y2*dR_neigh(:,j,itype,iat)
                dfeat2(k,itype,iat1,1,:)=dfeat2(k,itype,iat1,1,:)-y2*dR_neigh(:,j,itype,iat)
                    ! Note, (k+1,itype) is the feature inde
            endif


            !cccccccccccc So, one Rij will always have two features k, k+1  (1,2)
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
            1000  continue

        !endif
    enddo

    !cccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccc
    !   Now, we collect everything together, collapse the index (k,itype)
    !   feat2, into a single feature. 

    nfeat_atom_tmp=0

    iat1=0
    
    do iat=1,natom
        !if(mod(iat-1,nnodes).eq.inode-1) then
            iat1=iat1+1

            itype0=itype_atom(iat)
            nneigh=num_neigh_alltype(iat)

            num=0
            do itype=1,ntype
                do k=1,n2b(itype0)
                    num=num+1
                    feat_all(num,iat1)=feat2(k,itype,iat1)
                    dfeat_all(num,iat1,1:nneigh,:)=dfeat2(k,itype,iat1,1:nneigh,:)
                enddo
            enddo

            nfeat_atom(iat)=num
            
            if(num.gt.nfeat0m) then
                write(6,*) "num.gt.nfeat0m,stop",num,nfeat0m
                stop
            endif

        !endif
    enddo 


!call mpi_allreduce(nfeat_atom_tmp,nfeat_atom,natom,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

!ccccccccccccccccccccccccccccccccccc
!  Now, we have to redefine the dfeat_all in another way. 
!  dfeat_all(i,iat,jneigh,3) means:
!  d_ith_feat_of_iat/d_R(jth_neigh_of_iat)
!  dfeat_allR(i,iat,jneigh,3) means:
!  d_ith_feat_of_jth_neigh/d_R(iat)
!  Now, we just output dfeat_allR
!cccccccccccccccccccccccccccccccccccccc

      return
      end subroutine find_feature_2b_type2



