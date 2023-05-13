      subroutine find_feature_deepMD2(natom,itype_atom,Rc,n2b,weight_rterm, &
        num_neigh,list_neigh, &
        dR_neigh,iat_neigh,ntype,grid2,wgauss, &
        feat_all,dfeat_all,nfeat0m,m_neigh,nfeat_atom)
        use mod_mpi
        implicit none
      integer ntype
      integer natom,n2b(ntype)
      real*8 grid2(200,50),wgauss(200,50)
      integer m_neigh
      integer itype_atom(natom)
      real*8 Rc(ntype),Rc2(ntype),Rm(ntype)
      integer M_type(ntype)
      real*8 weight_rterm(ntype)
      real*8 dR_neigh(3,m_neigh,ntype,natom)
      real*8 dR_neigh_alltype(3,m_neigh,natom)
      integer iat_neigh(m_neigh,ntype,natom),list_neigh(m_neigh,ntype,natom)
      integer num_neigh(ntype,natom)
      integer num_neigh_alltype(natom)
      integer nperiod(3)
      integer iflag,i,j,num,iat,itype
      integer i1,i2,i3,itype1,itype2,j1,j2,iat1,iat2
      real*8 d,dx1,dx2,dx3,dd
      real*8 pi,pi2,x,f1
      real*8 Rt,ff,dt
      integer iflag_grid
      integer itype0,nfeat0m

      integer ind_f(2,m_neigh,ntype,natom)
      real*8 f32(2),df32(2,2,3)
      integer inf_f32(2),k,k1,k2,k12,j12,ii_f,jj,jj1,jj2,nneigh,ii
      real*8 y2
      integer itype12,ind_f32(2)
      integer ind_all_neigh(m_neigh,ntype,natom),list_neigh_alltype(m_neigh,natom)

      real*8 feat_all(nfeat0m,natom_n)
      real*8 dfeat_all(nfeat0m,natom_n,m_neigh,3)
      integer nfeat_atom(natom)
      integer nfeat_atom_tmp(natom)
      integer ierr

      integer M1,ii1,ii2
      real*8 dx(3),df2,f2,y,s,sum,dff
      real*8 poly(100),dpoly(100)
      real*8 ww(0:3)

      real*8, allocatable, dimension (:,:) :: tensor
      real*8, allocatable, dimension (:,:,:,:) :: dtensor
      real*8, allocatable, dimension (:,:) :: dsum


      num_neigh_alltype=0
      do iat=1,natom
      if(mod(iat-1,nnodes).eq.inode-1) then
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
      endif
      enddo

!ccccccccccccccccccccccccccccccccccccccccc

      pi=4*datan(1.d0)
      pi2=2*pi

      iat1=0
      do 3000 iat=1,natom
      if(mod(iat-1,nnodes).eq.inode-1) then
      iat1=iat1+1   
       
      itype0=itype_atom(iat)
       M1=n2b(itype0)

       nneigh=num_neigh_alltype(iat)

      allocate(tensor(0:3,M1*ntype))
      allocate(dtensor(0:3,M1*ntype,nneigh,3))
      allocate(dsum(nneigh,3))

      tensor=0.d0
      dtensor=0.d0
      

      do 1000 itype=1,ntype
      do 1000 j=1,num_neigh(itype,iat)

      jj=ind_all_neigh(j,itype,iat)

      dx(1)=dR_neigh(1,j,itype,iat)
      dx(2)=dR_neigh(2,j,itype,iat)
      dx(3)=dR_neigh(3,j,itype,iat)
      dd=dx(1)**2+dx(2)**2+dx(3)**2
      d=dsqrt(dd)

      if(d.gt.Rc(itype0)) goto 1001


      do i=1,M1
      Rt=wgauss(i,itype0)
      dt=grid2(i,itype)
      f1=exp(-((d-dt)/Rt)**2)
      x=pi*d/Rc(itype0)
      f2=0.5*(cos(x)+1)
      ff=f1*f2/d
      dff=-2*f1*f2/d*(d-dt)/Rt**2
      dff=dff-0.5*sin(x)*pi/Rc(itype0)*f1/d
      dff=dff-f1*f2/d**2
      poly(i)=ff
      dpoly(i)=dff    ! need to add d_d/d_dx=dx(:)/d
      enddo


      do i=1,M1
      ii=i+(itype-1)*M1
      tensor(0,ii)=tensor(0,ii)+d*poly(i)
      tensor(1,ii)=tensor(1,ii)+dx(1)*poly(i)
      tensor(2,ii)=tensor(2,ii)+dx(2)*poly(i)
      tensor(3,ii)=tensor(3,ii)+dx(3)*poly(i)

      dtensor(0,ii,jj,:)=dtensor(0,ii,jj,:)+(poly(i)+d*dpoly(i))*dx(:)/d
      dtensor(1,ii,jj,:)=dtensor(1,ii,jj,:)+dx(1)*dpoly(i)*dx(:)/d
      dtensor(1,ii,jj,1)=dtensor(1,ii,jj,1)+poly(i)
      dtensor(2,ii,jj,:)=dtensor(2,ii,jj,:)+dx(2)*dpoly(i)*dx(:)/d
      dtensor(2,ii,jj,2)=dtensor(2,ii,jj,2)+poly(i)
      dtensor(3,ii,jj,:)=dtensor(3,ii,jj,:)+dx(3)*dpoly(i)*dx(:)/d
      dtensor(3,ii,jj,3)=dtensor(3,ii,jj,3)+poly(i)

      dtensor(0,ii,1,:)=dtensor(0,ii,1,:)-(poly(i)+d*dpoly(i))*dx(:)/d
      dtensor(1,ii,1,:)=dtensor(1,ii,1,:)-dx(1)*dpoly(i)*dx(:)/d
      dtensor(1,ii,1,1)=dtensor(1,ii,1,1)-poly(i)
      dtensor(2,ii,1,:)=dtensor(2,ii,1,:)-dx(2)*dpoly(i)*dx(:)/d
      dtensor(2,ii,1,2)=dtensor(2,ii,1,2)-poly(i)
      dtensor(3,ii,1,:)=dtensor(3,ii,1,:)-dx(3)*dpoly(i)*dx(:)/d
      dtensor(3,ii,1,3)=dtensor(3,ii,1,3)-poly(i)
      enddo

1001  continue
1000  continue

      ww=1.d0
      ww(0)=weight_rterm(itype0)

      ii=0
      do ii1=1,M1*ntype
      do ii2=1,ii1
      ii=ii+1
      sum=0.d0
      dsum=0.d0
      do k=0,3
      sum=sum+tensor(k,ii1)*tensor(k,ii2)*ww(k)
      dsum(:,1)=dsum(:,1)+(dtensor(k,ii1,:,1)*tensor(k,ii2)+tensor(k,ii1)*dtensor(k,ii2,:,1))*ww(k)
      dsum(:,2)=dsum(:,2)+(dtensor(k,ii1,:,2)*tensor(k,ii2)+tensor(k,ii1)*dtensor(k,ii2,:,2))*ww(k)
      dsum(:,3)=dsum(:,3)+(dtensor(k,ii1,:,3)*tensor(k,ii2)+tensor(k,ii1)*dtensor(k,ii2,:,3))*ww(k)
      enddo

      feat_all(ii,iat1)=sum
      dfeat_all(ii,iat1,1:nneigh,1)=dsum(1:nneigh,1)
      dfeat_all(ii,iat1,1:nneigh,2)=dsum(1:nneigh,2)
      dfeat_all(ii,iat1,1:nneigh,3)=dsum(1:nneigh,3)
      enddo
      enddo

      deallocate(tensor)
      deallocate(dtensor)
      deallocate(dsum)
      nfeat_atom_tmp(iat)=M1*ntype*(M1*ntype+1)/2

      endif
3000  continue

      call mpi_allreduce(nfeat_atom_tmp,nfeat_atom,natom,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

!ccccccccccccccccccccccccccccccccccccc

      return
      end subroutine find_feature_deepMD2



