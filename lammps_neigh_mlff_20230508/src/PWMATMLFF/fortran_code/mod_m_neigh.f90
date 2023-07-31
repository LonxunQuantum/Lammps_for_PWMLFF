module mod_m_neigh
   use mod_data, only : iflag_model
   use calc_ftype1, only : m_neigh1
   use calc_ftype2, only : m_neigh2
   use calc_deepMD_f, only : m_neigh_t5
   implicit none
   integer :: m_neigh                                  !模型所使用的最大近邻数

   contains
   
   subroutine load_m_neigh()
      
      if(iflag_model.eq.1) then
         ! linear model
         if(m_neigh1.eq.m_neigh2) then
            m_neigh=m_neigh1
         else
            write(*,*) 'm_neigh1 and m_neigh2 are not equal'
            stop
         endif
      endif
      
      if(iflag_model.eq.5) then
         ! dp model
         m_neigh=m_neigh_t5
      endif
          
   end subroutine load_m_neigh
   
end module mod_m_neigh
