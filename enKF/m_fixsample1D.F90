module m_fixsample1D
contains
subroutine fixsample1D(E,n,m)
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   real   , intent(inout) :: E(n,m)

   integer iens,i
   real, allocatable :: average(:), variance(:)
   real var

   allocate(average(n), variance(n))

   average=0.0
   do iens=1,m
      average(:)=average(:)+E(:,iens)
   enddo
   average=average/float(m)

   do iens=1,m
      E(:,iens)=E(:,iens)-average(:)
   enddo

   variance=0.0
   do iens=1,m
      variance(:) = variance(:) + E(:,iens)**2.
   enddo

   print *,'variance '
   var=sum(variance)/real(m*n)
   print *,'1D var=',var
   print *

   do i=1,n
      variance(i)=1.0/sqrt( variance(i)/float(m) )
   enddo

   do iens=1,m
      do i=1,n
         E(i,iens)=variance(i)*E(i,iens)
      enddo
   enddo

   deallocate(average,variance)

end subroutine
end module
