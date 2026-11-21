subroutine Output
 use arrays 
!
! this subroutine to prepare the Input File
!
   nint=0
   write(ipr,'(3i7)')nnod,nelm,nint
  do nnum=1,nnod
   indicator=' '
  do i=1,nn1
  do j=1,nn2
   if(indicator==' '.and.node(i,j)%nnum==nnum)then
   write(ipr,'(3f12.5)')node(i,j)%x
   write(13,'(i7,1x,3f12.5)')nnum,node(i,j)%x
   indicator='*'
   endif
  enddo ; enddo ; enddo
   do i=1,nelm
   select case(nod)
   case(3)
   write(ipr,'(4i7,1x,f10.5)')(element(i)%node(j),j=1,nod),element(i)%kode,&
   element(i)%value
   case(4)
   write(ipr,'(5i7,1x,f10.5)')(element(i)%node(j),j=1,nod),element(i)%kode,&
   element(i)%value
   end select
   enddo
   write(14,'(i7,1x,3f10.5)')(i,element(i)%x,i=1,nelm)

end subroutine Output   
