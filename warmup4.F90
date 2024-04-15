  module mathOps
contains
  attributes(global) subroutine bigmat(n,x, y, bigz)
    implicit none
    integer :: n
    real :: x(n),y(n)
    real,dimension(n*n) :: bigz
    integer :: i,j,k
    k=1
    do i=1,n
do j=1,n
   bigz(k) = (((x(i)-x(j))**2 + (y(i)-y(j))**2))**0.5
   k = k+1
           enddo
           enddo



      end subroutine bigmat

  attributes(global) subroutine bigmatb(n,x, bigz)
    implicit none
    integer :: n
    real :: x(n)
    integer :: i,j,k

    real,dimension(n*n) :: bigz
    k = 1
    do i=1,n
       do j=1,n
          bigz(k) = (((x(i)-x(j))**2))**0.5
          k = k+1
           enddo
           enddo
           


      end subroutine bigmatb

  attributes(global) subroutine bigmatc(n,x,u,w, bigz,sta)
    implicit none
    integer :: n
    real :: x(n),y(n),sta,h
    real :: u(n),w(n)
    integer :: i,j,k
    real,dimension(n*n) :: bigz
    k = 1
    do i=1,n
       do j=1,n
              h = (((x(i)-x(j))**2)*sta) + ((u(i)-u(j))**2 + (w(i)-w(j))**2) 
              bigz(k) = h**0.5
              k = k+1
           enddo
           enddo
           


      end subroutine bigmatc


      
      attributes(global) subroutine type1(sill_d,rg_d,nug_d,mat1,bigsum_d)
        implicit none
        real :: sill_d,rg_d,nug_d
        integer :: i,j,k,n
        real :: mat1(:,:)
        real :: bigsum_d(:)
        n = size(mat1,1)
        k=1
         do i=1,n
            do j=1,n
             bigsum_d(k)=(sill_d)*(1 - exp(-mat1(i,j)/rg_d)) + nug_d
             k=k+1
             enddo
          enddo
        end subroutine type1

        attributes(global) subroutine type2(sill_d,rg_d,nug_d,mat1,bigsum)
          implicit none
        real :: sill_d,rg_d,nug_d
        integer :: i,j,k,n
        real :: mat1(:,:)
        real :: bigsum(:)
        n = size(mat1,1)
        k=1
         do i=1,n
            do j=1,n
             bigsum(k)= (sill_d)*(1 - exp(-(mat1(i,j)/rg_d)**2)) + nug_d   
             k=k+1
             enddo
          enddo
        end subroutine type2

        attributes(global) subroutine type3(sill_d,rg_d,nug_d,mat1,bigsum)
          implicit none
          real :: sill_d,rg_d,nug_d
          real :: xx,xy
        integer :: i,j,k,n
        real :: mat1(:,:)
        real :: bigsum(:)
        n = size(mat1,1)  
          k = 1
               do i=1,n
                  do j=1,n
                   if(mat1(i,j).gt.0)then
                   if(mat1(i,j).le.rg_d)then
                      xx = mat1(i,j)/rg_d
                      xy = 1.5*xx - 0.5*(xx**3) 
                      bigsum(k)= xy*(sill_d) + nug_d
                      
                      else
                            bigsum(k)=1.0*(sill_d) + nug_d
                      endif
                      else
                         bigsum(k) = 0.0
                   endif

                      k=k+1
                      enddo
                 enddo
          end subroutine type3
        
        
end module mathOps

subroutine warm4(bigsum,n1,n2)
     use mathOps
        use cudafor !has to go first        
        implicit none
        integer :: n1, n2
        integer :: i,j,k,n3,n4,n5,n4e
        character(len=3) :: st_s,st_t,st_st
        integer, device :: n1_d
        real,allocatable :: x(:),y(:),big1(:),big2(:),big3(:)
        real,allocatable :: xa(:),u(:),w(:)
        real, dimension(n2) :: bigsum,bigz,bigy,bigx
        real, dimension(n1,n1) :: mat1,mata,matb
        real, dimension(n2) :: bigsuma,bigsumb
        real, device, dimension(n1,n1) :: mat2
        real, device, dimension(n1) :: x_d, y_d, u_d,w_d
        real, device :: sill_d,rg_d,nug_d,st_d
        real, device, dimension(n2) :: bigsum_d,bigz_d
        real :: xx,xy,nug,v1
        real :: sill_s,rg_s,nug_s,sill_t,rg_t,nug_t
        real :: sill_st,rg_st,nug_st,stAni

        

        open (unit=1,file="file1.txt")
        read(1,200)st_s
        read(1,200)st_t
        read(1,200)st_st
        read(1,*)sill_s
        read(1,*)rg_s
        read(1,*)nug_s
        read(1,*)sill_t
        read(1,*)rg_t
        read(1,*)nug_t
        read(1,*)sill_st
        read(1,*)rg_st
        read(1,*)nug_st
        read(1,*)stAni
        close(unit=1)


200     format(a)
225     format(f10.2)
250     format(i5)

        allocate(x(n1),y(n1))
        open(unit=2,file="file2.txt")
        do i=1,n1
           read(2,*)x(i),y(i)
           enddo
        close(unit=2)   


        
        
        x_d = x
        y_d = y
        n1_d = n1
        bigz_d = 0
        call bigmat<<<512,4>>>(n1_d,x_d,y_d,bigz_d)
        bigz = bigz_d


        deallocate(x,y)


        
           sill_d = sill_s
           rg_d = rg_s
           nug_d = nug_s
           bigsum_d = 0.0
           mat2 = reshape(bigz,(/n1_d,n1_d/))
        
        if(st_s == "Exp")then
           call type1<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
           bigsum = bigsum_d
           
          else if(st_s == "Gau") then
          call type2<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
          bigsum = bigsum_d

        else if(st_s == "Sph") then


           call type3<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
           bigsum = bigsum_d

              endif
              
        write(*,*)"more",bigsum(1:10)
              
        allocate(xa(n1))
        open(unit=2,file="file11.txt")
        do i=1,n1
           read(2,*)xa(i)
           enddo
        close(unit=2)   

        x_d = xa
        n1_d = n1
        bigz_d = 0
        call bigmatb<<<512,4>>>(n1_d,x_d,bigz_d)
        bigy = bigz_d




           sill_d = sill_s
           rg_d = rg_s
           nug_d = nug_s
           bigsum_d = 0.0
           mat2 = reshape(bigy,(/n1_d,n1_d/))           

        
        if(st_s == "Exp")then
           call type1<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
           bigsuma = bigsum_d
           
          else if(st_s == "Gau") then
          call type2<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
          bigsuma = bigsum_d

        else if(st_s == "Sph") then

          
           call type3<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
           bigsuma = bigsum_d

              endif
              
              
              bigsum = bigsum + bigsuma
              write(*,*)bigsum(1:10)
              write(*,*)bigsuma(1:10)

        
        allocate(u(n1),w(n1))
        open(unit=8,file="file2.txt")
        do i=1,n1
           read(8,*)u(i),w(i)
           enddo
        close(unit=8)   

        x_d = xa
        n1_d = n1
        st_d = stAni
        u_d = u
        w_d = w
        bigz_d = 0.0
        call bigmatc<<<512,4>>>(n1_d,x_d,u_d,w_d,bigz_d,st_d)
        bigx = bigz_d

           sill_d = sill_s
           rg_d = rg_s
           nug_d = nug_s
           bigsum_d = 0.0
           mat2 = reshape(bigx,(/n1,n1/))
        
        if(st_s == "Exp")then
           call type1<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
           bigsumb = bigsum_d

          else if(st_s == "Gau") then
          call type2<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
          bigsumb = bigsum_d

        else if(st_s == "Sph") then



           call type3<<<512,4>>>(sill_d,rg_d,nug_d,mat2,bigsum_d)
           
           bigsumb = bigsum_d

        endif
        write(*,*)bigsumb(1:10)
        bigsum = bigsum + bigsumb
        


      end subroutine warm4
      
            


