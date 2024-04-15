subroutine cus2a(m,A_a,B_a)
!Compilation example: 'nvfortran cus2.F90 -o cus2 -cudalib=cusolver -cuda'
use cudafor !has to go first
use cusolverDn
    implicit none
integer :: info
 !   integer,parameter :: q2 = SELECTED_REAL_KIND(15,305)
    real, device, dimension(m,m) :: A_d,B_d
    real, dimension(m,m) :: A,B
    real, dimension(m*m) :: A_a, B_a
    !    real(q2), device, dimension(3) :: W_d
    integer, device, dimension(m) :: W_d
    real, dimension(m) :: W
    integer :: stat, lwork, m, lda,bs,ldb
    real, device, allocatable  :: work_d(:)
    integer, device :: devInfo, P
    type(cusolverDnHandle) :: h
    stat=cusolverDnCreate(h)
        W_d=0.0

        A_d = reshape(A_a,(/m,m/))
        B_d = reshape(B_a,(/m,m/))


    lda = m
    ldb = m
    bs = 1
!    A_d(1,1:3)=(/4,1,2/)
!    A_d(2,1:3)=(/1,-1,1/)
!    A_d(3,1:3)=(/2,1,3/)    !eigenvalues are 5.84947, 1.44865, -1.29812
!   A_d(1,1:3)=(/1,0,0/)
!   A_d(2,1:3)=(/0,1,0/)
!   A_d(3,1:3)=(/0,0,1/)
    stat=cusolverDnSgetrf_bufferSize(h, m,  m, A_d, lda, lwork)

    allocate(work_d(lwork))
    stat=cusolverDnSgetrf(h, m, m, A_d, lda, work_d, W_d, devInfo)



stat=cusolverDnSgetrs(h,CUBLAS_OP_N,m,m,A_d,lda,W_d,B_d,lda,devInfo)

info=devInfo
    stat=cudaDeviceSynchronize()

    A=A_d
    B=B_d

    
    A_a = reshape(A,(/m*m/))
    B_a = reshape(B,(/m*m/))

    
    deallocate(work_d)
    stat=cusolverDnDestroy(h)
    
end subroutine cus2a
