program Main

use lapack

real :: A(10,10), B(10,10), C(5,5), L(5,5), U(5,5)
real, allocatable :: work(:)
integer :: ipiv(5), info, flag, lwork

call random_number(A)

B = A

call sgetrf(5, 5, A(6,1), 10, ipiv, info)

call slaset('L', 5, 5, 0.0, 1.0, U, 5)
call slacpy('L', 5, 5, A(6,1), 10, L, 5)
call slacpy('U', 5, 5, A(6,1), 10, U, 5)
call slaset('U', 5, 5, 0.0, 1.0, L, 5)

call sgemm('N','N', 5, 5, 5, 1.0, L, 5, U, 5, 0.0, C, 5)
call slaswp(5, C, 5, 1, 5, ipiv, -1)

print '(5f5.2)', transpose(C)
print *
print '(5f5.2)', transpose(B(6:10,1:5))
print *

allocate(work(1), STAT=flag)
call sgetri(5, A(6,1), 10, ipiv, work, -1, info)

lwork = int(work(1))
allocate(work(lwork), STAT=flag)
call sgetri(5, A(6,1), 10, ipiv, work, lwork, info)

deallocate(work)

call sgemm('N','N',5,5,5,1.0,B(6,1),10,A(6,1),10,0.0,L,5)

print '(5f5.2)',transpose(L)

end program