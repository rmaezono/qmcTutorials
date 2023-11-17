module slaarnabt
use dsp
use slaarnaag,only : czero,pi,min_exp,max_exp,minimum_exp_arg,maximum_&
&exp_arg,largest_representable_number,smallest_representable_number,ma&
&x_rep_int,maximum_exp_arg_actual,minimum_exp_arg_actual,third,one_ove&
&r_pi
use format_utils,only : i2s
use run_control, only : errstop,errwarn,check_alloc
implicit none
private
public daxpy,dcopy,ddot,dnrm2,dscal,zaxpy,zscal,zdotu,zcopy,dgemv,zgem&
&v,dgemm,dsymm,dsymv
public dcopy3,dcopy6,dcopy_ee,dcopy_ee_flip,dsum,dsum3,dadd,zadd,ddot_&
&s,ddot_alt,dmatmul,dswap,iswap,swap1,dzscal,mxmb,mxmbc,multi_ddot_s,m&
&ulti_zdotu_s,zdotu_alt,zdot_rc,d_sparse_matrix_vector,d_sparse_vector&
&_matrix,d2_sparse_matrix_vector,d2_sparse_vector_matrix
public zdotu_cc,zaxpy_cc,zscal_cc,zgemm_cc,zgemv_cc,zger_cc
public d_or_z_gemm,d_or_z_gemv,d_or_z_ger,d_or_z_scal
public sv_decom,svd_inverse,inverse,eigenproblem,lu_decom,lu_solve_onc&
&e,lu_solve_n,lu_decom_cmplx_dz,lu_solve_once_cmplx_dz,lu_solve_n_cmpl&
&x_dz,lu_decom_cmplx,lu_solve_once_cmplx,lu_solve_n_cmplx,gmult,g3mult&
&,g6mult
public overflow_protection,exp_protect,exp_limit
public quicksort,partial_rank_int,sort_matrix_rect,sort_matrix_symm,am&
&b2cand_sort_matrix_rect,amb2cand_sort_matrix_symm
public correlation_time,reaverage,median,correlation_time_alt
public erfc,bessel_j0,bessel_j1,sph_bessel,chebyshev,exponential
public interp_nev,interp_nev_with_derivs,lookup,interp_nev_2d_with_der&
&ivs
public lu_logdet,lu_logdet_cmplx_dz,lu_logdet_cmplx
public resize_pointer
public param_limiting,param_unlimiting,choose,reduced_echelon,gauss_qu&
&adrature,make_representable,hypothenuse,next_permutation,brentmin,gsb&
&racket,bracket,inverse3,parabolic_min,exp_int,gamma_ser,get_numerical&
&_orbmask,get_numerical_orbrange,isdiag_int,approx_equal
interface
subroutine daxpy(n,da,dx,incx,dy,incy)
implicit none
integer,intent(in) :: incx,incy,n
real(kind(1.d0)),intent(in) :: da
real(kind(1.d0)),intent(in) :: dx(*)
real(kind(1.d0)),intent(inout) :: dy(*)
end subroutine daxpy
subroutine dcopy(n,dx,incx,dy,incy)
implicit none
integer,intent(in) :: incx,incy,n
real(kind(1.d0)),intent(in) :: dx(*)
real(kind(1.d0)),intent(out) :: dy(*)
end subroutine dcopy
subroutine dscal(n,da,dx,incx)
implicit none
integer,intent(in) :: incx,n
real(kind(1.d0)),intent(in) :: da
real(kind(1.d0)),intent(inout) :: dx(*)
end subroutine dscal
subroutine dswap(n,dx,incx,dy,incy)
implicit none
integer,intent(in) :: incx,incy,n
real(kind(1.d0)),intent(inout) :: dx(*),dy(*)
end subroutine dswap
real(kind(1.d0)) function ddot(n,dx,incx,dy,incy)
implicit none
integer,intent(in) :: incx,incy,n
real(kind(1.d0)),intent(in) :: dx(*),dy(*)
end function ddot
real(kind(1.d0)) function dnrm2(n,x,incx)
implicit none
integer,intent(in) :: incx,n
real(kind(1.d0)),intent(in) :: x(*)
end function dnrm2
subroutine zaxpy(n,za,zx,incx,zy,incy)
implicit none
integer,intent(in) :: incx,incy,n
complex(kind(1.d0)),intent(in) :: za
complex(kind(1.d0)),intent(in) :: zx(*)
complex(kind(1.d0)),intent(inout) :: zy(*)
end subroutine zaxpy
subroutine zcopy(n,zx,incx,zy,incy)
implicit none
integer,intent(in) :: incx,incy,n
complex(kind(1.d0)),intent(in) :: zx(*)
complex(kind(1.d0)),intent(out) :: zy(*)
end subroutine zcopy
subroutine zscal(n,za,zx,incx)
implicit none
integer,intent(in) :: incx,n
complex(kind(1.d0)),intent(in) :: za
complex(kind(1.d0)),intent(inout) :: zx(*)
end subroutine zscal
complex(kind(1.d0)) function zdotu(n,zx,incx,zy,incy)
implicit none
integer,intent(in) :: incx,incy,n
complex(kind(1.d0)),intent(in) :: zx(*),zy(*)
end function zdotu
subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
implicit none
integer,intent(in) :: k,lda,ldb,ldc,m,n
real(kind(1.d0)),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
character,intent(in) :: transa,transb
real(kind(1.d0)),intent(inout) :: c(ldc,*)
end subroutine dgemm
subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
implicit none
integer,intent(in) :: incx,incy,lda,m,n
real(kind(1.d0)),intent(in) :: alpha,beta,a(lda,*),x(*)
character,intent(in) :: trans
real(kind(1.d0)),intent(inout) :: y(*)
end subroutine dgemv
subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
implicit none
integer,intent(in) :: incx,incy,lda,m,n
complex(kind(1.d0)),intent(in) :: alpha,beta,a(lda,*),x(*)
character,intent(in) :: trans
complex(kind(1.d0)),intent(inout) :: y(*)
end subroutine zgemv
subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
implicit none
real(kind(1.d0)),intent(in) :: alpha,beta,x(*)
integer,intent(in) :: incx,incy,lda,n
character,intent(in) :: uplo
real(kind(1.d0)),intent(inout) :: a(lda,*),y(*)
end subroutine dsymv
subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
implicit none
integer,intent(in) :: lda,ldb,ldc,m,n
real(kind(1.d0)),intent(in) :: alpha,beta,b(ldb,*)
character,intent(in) :: side,uplo
real(kind(1.d0)),intent(inout) :: a(lda,*),c(ldc,*)
end subroutine dsymm
end interface
interface
subroutine dgetrf(m,n,a,lda,ipiv,info)
implicit none
integer,intent(in) :: m,n,lda
real(kind(1.d0)),intent(inout) :: a(lda,*)
integer,intent(out) :: info,ipiv(*)
end subroutine dgetrf
subroutine dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
implicit none
integer,intent(in) :: lda,ldb,n,nrhs,ipiv(*)
real(kind(1.d0)),intent(in) :: a(lda,*)
character,intent(in) :: trans
real(kind(1.d0)),intent(inout) :: b(ldb,*)
integer,intent(out) :: info
end subroutine dgetrs
subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,&
&info)
implicit none
integer,intent(in) :: lda,ldvl,ldvr,lwork,n
character,intent(in) :: jobvl,jobvr
real(kind(1.d0)),intent(inout) :: a(lda,*),work(*)
integer,intent(out) :: info
real(kind(1.d0)),intent(out) :: vl(ldvl,*),vr(ldvr,*),wi(*),wr(*)
end subroutine dgeev
subroutine dgesdd(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info&
&)
implicit none
integer,intent(in) :: lda,ldu,ldvt,lwork,m,n
character,intent(in) :: jobz
integer,intent(inout) :: iwork(*)
real(kind(1.d0)),intent(inout) :: a(lda,*),work(*)
integer,intent(out) :: info
real(kind(1.d0)),intent(out) :: s(*),u(ldu,*),vt(ldvt,*)
end subroutine dgesdd
subroutine zgetrf(m,n,a,lda,ipiv,info)
implicit none
integer,intent(in) :: lda,m,n
complex(kind(1.d0)),intent(inout) :: a(lda,*)
integer,intent(out) :: info,ipiv(*)
end subroutine zgetrf
subroutine zgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
implicit none
integer,intent(in) :: lda,ldb,n,nrhs,ipiv(*)
complex(kind(1.d0)),intent(in) :: a(lda,*)
character,intent(in) :: trans
complex(kind(1.d0)),intent(inout) :: b(ldb,*)
integer,intent(out) :: info
end subroutine zgetrs
end interface
interface swap1
module procedure xyzzyaaah1,xyzzyaaaf1,xyzzyaaag1,xyzzyaaai1
end interface swap1
interface resize_pointer
module procedure resize_pointer_dble1,resize_pointer_dble2,resize_poin&
&ter_int1,resize_pointer_int2,resize_pointer_int3,resize_pointer_bool1&
&,resize_pointer_bool2,resize_pointer_bool3,resize_pointer_char1
end interface resize_pointer
interface approx_equal
module procedure approx_equal_sp,approx_equal_dp
end interface approx_equal
contains
subroutine overflow_protection
implicit none
largest_representable_number=huge(1.d0)
smallest_representable_number=tiny(1.d0)
max_rep_int=dble(huge(1))
maximum_exp_arg_actual=log(largest_representable_number)
minimum_exp_arg_actual=log(smallest_representable_number)
maximum_exp_arg=log(max_exp)
minimum_exp_arg=log(min_exp)
if(maximum_exp_arg>=maximum_exp_arg_actual)call errstop('OVERFLOW_PROT&
&ECTION','Largest representable double precision number on this machin&
&e is too small for the inbuilt exp() overflow protection to work.')
if(minimum_exp_arg<=minimum_exp_arg_actual)call errstop('OVERFLOW_PROT&
&ECTION','Smallest representable double precision number on this machi&
&ne is too large for the inbuilt exp() underflow protection to work.')
end subroutine overflow_protection
real(dp) function exp_protect(exparg)
implicit none
real(dp),intent(in) :: exparg
if(exparg<minimum_exp_arg)then
exp_protect=0.d0
elseif(exparg<maximum_exp_arg)then
exp_protect=exp(exparg)
else
exp_protect=max_exp
endif
end function exp_protect
real(dp) function exp_limit(exparg)
implicit none
real(dp),intent(in) :: exparg
if(exparg<minimum_exp_arg)then
exp_limit=min_exp
elseif(exparg<maximum_exp_arg)then
exp_limit=exp(exparg)
else
exp_limit=max_exp
endif
end function exp_limit
function inverse3(a,is_singular,determinant)
implicit none
real(dp),intent(in) :: a(3,3)
real(dp),intent(out),optional :: determinant
logical,intent(out),optional :: is_singular
real(dp) xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29
real(dp) :: inverse3(3,3)
real(dp),parameter :: xyzzyaaaf29=1.d-16
xyzzyaaaa29=a(2,2)*a(3,3)-a(2,3)*a(3,2)
xyzzyaaab29=a(2,3)*a(3,1)-a(2,1)*a(3,3)
xyzzyaaac29=a(2,1)*a(3,2)-a(2,2)*a(3,1)
xyzzyaaad29=a(1,1)*xyzzyaaaa29+a(1,2)*xyzzyaaab29+a(1,3)*xyzzyaaac29
if(present(is_singular))is_singular=(abs(xyzzyaaad29)<xyzzyaaaf29)
if(present(determinant))determinant=xyzzyaaad29
if(xyzzyaaad29/=0.d0)then
xyzzyaaae29=1.d0/xyzzyaaad29
inverse3(1,1)=xyzzyaaaa29*xyzzyaaae29
inverse3(2,1)=xyzzyaaab29*xyzzyaaae29
inverse3(3,1)=xyzzyaaac29*xyzzyaaae29
xyzzyaaaa29=a(1,1)*xyzzyaaae29
xyzzyaaab29=a(1,2)*xyzzyaaae29
xyzzyaaac29=a(1,3)*xyzzyaaae29
inverse3(1,2)=xyzzyaaac29*a(3,2)-xyzzyaaab29*a(3,3)
inverse3(2,2)=xyzzyaaaa29*a(3,3)-xyzzyaaac29*a(3,1)
inverse3(3,2)=xyzzyaaab29*a(3,1)-xyzzyaaaa29*a(3,2)
inverse3(1,3)=xyzzyaaab29*a(2,3)-xyzzyaaac29*a(2,2)
inverse3(2,3)=xyzzyaaac29*a(2,1)-xyzzyaaaa29*a(2,3)
inverse3(3,3)=xyzzyaaaa29*a(2,2)-xyzzyaaab29*a(2,1)
elseif(.not.(present(is_singular).or.present(determinant)))then
call errstop("INVERSE3","Singular matrix found")
endif
end function inverse3
subroutine quicksort(n,x,indx)
implicit none
integer,intent(in) :: n
integer,intent(inout) :: indx(n)
real(dp),intent(in) :: x(n)
integer,parameter :: xyzzyaaaa30=30
integer,parameter :: xyzzyaaab30=128
integer xyzzyaaac30,xyzzyaaad30,xyzzyaaae30,xyzzyaaaf30,xyzzyaaag30,xy&
&zzyaaah30(xyzzyaaab30),xyzzyaaai30(xyzzyaaab30),xyzzyaaaj30
real(dp) xyzzyaaak30
do xyzzyaaae30=1,n
indx(xyzzyaaae30)=xyzzyaaae30
enddo
if(n>xyzzyaaaa30)then
xyzzyaaaj30=0
xyzzyaaac30=1
xyzzyaaad30=n
do
xyzzyaaae30=xyzzyaaac30
xyzzyaaaf30=xyzzyaaad30
xyzzyaaag30=(xyzzyaaac30+xyzzyaaad30)/2
if(x(indx(xyzzyaaac30))>x(indx(xyzzyaaag30)))call swap1(indx(xyzzyaaag&
&30),indx(xyzzyaaac30))
if(x(indx(xyzzyaaag30))>x(indx(xyzzyaaad30)))call swap1(indx(xyzzyaaag&
&30),indx(xyzzyaaad30))
if(x(indx(xyzzyaaac30))>x(indx(xyzzyaaag30)))call swap1(indx(xyzzyaaag&
&30),indx(xyzzyaaac30))
xyzzyaaak30=x(indx(xyzzyaaag30))
do
xyzzyaaae30=xyzzyaaae30+1
do while(x(indx(xyzzyaaae30))<xyzzyaaak30)
xyzzyaaae30=xyzzyaaae30+1
enddo
xyzzyaaaf30=xyzzyaaaf30-1
do while(x(indx(xyzzyaaaf30))>xyzzyaaak30)
xyzzyaaaf30=xyzzyaaaf30-1
enddo
if(xyzzyaaae30>=xyzzyaaaf30)exit
call swap1(indx(xyzzyaaae30),indx(xyzzyaaaf30))
enddo
if(xyzzyaaad30-xyzzyaaaf30>=xyzzyaaae30-xyzzyaaac30.and.xyzzyaaae30-xy&
&zzyaaac30>xyzzyaaaa30)then
xyzzyaaaj30=xyzzyaaaj30+1
if(xyzzyaaaj30>xyzzyaaab30)call errstop('QUICKSORT','Stack size too sm&
&all <1>.')
xyzzyaaah30(xyzzyaaaj30)=xyzzyaaaf30+1
xyzzyaaai30(xyzzyaaaj30)=xyzzyaaad30
xyzzyaaad30=xyzzyaaae30-1
elseif(xyzzyaaae30-xyzzyaaac30>xyzzyaaad30-xyzzyaaaf30.and.xyzzyaaad30&
&-xyzzyaaaf30>xyzzyaaaa30)then
xyzzyaaaj30=xyzzyaaaj30+1
if(xyzzyaaaj30>xyzzyaaab30)call errstop('QUICKSORT','Stack size too sm&
&all <2>.')
xyzzyaaah30(xyzzyaaaj30)=xyzzyaaac30
xyzzyaaai30(xyzzyaaaj30)=xyzzyaaae30-1
xyzzyaaac30=xyzzyaaaf30+1
elseif(xyzzyaaad30-xyzzyaaaf30>xyzzyaaaa30)then
xyzzyaaac30=xyzzyaaaf30+1
elseif(xyzzyaaae30-xyzzyaaac30>xyzzyaaaa30)then
xyzzyaaad30=xyzzyaaae30-1
else
if(xyzzyaaaj30<1)exit
xyzzyaaac30=xyzzyaaah30(xyzzyaaaj30)
xyzzyaaad30=xyzzyaaai30(xyzzyaaaj30)
xyzzyaaaj30=xyzzyaaaj30-1
endif
enddo
endif
do xyzzyaaae30=2,n
if(x(indx(xyzzyaaae30-1))<=x(indx(xyzzyaaae30)))cycle
call swap1(indx(xyzzyaaae30-1),indx(xyzzyaaae30))
do xyzzyaaaf30=xyzzyaaae30-1,2,-1
if(x(indx(xyzzyaaaf30-1))<=x(indx(xyzzyaaaf30)))exit
call swap1(indx(xyzzyaaaf30-1),indx(xyzzyaaaf30))
enddo
enddo
end subroutine quicksort
subroutine sort_matrix_rect(n,m,ns,matrix,row_indx,col_indx,namb_range&
&,amb_range)
implicit none
integer,intent(in) :: n,m,matrix(n,m)
integer,intent(inout) :: ns,row_indx(n,*),col_indx(m,*),namb_range,amb&
&_range(2,n)
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31,xyzzyaaad31,xyzzyaaae31,xy&
&zzyaaaf31,xyzzyaaag31,xyzzyaaah31,xyzzyaaai31(n),xyzzyaaaj31(n),xyzzy&
&aaak31(m),xyzzyaaal31,xyzzyaaam31(m),xyzzyaaan31,xyzzyaaao31(2,n),xyz&
&zyaaap31,xyzzyaaaq31,xyzzyaaar31,xyzzyaaas31,xyzzyaaat31,xyzzyaaau31,&
&xyzzyaaav31,xyzzyaaaw31
logical xyzzyaaax31,xyzzyaaay31,xyzzyaaaz31,xyzzyaaba31,xyzzyaabb31,xy&
&zzyaabc31
integer xyzzyaabd31(2),xyzzyaabe31(2,n,2),xyzzyaabf31
integer,parameter :: xyzzyaabg31=1024
integer xyzzyaabh31(2),xyzzyaabi31(n,xyzzyaabg31,2),xyzzyaabj31(m,xyzz&
&yaabg31,2),xyzzyaabk31
xyzzyaaaa31=2
xyzzyaaab31=1
xyzzyaabh31(xyzzyaaab31)=ns
xyzzyaabi31(1:n,1:ns,xyzzyaaab31)=row_indx(1:n,1:ns)
xyzzyaabj31(1:m,1:ns,xyzzyaaab31)=col_indx(1:m,1:ns)
xyzzyaabd31(xyzzyaaab31)=0
do xyzzyaabf31=1,namb_range
if(amb_range(2,xyzzyaabf31)>amb_range(1,xyzzyaabf31))then
xyzzyaabd31(xyzzyaaab31)=xyzzyaabd31(xyzzyaaab31)+1
xyzzyaabe31(1:2,xyzzyaabd31(xyzzyaaab31),xyzzyaaab31)=amb_range(1:2,xy&
&zzyaabf31)
endif
enddo
do xyzzyaaap31=1,m
call swap1(xyzzyaaaa31,xyzzyaaab31)
xyzzyaaax31=.false.
xyzzyaabh31(xyzzyaaab31)=0
do xyzzyaabk31=1,xyzzyaabh31(xyzzyaaaa31)
xyzzyaaak31(1:m)=xyzzyaabj31(1:m,xyzzyaabk31,xyzzyaaaa31)
xyzzyaaar31=m-xyzzyaaap31+1
xyzzyaaan31=0
do xyzzyaabf31=0,xyzzyaabd31(xyzzyaaaa31)
if(xyzzyaabf31>0)then
xyzzyaaag31=xyzzyaabe31(1,xyzzyaabf31,xyzzyaaaa31)
xyzzyaaah31=xyzzyaabe31(2,xyzzyaabf31,xyzzyaaaa31)
xyzzyaaas31=0
xyzzyaaba31=.false.
xyzzyaaaw31=0
do
xyzzyaabb31=.false.
xyzzyaaat31=0
xyzzyaaav31=0
xyzzyaaal31=0
do xyzzyaaae31=xyzzyaaap31,xyzzyaaap31+xyzzyaaar31-1
xyzzyaabc31=.false.
xyzzyaaau31=0
do xyzzyaaac31=xyzzyaaag31,xyzzyaaah31
xyzzyaaaf31=matrix(xyzzyaabi31(xyzzyaaac31,xyzzyaabk31,xyzzyaaaa31),xy&
&zzyaaak31(xyzzyaaae31))
if(xyzzyaaba31.and.xyzzyaaaf31<=xyzzyaaaw31)cycle
if(xyzzyaabb31.and.xyzzyaaaf31>xyzzyaaat31)cycle
if(.not.xyzzyaabb31.or.xyzzyaaaf31<xyzzyaaat31)then
xyzzyaabc31=.true.
xyzzyaaat31=xyzzyaaaf31
xyzzyaaau31=1
else
xyzzyaaau31=xyzzyaaau31+1
endif
xyzzyaabb31=.true.
enddo
if(xyzzyaabc31.or.xyzzyaaau31>xyzzyaaav31)then
xyzzyaaal31=1
xyzzyaaam31(1)=xyzzyaaae31
xyzzyaaav31=xyzzyaaau31
elseif(xyzzyaaau31==xyzzyaaav31)then
xyzzyaaal31=xyzzyaaal31+1
xyzzyaaam31(xyzzyaaal31)=xyzzyaaae31
endif
enddo
xyzzyaaba31=.true.
xyzzyaaaw31=xyzzyaaat31
if(xyzzyaaal31<xyzzyaaar31)then
do xyzzyaaae31=1,xyzzyaaal31
call swap1(xyzzyaaak31(xyzzyaaap31+xyzzyaaae31-1),xyzzyaaak31(xyzzyaaa&
&m31(xyzzyaaae31)))
enddo
xyzzyaaar31=xyzzyaaal31
endif
xyzzyaaaj31(xyzzyaaag31+xyzzyaaas31:xyzzyaaag31+xyzzyaaas31+xyzzyaaav3&
&1-1)=xyzzyaaat31
if(xyzzyaaav31>1)then
xyzzyaaan31=xyzzyaaan31+1
xyzzyaaao31(1:2,xyzzyaaan31)=(/xyzzyaaag31+xyzzyaaas31,xyzzyaaag31+xyz&
&zyaaas31+xyzzyaaav31-1/)
endif
xyzzyaaas31=xyzzyaaas31+xyzzyaaav31
if(xyzzyaaag31+xyzzyaaas31>xyzzyaaah31)exit
enddo
endif
xyzzyaaag31=1
if(xyzzyaabf31>0)xyzzyaaag31=xyzzyaabe31(2,xyzzyaabf31,xyzzyaaaa31)+1
xyzzyaaah31=n
if(xyzzyaabf31<xyzzyaabd31(xyzzyaaaa31))xyzzyaaah31=xyzzyaabe31(1,xyzz&
&yaabf31+1,xyzzyaaaa31)-1
if(xyzzyaaah31>=xyzzyaaag31)then
do xyzzyaaac31=xyzzyaaag31,xyzzyaaah31
if(xyzzyaaar31==1)then
do xyzzyaaad31=xyzzyaaac31,xyzzyaaah31
xyzzyaaaj31(xyzzyaaad31)=matrix(xyzzyaabi31(xyzzyaaad31,xyzzyaabk31,xy&
&zzyaaaa31),xyzzyaaak31(xyzzyaaap31))
enddo
exit
endif
xyzzyaabb31=.false.
xyzzyaaat31=0
xyzzyaaau31=0
xyzzyaaal31=0
do xyzzyaaae31=xyzzyaaap31,xyzzyaaap31+xyzzyaaar31-1
xyzzyaaaf31=matrix(xyzzyaabi31(xyzzyaaac31,xyzzyaabk31,xyzzyaaaa31),xy&
&zzyaaak31(xyzzyaaae31))
if(xyzzyaabb31.and.xyzzyaaaf31>xyzzyaaat31)cycle
if(.not.xyzzyaabb31.or.xyzzyaaaf31<xyzzyaaat31)then
xyzzyaaat31=xyzzyaaaf31
xyzzyaaal31=1
xyzzyaaam31(1)=xyzzyaaae31
else
xyzzyaaau31=xyzzyaaau31+1
xyzzyaaal31=xyzzyaaal31+1
xyzzyaaam31(xyzzyaaal31)=xyzzyaaae31
endif
xyzzyaabb31=.true.
enddo
xyzzyaaaj31(xyzzyaaac31)=xyzzyaaat31
if(xyzzyaaal31<xyzzyaaar31)then
do xyzzyaaae31=1,xyzzyaaal31
call swap1(xyzzyaaak31(xyzzyaaap31+xyzzyaaae31-1),xyzzyaaak31(xyzzyaaa&
&m31(xyzzyaaae31)))
enddo
xyzzyaaar31=xyzzyaaal31
endif
enddo
endif
enddo
xyzzyaaay31=.false.
xyzzyaaaz31=.false.
if(.not.xyzzyaaax31)then
xyzzyaaay31=.true.
else
do xyzzyaaac31=1,n
if(xyzzyaaai31(xyzzyaaac31)<xyzzyaaaj31(xyzzyaaac31))then
xyzzyaaaz31=.true.
exit
elseif(xyzzyaaai31(xyzzyaaac31)>xyzzyaaaj31(xyzzyaaac31))then
xyzzyaaay31=.true.
exit
endif
enddo
endif
if(xyzzyaaaz31)cycle
if(xyzzyaaay31)then
xyzzyaaai31=xyzzyaaaj31
xyzzyaabd31(xyzzyaaab31)=xyzzyaaan31
if(xyzzyaaan31>0)xyzzyaabe31(1:2,1:xyzzyaaan31,xyzzyaaab31)=xyzzyaaao3&
&1(1:2,1:xyzzyaaan31)
xyzzyaabh31(xyzzyaaab31)=0
xyzzyaaax31=.true.
endif
do xyzzyaaae31=xyzzyaaap31,xyzzyaaap31+xyzzyaaar31-1
xyzzyaabh31(xyzzyaaab31)=xyzzyaabh31(xyzzyaaab31)+1
if(xyzzyaabh31(xyzzyaaab31)>xyzzyaabg31)call errstop('SORT_MATRIX_RECT&
&','Size of buffer too small.')
xyzzyaabi31(1:n,xyzzyaabh31(xyzzyaaab31),xyzzyaaab31)=xyzzyaabi31(1:n,&
&xyzzyaabk31,xyzzyaaaa31)
xyzzyaabj31(1:m,xyzzyaabh31(xyzzyaaab31),xyzzyaaab31)=xyzzyaaak31(1:m)
call swap1(xyzzyaabj31(xyzzyaaap31,xyzzyaabh31(xyzzyaaab31),xyzzyaaab3&
&1),xyzzyaabj31(xyzzyaaae31,xyzzyaabh31(xyzzyaaab31),xyzzyaaab31))
xyzzyaaaq31=xyzzyaabj31(xyzzyaaap31,xyzzyaabh31(xyzzyaaab31),xyzzyaaab&
&31)
do xyzzyaabf31=1,xyzzyaabd31(xyzzyaaaa31)
call xyzzyaaaa1(n,xyzzyaabe31(1,xyzzyaabf31,xyzzyaaaa31),xyzzyaabe31(2&
&,xyzzyaabf31,xyzzyaaaa31),matrix(1,xyzzyaaaq31),xyzzyaabi31(1,xyzzyaa&
&bh31(xyzzyaaab31),xyzzyaaab31))
enddo
enddo
enddo
enddo
ns=xyzzyaabh31(xyzzyaaab31)
row_indx(1:n,1:ns)=xyzzyaabi31(1:n,1:ns,xyzzyaaab31)
col_indx(1:m,1:ns)=xyzzyaabj31(1:m,1:ns,xyzzyaaab31)
namb_range=xyzzyaabd31(xyzzyaaab31)
if(namb_range>0)amb_range(1:2,1:namb_range)=xyzzyaabe31(1:2,1:namb_ran&
&ge,xyzzyaaab31)
end subroutine sort_matrix_rect
subroutine sort_matrix_symm(n,ns,matrix,indx,namb_range,amb_range)
implicit none
integer,intent(in) :: n,matrix(n,n)
integer,intent(inout) :: ns,indx(n,*),namb_range,amb_range(2,n)
integer xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32,xyzzyaaad32,xyzzyaaae32,xy&
&zzyaaaf32(n),xyzzyaaag32,xyzzyaaah32,xyzzyaaai32,xyzzyaaaj32,xyzzyaaa&
&k32,xyzzyaaal32,xyzzyaaam32,xyzzyaaan32
logical xyzzyaaao32,xyzzyaaap32,xyzzyaaaq32
integer xyzzyaaar32(2),xyzzyaaas32(2,n,2),xyzzyaaat32,xyzzyaaau32
integer,parameter :: xyzzyaaav32=1024
integer xyzzyaaaw32(2),xyzzyaaax32(n,xyzzyaaav32,2),xyzzyaaay32
xyzzyaaaj32=2
xyzzyaaak32=1
xyzzyaaaw32(xyzzyaaak32)=ns
xyzzyaaax32(1:n,1:ns,xyzzyaaak32)=indx(1:n,1:ns)
xyzzyaaar32(xyzzyaaak32)=0
do xyzzyaaat32=1,namb_range
if(amb_range(2,xyzzyaaat32)>amb_range(1,xyzzyaaat32))then
xyzzyaaar32(xyzzyaaak32)=xyzzyaaar32(xyzzyaaak32)+1
xyzzyaaas32(1:2,xyzzyaaar32(xyzzyaaak32),xyzzyaaak32)=amb_range(1:2,xy&
&zzyaaat32)
endif
enddo
do xyzzyaaaa32=1,n
call swap1(xyzzyaaaj32,xyzzyaaak32)
do xyzzyaaau32=1,xyzzyaaar32(xyzzyaaaj32)
if(xyzzyaaas32(1,xyzzyaaau32,xyzzyaaaj32)==xyzzyaaaa32)exit
enddo
if(xyzzyaaau32>xyzzyaaar32(xyzzyaaaj32))xyzzyaaau32=0
xyzzyaaaf32=0
xyzzyaaao32=.false.
xyzzyaaaw32(xyzzyaaak32)=0
do xyzzyaaay32=1,xyzzyaaaw32(xyzzyaaaj32)
if(xyzzyaaau32==0)then
xyzzyaaal32=0
xyzzyaaam32=0
else
xyzzyaaal32=xyzzyaaas32(1,xyzzyaaau32,xyzzyaaaj32)
xyzzyaaam32=xyzzyaaas32(2,xyzzyaaau32,xyzzyaaaj32)
endif
do xyzzyaaac32=xyzzyaaal32,xyzzyaaam32
xyzzyaaaw32(xyzzyaaak32)=xyzzyaaaw32(xyzzyaaak32)+1
if(xyzzyaaaw32(xyzzyaaak32)>xyzzyaaav32)call errstop('SORT_MATRIX_SYMM&
&','Size of buffer too small.')
xyzzyaaax32(1:n,xyzzyaaaw32(xyzzyaaak32),xyzzyaaak32)=xyzzyaaax32(1:n,&
&xyzzyaaay32,xyzzyaaaj32)
if(xyzzyaaau32>0)then
call swap1(xyzzyaaax32(xyzzyaaaa32,xyzzyaaaw32(xyzzyaaak32),xyzzyaaak3&
&2),xyzzyaaax32(xyzzyaaac32,xyzzyaaaw32(xyzzyaaak32),xyzzyaaak32))
endif
xyzzyaaan32=xyzzyaaax32(xyzzyaaaa32,xyzzyaaaw32(xyzzyaaak32),xyzzyaaak&
&32)
do xyzzyaaat32=1,xyzzyaaar32(xyzzyaaaj32)
xyzzyaaah32=xyzzyaaas32(1,xyzzyaaat32,xyzzyaaaj32)
xyzzyaaai32=xyzzyaaas32(2,xyzzyaaat32,xyzzyaaaj32)
if(xyzzyaaat32==xyzzyaaau32)xyzzyaaah32=xyzzyaaah32+1
if(xyzzyaaai32<=xyzzyaaah32)cycle
call xyzzyaaaa1(n,xyzzyaaah32,xyzzyaaai32,matrix(1,xyzzyaaan32),xyzzya&
&aax32(1,xyzzyaaaw32(xyzzyaaak32),xyzzyaaak32))
enddo
xyzzyaaaq32=.false.
if(xyzzyaaao32)then
do xyzzyaaab32=xyzzyaaaa32,n
xyzzyaaad32=matrix(xyzzyaaax32(xyzzyaaab32,xyzzyaaaw32(xyzzyaaak32),xy&
&zzyaaak32),xyzzyaaan32)
if(xyzzyaaad32==xyzzyaaaf32(xyzzyaaab32))cycle
if(xyzzyaaad32<xyzzyaaaf32(xyzzyaaab32))then
xyzzyaaaq32=.true.
xyzzyaaax32(:,1,xyzzyaaak32)=xyzzyaaax32(:,xyzzyaaaw32(xyzzyaaak32),xy&
&zzyaaak32)
xyzzyaaaw32(xyzzyaaak32)=1
else
xyzzyaaaw32(xyzzyaaak32)=xyzzyaaaw32(xyzzyaaak32)-1
endif
exit
enddo
else
xyzzyaaaq32=.true.
endif
if(xyzzyaaaq32)then
do xyzzyaaab32=1,n
xyzzyaaaf32(xyzzyaaab32)=matrix(xyzzyaaax32(xyzzyaaab32,xyzzyaaaw32(xy&
&zzyaaak32),xyzzyaaak32),xyzzyaaan32)
enddo
xyzzyaaao32=.true.
endif
enddo
enddo
xyzzyaaar32(xyzzyaaak32)=0
do xyzzyaaat32=1,xyzzyaaar32(xyzzyaaaj32)
xyzzyaaae32=0
xyzzyaaap32=.false.
xyzzyaaag32=0
do xyzzyaaab32=max(xyzzyaaas32(1,xyzzyaaat32,xyzzyaaaj32),xyzzyaaaa32+&
&1),xyzzyaaas32(2,xyzzyaaat32,xyzzyaaaj32)
if(.not.xyzzyaaap32.or.xyzzyaaae32/=xyzzyaaaf32(xyzzyaaab32))then
if(xyzzyaaag32>1)then
xyzzyaaar32(xyzzyaaak32)=xyzzyaaar32(xyzzyaaak32)+1
xyzzyaaas32(1:2,xyzzyaaar32(xyzzyaaak32),xyzzyaaak32)=(/xyzzyaaab32-xy&
&zzyaaag32,xyzzyaaab32-1/)
endif
xyzzyaaag32=1
else
xyzzyaaag32=xyzzyaaag32+1
endif
xyzzyaaae32=xyzzyaaaf32(xyzzyaaab32)
xyzzyaaap32=.true.
enddo
if(xyzzyaaag32>1)then
xyzzyaaar32(xyzzyaaak32)=xyzzyaaar32(xyzzyaaak32)+1
xyzzyaaas32(1:2,xyzzyaaar32(xyzzyaaak32),xyzzyaaak32)=(/xyzzyaaas32(2,&
&xyzzyaaat32,xyzzyaaaj32)-xyzzyaaag32+1,xyzzyaaas32(2,xyzzyaaat32,xyzz&
&yaaaj32)/)
endif
enddo
if(xyzzyaaar32(xyzzyaaak32)==0.and.xyzzyaaaw32(xyzzyaaak32)==1)exit
enddo
ns=xyzzyaaaw32(xyzzyaaak32)
indx(1:n,1:ns)=xyzzyaaax32(1:n,1:ns,xyzzyaaak32)
namb_range=xyzzyaaar32(xyzzyaaak32)
if(namb_range>0)amb_range(1:2,1:namb_range)=xyzzyaaas32(1:2,1:namb_ran&
&ge,xyzzyaaak32)
end subroutine sort_matrix_symm
subroutine xyzzyaaaa1(n,i1,i2,x,indx)
implicit none
integer,intent(in) :: n,i1,i2,x(n)
integer,intent(inout) :: indx(n)
integer,parameter :: xyzzyaaaa33=30,xyzzyaaab33=128
integer xyzzyaaac33(xyzzyaaab33),xyzzyaaad33(xyzzyaaab33),xyzzyaaae33
integer xyzzyaaaf33,xyzzyaaag33,xyzzyaaah33,xyzzyaaai33,xyzzyaaaj33,xy&
&zzyaaak33,xyzzyaaal33,xyzzyaaam33
if(i2-i1>=xyzzyaaaa33)then
xyzzyaaae33=0
xyzzyaaaf33=i1
xyzzyaaag33=i2
do
xyzzyaaah33=xyzzyaaaf33-1
xyzzyaaai33=xyzzyaaag33
xyzzyaaam33=x(indx(xyzzyaaag33))
xyzzyaaak33=xyzzyaaaf33-1
xyzzyaaal33=xyzzyaaag33
do
xyzzyaaah33=xyzzyaaah33+1
do while(x(indx(xyzzyaaah33))<xyzzyaaam33)
xyzzyaaah33=xyzzyaaah33+1
enddo
xyzzyaaai33=xyzzyaaai33-1
do while(x(indx(xyzzyaaai33))>xyzzyaaam33)
if(xyzzyaaai33==xyzzyaaaf33)exit
xyzzyaaai33=xyzzyaaai33-1
enddo
if(xyzzyaaah33>=xyzzyaaai33)exit
call swap1(indx(xyzzyaaah33),indx(xyzzyaaai33))
if(x(indx(xyzzyaaah33))==xyzzyaaam33)then
xyzzyaaak33=xyzzyaaak33+1
call swap1(indx(xyzzyaaak33),indx(xyzzyaaah33))
endif
if(xyzzyaaam33==x(indx(xyzzyaaai33)))then
xyzzyaaal33=xyzzyaaal33-1
call swap1(indx(xyzzyaaai33),indx(xyzzyaaal33))
endif
enddo
call swap1(indx(xyzzyaaah33),indx(xyzzyaaag33))
xyzzyaaai33=xyzzyaaah33+1
xyzzyaaah33=xyzzyaaah33-1
do xyzzyaaaj33=xyzzyaaaf33,xyzzyaaak33-1
call swap1(indx(xyzzyaaaj33),indx(xyzzyaaah33))
xyzzyaaah33=xyzzyaaah33-1
enddo
do xyzzyaaaj33=xyzzyaaag33-1,xyzzyaaal33+1,-1
call swap1(indx(xyzzyaaai33),indx(xyzzyaaaj33))
xyzzyaaai33=xyzzyaaai33+1
enddo
if(xyzzyaaag33-xyzzyaaai33>=xyzzyaaah33-xyzzyaaaf33.and.xyzzyaaah33-xy&
&zzyaaaf33>=xyzzyaaaa33)then
xyzzyaaae33=xyzzyaaae33+1
if(xyzzyaaae33>xyzzyaaab33)call errstop('QSORT3_INT_PARTIAL_PREINIT','&
&Stack size too small <1>.')
xyzzyaaac33(xyzzyaaae33)=xyzzyaaai33
xyzzyaaad33(xyzzyaaae33)=xyzzyaaag33
xyzzyaaag33=xyzzyaaah33
elseif(xyzzyaaah33-xyzzyaaaf33>xyzzyaaag33-xyzzyaaai33.and.xyzzyaaag33&
&-xyzzyaaai33>=xyzzyaaaa33)then
xyzzyaaae33=xyzzyaaae33+1
if(xyzzyaaae33>xyzzyaaab33)call errstop('QSORT3_INT_PARTIAL_PREINIT','&
&Stack size too small <2>.')
xyzzyaaac33(xyzzyaaae33)=xyzzyaaaf33
xyzzyaaad33(xyzzyaaae33)=xyzzyaaah33
xyzzyaaaf33=xyzzyaaai33
elseif(xyzzyaaah33-xyzzyaaaf33>=xyzzyaaaa33)then
xyzzyaaag33=xyzzyaaah33
elseif(xyzzyaaag33-xyzzyaaai33>=xyzzyaaaa33)then
xyzzyaaaf33=xyzzyaaai33
else
if(xyzzyaaae33<1)exit
xyzzyaaaf33=xyzzyaaac33(xyzzyaaae33)
xyzzyaaag33=xyzzyaaad33(xyzzyaaae33)
xyzzyaaae33=xyzzyaaae33-1
endif
enddo
endif
do xyzzyaaah33=i1+1,i2
if(x(indx(xyzzyaaah33-1))<=x(indx(xyzzyaaah33)))cycle
call swap1(indx(xyzzyaaah33-1),indx(xyzzyaaah33))
do xyzzyaaai33=xyzzyaaah33-1,i1+1,-1
if(x(indx(xyzzyaaai33-1))<=x(indx(xyzzyaaai33)))exit
call swap1(indx(xyzzyaaai33-1),indx(xyzzyaaai33))
enddo
enddo
end subroutine xyzzyaaaa1
subroutine amb2cand_sort_matrix_rect(n,m,ns,row_indx,col_indx,namb_ran&
&ge,amb_range)
implicit none
integer,intent(in) :: n,m
integer,intent(inout) :: ns,row_indx(n,*),col_indx(m,*),namb_range,amb&
&_range(2,n)
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34,xyzzyaaae34,xy&
&zzyaaaf34,xyzzyaaag34(n),xyzzyaaah34
logical xyzzyaaai34,xyzzyaaaj34
if(namb_range<1)return
do xyzzyaaaa34=1,namb_range
xyzzyaaad34=amb_range(1,xyzzyaaaa34)
xyzzyaaae34=amb_range(2,xyzzyaaaa34)
xyzzyaaaf34=xyzzyaaae34-xyzzyaaad34+1
if(xyzzyaaaf34<=1)cycle
xyzzyaaac34=ns
do xyzzyaaab34=1,xyzzyaaac34
xyzzyaaai34=.true.
xyzzyaaaj34=xyzzyaaab1(xyzzyaaaf34,xyzzyaaag34,xyzzyaaai34)
do while(xyzzyaaab1(xyzzyaaaf34,xyzzyaaag34,xyzzyaaai34))
ns=ns+1
col_indx(1:m,ns)=col_indx(1:m,xyzzyaaab34)
if(xyzzyaaad34>1)row_indx(1:xyzzyaaad34-1,ns)=row_indx(1:xyzzyaaad34-1&
&,xyzzyaaab34)
do xyzzyaaah34=1,xyzzyaaaf34
row_indx(xyzzyaaah34+xyzzyaaad34-1,ns)=row_indx(xyzzyaaag34(xyzzyaaah3&
&4)+xyzzyaaad34-1,xyzzyaaab34)
enddo
if(xyzzyaaae34<n)row_indx(xyzzyaaae34+1:n,ns)=row_indx(xyzzyaaae34+1:n&
&,xyzzyaaab34)
enddo
enddo
enddo
namb_range=0
end subroutine amb2cand_sort_matrix_rect
subroutine amb2cand_sort_matrix_symm(n,ns,indx,namb_range,amb_range)
implicit none
integer,intent(in) :: n
integer,intent(inout) :: ns,indx(n,*),namb_range,amb_range(2,n)
integer xyzzyaaaa35,xyzzyaaab35,xyzzyaaac35,xyzzyaaad35,xyzzyaaae35,xy&
&zzyaaaf35,xyzzyaaag35(n),xyzzyaaah35
logical xyzzyaaai35,xyzzyaaaj35
if(namb_range<1)return
do xyzzyaaaa35=1,namb_range
xyzzyaaad35=amb_range(1,xyzzyaaaa35)
xyzzyaaae35=amb_range(2,xyzzyaaaa35)
xyzzyaaaf35=xyzzyaaae35-xyzzyaaad35+1
if(xyzzyaaaf35<=1)cycle
xyzzyaaac35=ns
do xyzzyaaab35=1,xyzzyaaac35
xyzzyaaai35=.true.
xyzzyaaaj35=xyzzyaaab1(xyzzyaaaf35,xyzzyaaag35,xyzzyaaai35)
do while(xyzzyaaab1(xyzzyaaaf35,xyzzyaaag35,xyzzyaaai35))
ns=ns+1
if(xyzzyaaad35>1)indx(1:xyzzyaaad35-1,ns)=indx(1:xyzzyaaad35-1,xyzzyaa&
&ab35)
do xyzzyaaah35=1,xyzzyaaaf35
indx(xyzzyaaah35+xyzzyaaad35-1,ns)=indx(xyzzyaaag35(xyzzyaaah35)+xyzzy&
&aaad35-1,xyzzyaaab35)
enddo
if(xyzzyaaae35<n)indx(xyzzyaaae35+1:n,ns)=indx(xyzzyaaae35+1:n,xyzzyaa&
&ab35)
enddo
enddo
enddo
namb_range=0
end subroutine amb2cand_sort_matrix_symm
logical function xyzzyaaab1(n,indices,is_first)
integer,intent(in) :: n
integer,intent(inout) :: indices(n)
logical,intent(inout) :: is_first
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36,xyzzyaaad36
logical xyzzyaaae36(n,n)
if(is_first)then
do xyzzyaaaa36=1,n
indices(xyzzyaaaa36)=xyzzyaaaa36
enddo
xyzzyaaab1=.true.
is_first=.false.
return
endif
xyzzyaaae36=.true.
do xyzzyaaaa36=1,n
xyzzyaaae36(indices(xyzzyaaaa36),xyzzyaaaa36:n)=.false.
enddo
do xyzzyaaaa36=n,1,-1
xyzzyaaae36(indices(xyzzyaaaa36),xyzzyaaaa36:n)=.true.
do xyzzyaaab36=indices(xyzzyaaaa36)+1,n
if(xyzzyaaae36(xyzzyaaab36,xyzzyaaaa36))then
indices(xyzzyaaaa36)=xyzzyaaab36
xyzzyaaae36(xyzzyaaab36,xyzzyaaaa36:n)=.false.
do xyzzyaaac36=xyzzyaaaa36+1,n
do xyzzyaaad36=1,n
if(xyzzyaaae36(xyzzyaaad36,xyzzyaaac36))then
indices(xyzzyaaac36)=xyzzyaaad36
xyzzyaaae36(xyzzyaaad36,xyzzyaaac36:n)=.false.
exit
endif
enddo
enddo
xyzzyaaab1=.true.
return
endif
enddo
enddo
xyzzyaaab1=.false.
end function xyzzyaaab1
subroutine gmult(v,n,m,ao,phase,norb)
implicit none
integer,intent(in) :: n,m,norb
real(dp),intent(in) :: phase(m),ao(norb)
real(dp),intent(inout) :: v(n,*)
integer xyzzyaaaa37
do xyzzyaaaa37=1,m
call daxpy(norb,phase(xyzzyaaaa37),ao(1),1,v(1,xyzzyaaaa37),1)
enddo
end subroutine gmult
subroutine g3mult(v1,v2,v3,n,m,agra,phase,norb)
implicit none
integer,intent(in) :: n,m,norb
real(dp),intent(in) :: phase(m),agra(3,norb)
real(dp),intent(inout) :: v1(n,*),v2(n,*),v3(n,*)
integer xyzzyaaaa38
do xyzzyaaaa38=1,m
call daxpy(norb,phase(xyzzyaaaa38),agra(1,1),3,v1(1,xyzzyaaaa38),1)
call daxpy(norb,phase(xyzzyaaaa38),agra(2,1),3,v2(1,xyzzyaaaa38),1)
call daxpy(norb,phase(xyzzyaaaa38),agra(3,1),3,v3(1,xyzzyaaaa38),1)
enddo
end subroutine g3mult
subroutine g6mult(v1,v2,v3,v4,v5,v6,n,m,asderiv,phase,norb)
implicit none
integer,intent(in) :: n,m,norb
real(dp),intent(in) :: phase(m),asderiv(6,norb)
real(dp),intent(inout) :: v1(n,*),v2(n,*),v3(n,*),v4(n,*),v5(n,*),v6(n&
&,*)
integer xyzzyaaaa39
do xyzzyaaaa39=1,m
call daxpy(norb,phase(xyzzyaaaa39),asderiv(1,1),6,v1(1,xyzzyaaaa39),1)
call daxpy(norb,phase(xyzzyaaaa39),asderiv(2,1),6,v2(1,xyzzyaaaa39),1)
call daxpy(norb,phase(xyzzyaaaa39),asderiv(3,1),6,v3(1,xyzzyaaaa39),1)
call daxpy(norb,phase(xyzzyaaaa39),asderiv(4,1),6,v4(1,xyzzyaaaa39),1)
call daxpy(norb,phase(xyzzyaaaa39),asderiv(5,1),6,v5(1,xyzzyaaaa39),1)
call daxpy(norb,phase(xyzzyaaaa39),asderiv(6,1),6,v6(1,xyzzyaaaa39),1)
enddo
end subroutine g6mult
subroutine mxmb(a,size_a,nca,nra,b,ncb,nrb,r,ncr,nrr,ncol,nlink,nrow)
implicit none
integer,intent(in) :: nca,nra,ncb,nrb,ncr,nrr,ncol,nlink,nrow,size_a
real(dp),intent(in) :: a(size_a),b(*)
real(dp),intent(inout) :: r(*)
integer xyzzyaaaa40,xyzzyaaab40,xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xy&
&zzyaaaf40,xyzzyaaag40,xyzzyaaah40,xyzzyaaai40,xyzzyaaaj40,xyzzyaaak40&
&,xyzzyaaal40,xyzzyaaam40,xyzzyaaan40,xyzzyaaao40,xyzzyaaap40,xyzzyaaa&
&q40,xyzzyaaar40,xyzzyaaas40,xyzzyaaat40,xyzzyaaau40,xyzzyaaav40,xyzzy&
&aaaw40,xyzzyaaax40,xyzzyaaay40
real(dp) xyzzyaaaz40,xyzzyaaba40,xyzzyaabb40,xyzzyaabc40,xyzzyaabd40,x&
&yzzyaabe40,xyzzyaabf40
xyzzyaaan40=ncol/2
xyzzyaaam40=nlink*nra+1
xyzzyaaao40=nca+nca
xyzzyaaap40=ncr+ncr
xyzzyaaad40=1
xyzzyaaac40=1
if((xyzzyaaan40+xyzzyaaan40)==ncol)then
do xyzzyaaaa40=1,nrow
xyzzyaaak40=xyzzyaaac40
xyzzyaaab40=1
1   if(xyzzyaaab40==xyzzyaaam40)goto 15
xyzzyaaaz40=b(xyzzyaaak40)
xyzzyaaaf40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaaaz40==0.d0)goto 1
xyzzyaaaw40=xyzzyaaad40+ncr
xyzzyaaaq40=xyzzyaaaf40+nca
3   if(xyzzyaaab40==xyzzyaaam40)goto 14
xyzzyaaba40=b(xyzzyaaak40)
xyzzyaaag40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaaba40==0.d0)goto 3
xyzzyaaar40=xyzzyaaag40+nca
5   if(xyzzyaaab40==xyzzyaaam40)goto 13
xyzzyaabb40=b(xyzzyaaak40)
xyzzyaaah40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaabb40==0.d0)goto 5
xyzzyaaas40=xyzzyaaah40+nca
7   if(xyzzyaaab40==xyzzyaaam40)goto 12
xyzzyaabc40=b(xyzzyaaak40)
xyzzyaaai40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaabc40==0.d0)goto 7
xyzzyaaat40=xyzzyaaai40+nca
9   if(xyzzyaaab40==xyzzyaaam40)goto 11
xyzzyaabd40=b(xyzzyaaak40)
xyzzyaaaj40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaabd40==0.d0)goto 9
xyzzyaaau40=xyzzyaaaj40+nca
xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&h40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaai40+xyzzyaaae40)+xyzzyaabd40*a&
&(xyzzyaaaj40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&s40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaat40+xyzzyaaae40)+xyzzyaabd40*a&
&(xyzzyaaau40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
goto 1
11  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&h40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaai40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&s40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaat40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
goto 15
12  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&h40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&s40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
goto 15
13  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
goto 15
14  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
15  xyzzyaaad40=xyzzyaaad40+nrr
xyzzyaaac40=xyzzyaaac40+nrb
enddo
else
xyzzyaaax40=xyzzyaaan40*xyzzyaaap40
xyzzyaaay40=xyzzyaaan40*xyzzyaaao40
do xyzzyaaaa40=1,nrow
xyzzyaaak40=xyzzyaaac40
xyzzyaaab40=1
16  if(xyzzyaaab40==xyzzyaaam40)goto 29
xyzzyaaaz40=b(xyzzyaaak40)
xyzzyaaaf40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaaaz40==0.d0)goto 16
17  xyzzyaaaw40=xyzzyaaad40+ncr
xyzzyaaaq40=xyzzyaaaf40+nca
if(xyzzyaaab40==xyzzyaaam40)goto 28
xyzzyaaba40=b(xyzzyaaak40)
xyzzyaaag40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaaba40==0.d0)goto 17
xyzzyaaar40=xyzzyaaag40+nca
19  if(xyzzyaaab40==xyzzyaaam40)goto 27
xyzzyaabb40=b(xyzzyaaak40)
xyzzyaaah40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaabb40==0.d0)goto 19
xyzzyaaas40=xyzzyaaah40+nca
21  if(xyzzyaaab40==xyzzyaaam40)goto 26
xyzzyaabc40=b(xyzzyaaak40)
xyzzyaaai40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaabc40==0.d0)goto 21
xyzzyaaat40=xyzzyaaai40+nca
23  if(xyzzyaaab40==xyzzyaaam40)goto 25
xyzzyaabd40=b(xyzzyaaak40)
xyzzyaaaj40=xyzzyaaab40
xyzzyaaak40=xyzzyaaak40+ncb
xyzzyaaab40=xyzzyaaab40+nra
if(xyzzyaabd40==0.d0)goto 23
xyzzyaaau40=xyzzyaaaj40+nca
xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&h40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaai40+xyzzyaaae40)+xyzzyaabd40*a&
&(xyzzyaaaj40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&s40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaat40+xyzzyaaae40)+xyzzyaabd40*a&
&(xyzzyaaau40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
r(xyzzyaaad40+xyzzyaaax40)=r(xyzzyaaad40+xyzzyaaax40)+xyzzyaaaz40*a(xy&
&zzyaaaf40+xyzzyaaay40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaay40)+xyzzyaa&
&bb40*a(xyzzyaaah40+xyzzyaaay40)+xyzzyaabc40*a(xyzzyaaai40+xyzzyaaay40&
&)+xyzzyaabd40*a(xyzzyaaaj40+xyzzyaaay40)
goto 16
25  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&h40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaai40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&s40+xyzzyaaae40)+xyzzyaabc40*a(xyzzyaaat40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
r(xyzzyaaad40+xyzzyaaax40)=r(xyzzyaaad40+xyzzyaaax40)+xyzzyaaaz40*a(xy&
&zzyaaaf40+xyzzyaaay40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaay40)+xyzzyaa&
&bb40*a(xyzzyaaah40+xyzzyaaay40)+xyzzyaabc40*a(xyzzyaaai40+xyzzyaaay40&
&)
goto 29
26  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&h40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)+xyzzyaabb40*a(xyzzyaaa&
&s40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
r(xyzzyaaad40+xyzzyaaax40)=r(xyzzyaaad40+xyzzyaaax40)+xyzzyaaaz40*a(xy&
&zzyaaaf40+xyzzyaaay40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaay40)+xyzzyaa&
&bb40*a(xyzzyaaah40+xyzzyaaay40)
goto 29
27  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)+xyzzyaaba40*a(xyzzyaaar40+xyzzyaaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
r(xyzzyaaad40+xyzzyaaax40)=r(xyzzyaaad40+xyzzyaaax40)+xyzzyaaaz40*a(xy&
&zzyaaaf40+xyzzyaaay40)+xyzzyaaba40*a(xyzzyaaag40+xyzzyaaay40)
goto 29
28  xyzzyaaal40=0
xyzzyaaae40=0
do xyzzyaaav40=1,xyzzyaaan40
xyzzyaabe40=r(xyzzyaaad40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaf40+xyzzy&
&aaae40)
xyzzyaabf40=r(xyzzyaaaw40+xyzzyaaal40)+xyzzyaaaz40*a(xyzzyaaaq40+xyzzy&
&aaae40)
r(xyzzyaaad40+xyzzyaaal40)=xyzzyaabe40
r(xyzzyaaaw40+xyzzyaaal40)=xyzzyaabf40
xyzzyaaal40=xyzzyaaal40+xyzzyaaap40
xyzzyaaae40=xyzzyaaae40+xyzzyaaao40
enddo
r(xyzzyaaad40+xyzzyaaax40)=r(xyzzyaaad40+xyzzyaaax40)+xyzzyaaaz40*a(xy&
&zzyaaaf40+xyzzyaaay40)
29  xyzzyaaad40=xyzzyaaad40+nrr
xyzzyaaac40=xyzzyaaac40+nrb
enddo
endif
end subroutine mxmb
subroutine mxmbc(a,size_a,nca,nra,b,ncb,nrb,r,ncr,nrr,ncol,nlink,nrow)
implicit none
integer,intent(in) :: size_a,nca,nra,ncb,nrb,ncr,nrr,ncol,nlink,nrow
complex(dp),intent(in) :: a(size_a),b(*)
complex(dp),intent(inout) :: r(*)
integer xyzzyaaaa41,xyzzyaaab41,xyzzyaaac41,xyzzyaaad41,xyzzyaaae41,xy&
&zzyaaaf41,xyzzyaaag41,xyzzyaaah41,xyzzyaaai41,xyzzyaaaj41,xyzzyaaak41&
&,xyzzyaaal41,xyzzyaaam41,xyzzyaaan41,xyzzyaaao41,xyzzyaaap41,xyzzyaaa&
&q41,xyzzyaaar41,xyzzyaaas41,xyzzyaaat41,xyzzyaaau41,xyzzyaaav41,xyzzy&
&aaaw41,xyzzyaaax41,xyzzyaaay41
complex(dp) xyzzyaaaz41,xyzzyaaba41,xyzzyaabb41,xyzzyaabc41,xyzzyaabd4&
&1,xyzzyaabe41,xyzzyaabf41
xyzzyaaan41=ncol/2
xyzzyaaam41=nlink*nra+1
xyzzyaaao41=nca+nca
xyzzyaaap41=ncr+ncr
xyzzyaaad41=1
xyzzyaaac41=1
if((xyzzyaaan41+xyzzyaaan41)==ncol)then
do xyzzyaaaa41=1,nrow
xyzzyaaak41=xyzzyaaac41
xyzzyaaab41=1
1   if(xyzzyaaab41==xyzzyaaam41)goto 15
xyzzyaaaz41=b(xyzzyaaak41)
xyzzyaaaf41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaaaz41==czero)goto 1
xyzzyaaaw41=xyzzyaaad41+ncr
xyzzyaaaq41=xyzzyaaaf41+nca
3   if(xyzzyaaab41==xyzzyaaam41)goto 14
xyzzyaaba41=b(xyzzyaaak41)
xyzzyaaag41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaaba41==czero)goto 3
xyzzyaaar41=xyzzyaaag41+nca
5   if(xyzzyaaab41==xyzzyaaam41)goto 13
xyzzyaabb41=b(xyzzyaaak41)
xyzzyaaah41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaabb41==czero)goto 5
xyzzyaaas41=xyzzyaaah41+nca
7   if(xyzzyaaab41==xyzzyaaam41)goto 12
xyzzyaabc41=b(xyzzyaaak41)
xyzzyaaai41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaabc41==czero)goto 7
xyzzyaaat41=xyzzyaaai41+nca
9   if(xyzzyaaab41==xyzzyaaam41)goto 11
xyzzyaabd41=b(xyzzyaaak41)
xyzzyaaaj41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaabd41==czero)goto 9
xyzzyaaau41=xyzzyaaaj41+nca
xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&h41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaai41+xyzzyaaae41)+xyzzyaabd41*a&
&(xyzzyaaaj41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&s41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaat41+xyzzyaaae41)+xyzzyaabd41*a&
&(xyzzyaaau41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
goto 1
11  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&h41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaai41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&s41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaat41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
goto 15
12  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&h41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&s41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
goto 15
13  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
goto 15
14  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
15  xyzzyaaad41=xyzzyaaad41+nrr
xyzzyaaac41=xyzzyaaac41+nrb
enddo
else
xyzzyaaax41=xyzzyaaan41*xyzzyaaap41
xyzzyaaay41=xyzzyaaan41*xyzzyaaao41
do xyzzyaaaa41=1,nrow
xyzzyaaak41=xyzzyaaac41
xyzzyaaab41=1
16  if(xyzzyaaab41==xyzzyaaam41)goto 29
xyzzyaaaz41=b(xyzzyaaak41)
xyzzyaaaf41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaaaz41==czero)goto 16
17  xyzzyaaaw41=xyzzyaaad41+ncr
xyzzyaaaq41=xyzzyaaaf41+nca
if(xyzzyaaab41==xyzzyaaam41)goto 28
xyzzyaaba41=b(xyzzyaaak41)
xyzzyaaag41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaaba41==czero)goto 17
xyzzyaaar41=xyzzyaaag41+nca
19  if(xyzzyaaab41==xyzzyaaam41)goto 27
xyzzyaabb41=b(xyzzyaaak41)
xyzzyaaah41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaabb41==czero)goto 19
xyzzyaaas41=xyzzyaaah41+nca
21  if(xyzzyaaab41==xyzzyaaam41)goto 26
xyzzyaabc41=b(xyzzyaaak41)
xyzzyaaai41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaabc41==czero)goto 21
xyzzyaaat41=xyzzyaaai41+nca
23  if(xyzzyaaab41==xyzzyaaam41)goto 25
xyzzyaabd41=b(xyzzyaaak41)
xyzzyaaaj41=xyzzyaaab41
xyzzyaaak41=xyzzyaaak41+ncb
xyzzyaaab41=xyzzyaaab41+nra
if(xyzzyaabd41==czero)goto 23
xyzzyaaau41=xyzzyaaaj41+nca
xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&h41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaai41+xyzzyaaae41)+xyzzyaabd41*a&
&(xyzzyaaaj41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&s41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaat41+xyzzyaaae41)+xyzzyaabd41*a&
&(xyzzyaaau41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
r(xyzzyaaad41+xyzzyaaax41)=r(xyzzyaaad41+xyzzyaaax41)+xyzzyaaaz41*a(xy&
&zzyaaaf41+xyzzyaaay41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaay41)+xyzzyaa&
&bb41*a(xyzzyaaah41+xyzzyaaay41)+xyzzyaabc41*a(xyzzyaaai41+xyzzyaaay41&
&)+xyzzyaabd41*a(xyzzyaaaj41+xyzzyaaay41)
goto 16
25  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&h41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaai41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&s41+xyzzyaaae41)+xyzzyaabc41*a(xyzzyaaat41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
r(xyzzyaaad41+xyzzyaaax41)=r(xyzzyaaad41+xyzzyaaax41)+xyzzyaaaz41*a(xy&
&zzyaaaf41+xyzzyaaay41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaay41)+xyzzyaa&
&bb41*a(xyzzyaaah41+xyzzyaaay41)+xyzzyaabc41*a(xyzzyaaai41+xyzzyaaay41&
&)
goto 29
26  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&h41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)+xyzzyaabb41*a(xyzzyaaa&
&s41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
r(xyzzyaaad41+xyzzyaaax41)=r(xyzzyaaad41+xyzzyaaax41)+xyzzyaaaz41*a(xy&
&zzyaaaf41+xyzzyaaay41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaay41)+xyzzyaa&
&bb41*a(xyzzyaaah41+xyzzyaaay41)
goto 29
27  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)+xyzzyaaba41*a(xyzzyaaar41+xyzzyaaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
r(xyzzyaaad41+xyzzyaaax41)=r(xyzzyaaad41+xyzzyaaax41)+xyzzyaaaz41*a(xy&
&zzyaaaf41+xyzzyaaay41)+xyzzyaaba41*a(xyzzyaaag41+xyzzyaaay41)
goto 29
28  xyzzyaaal41=0
xyzzyaaae41=0
do xyzzyaaav41=1,xyzzyaaan41
xyzzyaabe41=r(xyzzyaaad41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaf41+xyzzy&
&aaae41)
xyzzyaabf41=r(xyzzyaaaw41+xyzzyaaal41)+xyzzyaaaz41*a(xyzzyaaaq41+xyzzy&
&aaae41)
r(xyzzyaaad41+xyzzyaaal41)=xyzzyaabe41
r(xyzzyaaaw41+xyzzyaaal41)=xyzzyaabf41
xyzzyaaal41=xyzzyaaal41+xyzzyaaap41
xyzzyaaae41=xyzzyaaae41+xyzzyaaao41
enddo
r(xyzzyaaad41+xyzzyaaax41)=r(xyzzyaaad41+xyzzyaaax41)+xyzzyaaaz41*a(xy&
&zzyaaaf41+xyzzyaaay41)
29  xyzzyaaad41=xyzzyaaad41+nrr
xyzzyaaac41=xyzzyaaac41+nrb
enddo
endif
end subroutine mxmbc
real(dp) function exp_int(k)
use slaarnaag,only : euler
implicit none
real(dp),intent(in) :: k
integer xyzzyaaaa42
integer,parameter :: xyzzyaaab42=500
real(dp) xyzzyaaac42,xyzzyaaad42,xyzzyaaae42,xyzzyaaaf42,xyzzyaaag42,x&
&yzzyaaah42,xyzzyaaai42
real(dp),parameter :: xyzzyaaaj42=1.d-16,xyzzyaaak42=1.d-30
if(k<=0)then
call errstop("EXP_INT","Cannot evaluate exponential integral at k<=0")
elseif(k>1)then
xyzzyaaad42=k+1.d0
xyzzyaaae42=1.d0/xyzzyaaak42
xyzzyaaaf42=1.d0/xyzzyaaad42
xyzzyaaai42=xyzzyaaaf42
do xyzzyaaaa42=1,xyzzyaaab42
xyzzyaaac42=dble(-xyzzyaaaa42*xyzzyaaaa42)
xyzzyaaad42=xyzzyaaad42+2.d0
xyzzyaaaf42=1.d0/(xyzzyaaac42*xyzzyaaaf42+xyzzyaaad42)
xyzzyaaae42=xyzzyaaad42+xyzzyaaac42/xyzzyaaae42
xyzzyaaag42=xyzzyaaae42*xyzzyaaaf42
xyzzyaaai42=xyzzyaaai42*xyzzyaaag42
if(abs(xyzzyaaag42-1.d0)<xyzzyaaaj42)then
exp_int=xyzzyaaai42*exp(-k)
return
endif
enddo
call errstop("EXP_INT","Evaluation of the exponential integral failed &
&for k>1")
else
exp_int=-log(k)-euler
xyzzyaaah42=1.d0
do xyzzyaaaa42=1,xyzzyaaab42
xyzzyaaah42=-xyzzyaaah42*k/dble(xyzzyaaaa42)
xyzzyaaag42=-xyzzyaaah42/dble(xyzzyaaaa42)
exp_int=exp_int+xyzzyaaag42
if(abs(xyzzyaaag42)<abs(exp_int)*xyzzyaaaj42)return
enddo
call errstop("EXP_INT","Evaluation of the exponential integral failed &
&for k<1")
endif
end function exp_int
real(dp) function erfc(x)
implicit none
real(dp),intent(in) :: x
real(dp) xyzzyaaaa43
xyzzyaaaa43=x*x
if(xyzzyaaaa43<1.5d0)then
erfc=1.d0-gamma_ser(xyzzyaaaa43)
else
erfc=xyzzyaaac1(xyzzyaaaa43)
endif
if(x<0.d0)erfc=2.d0-erfc
end function erfc
real(dp) function gamma_ser(x)
use slaarnaag,only : half_log_pi
implicit none
real(dp),intent(in) :: x
integer n
integer,parameter :: maxit=100
real(dp) xyzzyaaaa44,xyzzyaaab44,xyzzyaaac44
real(dp),parameter :: xyzzyaaad44=1.d-12
if(x<=0.d0)then
if(x<0.d0)call errstop('GAMMA_SER','Bad argument in series representat&
&ion of gamma function.')
gamma_ser=0.d0
return
endif
xyzzyaaaa44=0.5d0
xyzzyaaac44=2.d0
xyzzyaaab44=xyzzyaaac44
do n=1,maxit
xyzzyaaaa44=xyzzyaaaa44+1.d0
xyzzyaaab44=xyzzyaaab44*x/xyzzyaaaa44
xyzzyaaac44=xyzzyaaac44+xyzzyaaab44
if(abs(xyzzyaaab44)<abs(xyzzyaaac44)*xyzzyaaad44)exit
enddo
if(n==maxit+1)call errstop('GAMMA_SER','Series representation of gamma&
& function failed to converge.')
gamma_ser=xyzzyaaac44*exp_protect(-x+0.5d0*log(x)-half_log_pi)
end function gamma_ser
real(dp) function xyzzyaaac1(x)
use slaarnaag,only : half_log_pi
implicit none
real(dp),intent(in) :: x
integer i
integer,parameter :: maxit=100
real(dp) xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45,xyzzyaaad45,xyzzyaaae45,x&
&yzzyaaaf45
real(dp),parameter :: xyzzyaaag45=1.d-12,xyzzyaaah45=1.d-30
xyzzyaaab45=x+0.5d0
xyzzyaaac45=1.d0/xyzzyaaah45
xyzzyaaad45=1.d0/xyzzyaaab45
xyzzyaaaf45=xyzzyaaad45
do i=1,maxit
xyzzyaaaa45=dble(i)*(0.5d0-dble(i))
xyzzyaaab45=xyzzyaaab45+2.d0
xyzzyaaad45=xyzzyaaaa45*xyzzyaaad45+xyzzyaaab45
if(abs(xyzzyaaad45)<xyzzyaaah45)xyzzyaaad45=xyzzyaaah45
xyzzyaaac45=xyzzyaaab45+xyzzyaaaa45/xyzzyaaac45
if(abs(xyzzyaaac45)<xyzzyaaah45)xyzzyaaac45=xyzzyaaah45
xyzzyaaad45=1.d0/xyzzyaaad45
xyzzyaaae45=xyzzyaaad45*xyzzyaaac45
xyzzyaaaf45=xyzzyaaaf45*xyzzyaaae45
if(abs(xyzzyaaae45-1.d0)<xyzzyaaag45)exit
enddo
if(i==maxit+1)call errstop('GAMMA_CF','Continued fraction representati&
&on of gamma function failed to converge.')
xyzzyaaac1=exp_protect(-x+0.5d0*log(x)-half_log_pi)*xyzzyaaaf45
end function xyzzyaaac1
subroutine lu_decom(a,piv,np,n,issingular)
implicit none
integer,intent(in) :: n,np
integer,intent(inout) :: piv(n)
real(dp),intent(inout) :: a(np,np)
logical,intent(out) :: issingular
integer xyzzyaaaa46
call dgetrf(n,n,a,np,piv,xyzzyaaaa46)
if(xyzzyaaaa46<0)call errstop('LU_DECOM','DGETRF says parameter #'//tr&
&im(i2s(-xyzzyaaaa46))//' has an illegal value.')
issingular=xyzzyaaaa46>0
end subroutine lu_decom
subroutine lu_decom_cmplx(a,piv,n,issingular)
implicit none
integer,intent(in) :: n
complex(dp),intent(inout) :: a(n,n)
integer,intent(inout) :: piv(n)
logical,intent(out) :: issingular
integer xyzzyaaaa47
call zgetrf(n,n,a,n,piv,xyzzyaaaa47)
if(xyzzyaaaa47<0)call errstop('LU_DECOMP_CMPLX','ZGETRF says parameter&
& #'//trim(i2s(-xyzzyaaaa47))//' has an illegal value.')
issingular=xyzzyaaaa47>0
end subroutine lu_decom_cmplx
subroutine lu_decom_cmplx_dz(a,awork,piv,np,n,issingular)
implicit none
integer,intent(in) :: n,np
integer,intent(inout) :: piv(n)
real(dp),intent(inout) :: a(np,np,2)
complex(dp),intent(inout) :: awork(np,np)
logical,intent(out) :: issingular
integer xyzzyaaaa48
call xyzzyaaad1(np*np,a(1,1,1),a(1,1,2),1,awork(1,1),1)
call zgetrf(n,n,awork,np,piv,xyzzyaaaa48)
if(xyzzyaaaa48<0)call errstop('LU_DECOM_CMPLX_DZ','ZGETRF says paramet&
&er #'//trim(i2s(-xyzzyaaaa48))//' has an illegal value.')
call xyzzyaaae1(np*np,awork(1,1),1,a(1,1,1),a(1,1,2),1)
issingular=xyzzyaaaa48>0
end subroutine lu_decom_cmplx_dz
subroutine lu_solve_once(a,piv,x,np,n)
implicit none
integer,intent(in) :: np,n,piv(n)
real(dp),intent(in) :: a(np,n)
real(dp),intent(inout) :: x(n)
integer xyzzyaaaa49
call dgetrs('N',n,1,a,np,piv,x,n,xyzzyaaaa49)
if(xyzzyaaaa49<0)call errstop('LU_SOLVE_ONCE','DGETRS says parameter #&
&'//trim(i2s(-xyzzyaaaa49))//' has an illegal value.')
end subroutine lu_solve_once
subroutine lu_solve_n(a,piv,x,np,n,npx)
implicit none
integer,intent(in) :: np,n,npx,piv(n)
real(dp),intent(in) :: a(np,n)
real(dp),intent(inout) :: x(npx,n)
integer xyzzyaaaa50
call dgetrs('N',n,n,a,np,piv,x,npx,xyzzyaaaa50)
if(xyzzyaaaa50<0)call errstop('LU_SOLVE_N','DGETRS says parameter #'//&
&trim(i2s(-xyzzyaaaa50))//' has an illegal value.')
end subroutine lu_solve_n
subroutine lu_solve_once_cmplx(a,x,piv,np,n)
implicit none
integer,intent(in) :: np,n,piv(n)
complex(dp),intent(in) :: a(np,n)
complex(dp),intent(inout) :: x(n)
integer xyzzyaaaa51
call zgetrs('N',n,1,a,np,piv,x,n,xyzzyaaaa51)
if(xyzzyaaaa51<0)call errstop('LU_SOLVE_CMPLX_ONCE','ZGETRS says param&
&eter #'//trim(i2s(-xyzzyaaaa51))//' has an illegal value.')
end subroutine lu_solve_once_cmplx
subroutine lu_solve_n_cmplx(a,piv,x,np,n,npx)
implicit none
integer,intent(in) :: np,n,npx,piv(n)
complex(dp),intent(in) :: a(np,n)
complex(dp),intent(inout) :: x(npx,n)
integer xyzzyaaaa52
call zgetrs('N',n,n,a,np,piv,x,npx,xyzzyaaaa52)
if(xyzzyaaaa52<0)call errstop('LU_SOLVE_N_CMPLX','ZGETRS says paramete&
&r #'//trim(i2s(-xyzzyaaaa52))//' has an illegal value.')
end subroutine lu_solve_n_cmplx
subroutine lu_solve_once_cmplx_dz(a,awork,piv,x,xwork,np,n)
implicit none
integer,intent(in) :: np,n,piv(n)
real(dp),intent(in) :: a(np,np,2)
real(dp),intent(inout) :: x(np,2)
complex(dp),intent(inout) :: awork(np,n),xwork(n)
integer xyzzyaaaa53
call xyzzyaaad1(np*n,a(1,1,1),a(1,1,2),1,awork(1,1),1)
call xyzzyaaad1(n,x(1,1),x(1,2),1,xwork(1),1)
call zgetrs('N',n,1,awork,np,piv,xwork,n,xyzzyaaaa53)
if(xyzzyaaaa53<0)call errstop('LU_SOLVE_CMPLX_ONCE_DZ','ZGETRS says pa&
&rameter #'//trim(i2s(-xyzzyaaaa53))//' has an illegal value.')
call xyzzyaaae1(n,xwork,1,x(1,1),x(1,2),1)
end subroutine lu_solve_once_cmplx_dz
subroutine lu_solve_n_cmplx_dz(a,awork,piv,x,xwork,np,n,npx)
implicit none
integer,intent(in) :: np,n,npx,piv(n)
real(dp),intent(in) :: a(np,np,2)
complex(dp),intent(inout) :: awork(np,n),xwork(npx,n)
real(dp),intent(inout) :: x(npx,npx,2)
integer xyzzyaaaa54
call xyzzyaaad1(np*n,a(1,1,1),a(1,1,2),1,awork(1,1),1)
call xyzzyaaad1(npx*n,x(1,1,1),x(1,1,2),1,xwork(1,1),1)
call zgetrs('N',n,n,awork,np,piv,xwork,npx,xyzzyaaaa54)
if(xyzzyaaaa54<0)call errstop('LU_SOLVE_CMPLX_N_DZ','ZGETRS says param&
&eter #'//trim(i2s(-xyzzyaaaa54))//' has an illegal value.')
call xyzzyaaae1(npx*n,xwork,1,x(1,1,1),x(1,1,2),1)
end subroutine lu_solve_n_cmplx_dz
complex(dp) function lu_logdet(a,piv,np,n)
implicit none
integer,intent(in) :: n,np,piv(n)
real(dp),intent(in) :: a(np,n)
integer xyzzyaaaa55
real(dp) xyzzyaaab55,xyzzyaaac55
logical xyzzyaaad55
xyzzyaaad55=.false.
xyzzyaaac55=a(n,n)
xyzzyaaab55=log(abs(xyzzyaaac55))
if(xyzzyaaac55<0.d0)xyzzyaaad55=.true.
do xyzzyaaaa55=1,n-1
xyzzyaaac55=a(xyzzyaaaa55,xyzzyaaaa55)
xyzzyaaab55=xyzzyaaab55+log(abs(xyzzyaaac55))
if(xyzzyaaac55<0.d0.neqv.xyzzyaaaa55/=piv(xyzzyaaaa55))xyzzyaaad55=.no&
&t.xyzzyaaad55
enddo
if(xyzzyaaad55)then
lu_logdet=cmplx(xyzzyaaab55,pi,dp)
else
lu_logdet=cmplx(xyzzyaaab55,0.d0,dp)
endif
end function lu_logdet
complex(dp) function lu_logdet_cmplx(a,piv,np,n)
implicit none
integer,intent(in) :: n,np,piv(n)
complex(dp),intent(in) :: a(np,n)
integer xyzzyaaaa56
logical xyzzyaaab56
xyzzyaaab56=.false.
lu_logdet_cmplx=log(a(n,n))
do xyzzyaaaa56=1,n-1
lu_logdet_cmplx=lu_logdet_cmplx+log(a(xyzzyaaaa56,xyzzyaaaa56))
if(xyzzyaaaa56/=piv(xyzzyaaaa56))xyzzyaaab56=.not.xyzzyaaab56
enddo
if(xyzzyaaab56)lu_logdet_cmplx=lu_logdet_cmplx+cmplx(0.d0,pi,dp)
end function lu_logdet_cmplx
complex(dp) function lu_logdet_cmplx_dz(a,piv,np,n)
implicit none
integer,intent(in) :: n,np,piv(n)
real(dp),intent(in) :: a(np,np,2)
integer xyzzyaaaa57
logical xyzzyaaab57
xyzzyaaab57=.false.
lu_logdet_cmplx_dz=log(cmplx(a(n,n,1),a(n,n,2),dp))
do xyzzyaaaa57=1,n-1
lu_logdet_cmplx_dz=lu_logdet_cmplx_dz+log(cmplx(a(xyzzyaaaa57,xyzzyaaa&
&a57,1),a(xyzzyaaaa57,xyzzyaaaa57,2),dp))
if(xyzzyaaaa57/=piv(xyzzyaaaa57))xyzzyaaab57=.not.xyzzyaaab57
enddo
if(xyzzyaaab57)lu_logdet_cmplx_dz=lu_logdet_cmplx_dz+cmplx(0.d0,pi,dp)
end function lu_logdet_cmplx_dz
subroutine eigenproblem(a,n,eigenvals,eigenvecs,lerror)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: a(n,n)
complex(dp),intent(out) :: eigenvals(n),eigenvecs(n,n)
logical,intent(out) :: lerror
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58
real(dp) xyzzyaaae58(1,1)
real(dp),allocatable :: xyzzyaaaf58(:,:),xyzzyaaag58(:),xyzzyaaah58(:,&
&:),xyzzyaaai58(:,:)
xyzzyaaab58=0
allocate(xyzzyaaaf58(n,n),xyzzyaaah58(n,2),xyzzyaaai58(n,n),stat=xyzzy&
&aaac58)
call check_alloc(xyzzyaaac58,'EIGENVECTOR','a_work')
call dcopy(n*n,a(1,1),1,xyzzyaaaf58(1,1),1)
xyzzyaaah58=0.d0
xyzzyaaai58=0.d0
allocate(xyzzyaaag58(1),stat=xyzzyaaac58)
call check_alloc(xyzzyaaac58,'EIGENVALUES','work')
xyzzyaaag58=0.d0
xyzzyaaad58=-1
call dgeev('N','V',n,xyzzyaaaf58,n,xyzzyaaah58(1,1),xyzzyaaah58(1,2),x&
&yzzyaaae58,1,xyzzyaaai58,n,xyzzyaaag58,xyzzyaaad58,xyzzyaaab58)
if(xyzzyaaab58/=0)call errstop('EIGENVALUES','DGEEV says parameter #'/&
&/trim(i2s(-xyzzyaaab58))//' has an illegal value.')
xyzzyaaad58=nint(xyzzyaaag58(1))
deallocate(xyzzyaaag58)
allocate(xyzzyaaag58(xyzzyaaad58),stat=xyzzyaaac58)
call check_alloc(xyzzyaaac58,'EIGENVALUES','work')
xyzzyaaag58=0.d0
call dgeev('N','V',n,xyzzyaaaf58,n,xyzzyaaah58(1,1),xyzzyaaah58(1,2),x&
&yzzyaaae58,1,xyzzyaaai58,n,xyzzyaaag58,xyzzyaaad58,xyzzyaaab58)
lerror=xyzzyaaab58/=0
eigenvals(1:n)=cmplx(xyzzyaaah58(1:n,1),xyzzyaaah58(1:n,2),dp)
xyzzyaaaa58=1
do while(xyzzyaaaa58<=n)
if(xyzzyaaaa58<n.and.aimag(eigenvals(xyzzyaaaa58))/=0.d0)then
if(eigenvals(xyzzyaaaa58)==conjg(eigenvals(xyzzyaaaa58+1)))then
eigenvecs(1:n,xyzzyaaaa58)=cmplx(xyzzyaaai58(1:n,xyzzyaaaa58),xyzzyaaa&
&i58(1:n,xyzzyaaaa58+1),dp)
eigenvecs(1:n,xyzzyaaaa58+1)=cmplx(xyzzyaaai58(1:n,xyzzyaaaa58),-xyzzy&
&aaai58(1:n,xyzzyaaaa58+1),dp)
xyzzyaaaa58=xyzzyaaaa58+2
else
eigenvecs(1:n,xyzzyaaaa58)=cmplx(xyzzyaaai58(1:n,xyzzyaaaa58),0.d0,dp)
xyzzyaaaa58=xyzzyaaaa58+1
endif
else
eigenvecs(1:n,xyzzyaaaa58)=cmplx(xyzzyaaai58(1:n,xyzzyaaaa58),0.d0,dp)
xyzzyaaaa58=xyzzyaaaa58+1
endif
enddo
deallocate(xyzzyaaaf58,xyzzyaaag58,xyzzyaaah58,xyzzyaaai58)
end subroutine eigenproblem
subroutine sv_decom(m,n,a,w,v,lerror)
implicit none
integer,intent(in) :: m,n
real(dp),intent(inout) :: a(m,n)
real(dp),intent(out) :: w(n),v(n,n)
logical,intent(out) :: lerror
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59
integer,allocatable :: xyzzyaaad59(:)
real(dp) xyzzyaaae59(1,1)
real(dp),allocatable :: xyzzyaaaf59(:)
v=0.d0
w=0.d0
xyzzyaaac59=0
lerror=.false.
allocate(xyzzyaaad59(8*min(m,n)),xyzzyaaaf59(1),stat=xyzzyaaaa59)
call check_alloc(xyzzyaaaa59,'SV_DECOMPOSE','iwork')
xyzzyaaad59=0
xyzzyaaaf59=0.d0
xyzzyaaab59=-1
call dgesdd('O',m,n,a,m,w,xyzzyaaae59,1,v,n,xyzzyaaaf59,xyzzyaaab59,xy&
&zzyaaad59,xyzzyaaac59)
if(xyzzyaaac59/=0)call errstop('SV_DECOMPOSE','DGESDD says parameter #&
&'//trim(i2s(-xyzzyaaac59))//' has an illegal value.')
xyzzyaaab59=nint(xyzzyaaaf59(1))
deallocate(xyzzyaaaf59)
allocate(xyzzyaaaf59(xyzzyaaab59),stat=xyzzyaaaa59)
call check_alloc(xyzzyaaaa59,'SV_DECOMPOSE','work')
call dgesdd('O',m,n,a,m,w,xyzzyaaae59,1,v,n,xyzzyaaaf59,xyzzyaaab59,xy&
&zzyaaad59,xyzzyaaac59)
lerror=xyzzyaaac59/=0
v=transpose(v)
deallocate(xyzzyaaaf59,xyzzyaaad59)
end subroutine sv_decom
subroutine svd_inverse(m,n,u,w,v,ainv,maxw,nzeros)
implicit none
integer,intent(in) :: m,n
integer,intent(out) :: nzeros
real(dp),intent(in) :: u(m,n),v(n,n)
real(dp),intent(inout) :: w(n)
real(dp),intent(out) :: ainv(n,m),maxw
integer xyzzyaaaa60,xyzzyaaab60,xyzzyaaac60
real(dp) xyzzyaaad60(n),xyzzyaaae60,xyzzyaaaf60
real(dp),parameter :: xyzzyaaag60=1.d3
maxw=maxval(w)
xyzzyaaae60=xyzzyaaag60*epsilon(1.d0)*maxw
nzeros=0
do xyzzyaaaa60=1,n
if(w(xyzzyaaaa60)<xyzzyaaae60)then
w(xyzzyaaaa60)=0.d0
xyzzyaaad60(xyzzyaaaa60)=0.d0
nzeros=nzeros+1
else
xyzzyaaad60(xyzzyaaaa60)=1.d0/w(xyzzyaaaa60)
endif
enddo
do xyzzyaaab60=1,m
do xyzzyaaaa60=1,n
xyzzyaaaf60=0.d0
do xyzzyaaac60=1,n
xyzzyaaaf60=xyzzyaaaf60+u(xyzzyaaab60,xyzzyaaac60)*xyzzyaaad60(xyzzyaa&
&ac60)*v(xyzzyaaaa60,xyzzyaaac60)
enddo
ainv(xyzzyaaaa60,xyzzyaaab60)=xyzzyaaaf60
enddo
enddo
end subroutine svd_inverse
function inverse(m,n,a,is_singular,nzeros,lerror)
integer,intent(in) :: m,n
integer,intent(out),optional :: nzeros
real(dp),intent(in) :: a(m,n)
logical,intent(out),optional :: is_singular,lerror
integer xyzzyaaaa61,xyzzyaaab61
real(dp) inverse(n,m),xyzzyaaac61
real(dp),allocatable :: xyzzyaaad61(:,:),xyzzyaaae61(:),xyzzyaaaf61(:,&
&:)
logical xyzzyaaag61
allocate(xyzzyaaad61(m,n),xyzzyaaae61(n),xyzzyaaaf61(n,n),stat=xyzzyaa&
&ab61)
call check_alloc(xyzzyaaab61,'INVERSE','')
call dcopy(m*n,a(1,1),1,xyzzyaaad61(1,1),1)
call sv_decom(m,n,xyzzyaaad61(1,1),xyzzyaaae61(1),xyzzyaaaf61(1,1),xyz&
&zyaaag61)
if(present(lerror))then
lerror=xyzzyaaag61
else
if(xyzzyaaag61)call errstop('INVERSE','SV decomposition failed.')
endif
call svd_inverse(m,n,xyzzyaaad61(1,1),xyzzyaaae61(1),xyzzyaaaf61(1,1),&
&inverse,xyzzyaaac61,xyzzyaaaa61)
if(present(is_singular))is_singular=xyzzyaaaa61>0
if(present(nzeros))nzeros=xyzzyaaaa61
deallocate(xyzzyaaad61,xyzzyaaae61,xyzzyaaaf61)
end function inverse
subroutine dmatmul(b,c,a,m,p,n)
implicit none
integer,intent(in) :: m,n,p
real(dp),intent(in) :: b(m,p),c(p,n)
real(dp),intent(out) :: a(m,n)
call dgemm('N','N',m,n,p,1.d0,b,m,c,p,0.d0,a,m)
end subroutine dmatmul
real(dp) function hypothenuse(a,b)
implicit none
real(dp),intent(in) :: a,b
real(dp) xyzzyaaaa63,xyzzyaaab63,xyzzyaaac63,xyzzyaaad63,xyzzyaaae63
xyzzyaaaa63=abs(a)
xyzzyaaab63=abs(b)
xyzzyaaac63=max(xyzzyaaaa63,xyzzyaaab63)
if(xyzzyaaac63==0.d0)then
hypothenuse=0.d0
else
xyzzyaaad63=min(xyzzyaaaa63,xyzzyaaab63)
xyzzyaaae63=xyzzyaaad63/xyzzyaaac63
hypothenuse=xyzzyaaac63*sqrt(1.d0+xyzzyaaae63*xyzzyaaae63)
endif
end function hypothenuse
subroutine interp_nev(points,val,npoints,x,val_interp,err_interp)
implicit none
integer,intent(in) :: npoints
real(dp),intent(in) :: points(npoints),val(npoints),x
real(dp),intent(out) :: val_interp,err_interp
integer xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64
integer,parameter :: xyzzyaaad64=10
real(dp) xyzzyaaae64(xyzzyaaad64),xyzzyaaaf64(xyzzyaaad64),xyzzyaaag64&
&,xyzzyaaah64,xyzzyaaai64,xyzzyaaaj64
xyzzyaaag64=abs(x-points(1))
xyzzyaaac64=1
do xyzzyaaaa64=1,npoints
xyzzyaaah64=abs(x-points(xyzzyaaaa64))
if(xyzzyaaah64<xyzzyaaag64)then
xyzzyaaag64=xyzzyaaah64
xyzzyaaac64=xyzzyaaaa64
endif
xyzzyaaae64(xyzzyaaaa64)=val(xyzzyaaaa64)
xyzzyaaaf64(xyzzyaaaa64)=val(xyzzyaaaa64)
enddo
val_interp=val(xyzzyaaac64)
xyzzyaaac64=xyzzyaaac64-1
do xyzzyaaaa64=1,npoints-1
do xyzzyaaab64=1,npoints-xyzzyaaaa64
xyzzyaaag64=points(xyzzyaaab64)-x
xyzzyaaah64=points(xyzzyaaaa64+xyzzyaaab64)-x
xyzzyaaaj64=xyzzyaaae64(xyzzyaaab64+1)-xyzzyaaaf64(xyzzyaaab64)
xyzzyaaai64=xyzzyaaag64-xyzzyaaah64
xyzzyaaae64(xyzzyaaab64)=xyzzyaaag64*xyzzyaaaj64/xyzzyaaai64
xyzzyaaaf64(xyzzyaaab64)=xyzzyaaah64*xyzzyaaaj64/xyzzyaaai64
enddo
if(xyzzyaaac64*2<npoints-xyzzyaaaa64)then
err_interp=xyzzyaaae64(xyzzyaaac64+1)
else
err_interp=xyzzyaaaf64(xyzzyaaac64)
xyzzyaaac64=xyzzyaaac64-1
endif
val_interp=val_interp+err_interp
enddo
end subroutine interp_nev
subroutine interp_nev_with_derivs(points,val,npoints,x,val_interp,fd_i&
&nterp,sd_interp)
implicit none
integer,intent(in) :: npoints
real(dp),intent(in) :: points(npoints),val(npoints),x
real(dp),intent(out) :: val_interp,fd_interp,sd_interp
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65
integer,parameter :: xyzzyaaad65=10
real(dp) xyzzyaaae65,xyzzyaaaf65,xyzzyaaag65,xyzzyaaah65,xyzzyaaai65,x&
&yzzyaaaj65,xyzzyaaak65,xyzzyaaal65,xyzzyaaam65,xyzzyaaan65(xyzzyaaad6&
&5),xyzzyaaao65(xyzzyaaad65),xyzzyaaap65(xyzzyaaad65),xyzzyaaaq65(xyzz&
&yaaad65),xyzzyaaar65(xyzzyaaad65),xyzzyaaas65(xyzzyaaad65)
xyzzyaaae65=abs(x-points(1))
xyzzyaaac65=1
do xyzzyaaaa65=1,npoints
xyzzyaaaf65=abs(x-points(xyzzyaaaa65))
if(xyzzyaaaf65<xyzzyaaae65)then
xyzzyaaae65=xyzzyaaaf65
xyzzyaaac65=xyzzyaaaa65
endif
xyzzyaaan65(xyzzyaaaa65)=val(xyzzyaaaa65)
xyzzyaaao65(xyzzyaaaa65)=0.d0
xyzzyaaap65(xyzzyaaaa65)=0.d0
xyzzyaaaq65(xyzzyaaaa65)=val(xyzzyaaaa65)
xyzzyaaar65(xyzzyaaaa65)=0.d0
xyzzyaaas65(xyzzyaaaa65)=0.d0
enddo
val_interp=val(xyzzyaaac65)
fd_interp=0.d0
sd_interp=0.d0
xyzzyaaac65=xyzzyaaac65-1
do xyzzyaaaa65=1,npoints-1
do xyzzyaaab65=1,npoints-xyzzyaaaa65
xyzzyaaae65=points(xyzzyaaab65)-x
xyzzyaaaf65=points(xyzzyaaaa65+xyzzyaaab65)-x
xyzzyaaah65=xyzzyaaan65(xyzzyaaab65+1)-xyzzyaaaq65(xyzzyaaab65)
xyzzyaaai65=xyzzyaaao65(xyzzyaaab65+1)-xyzzyaaar65(xyzzyaaab65)
xyzzyaaaj65=xyzzyaaap65(xyzzyaaab65+1)-xyzzyaaas65(xyzzyaaab65)
xyzzyaaag65=xyzzyaaae65-xyzzyaaaf65
xyzzyaaan65(xyzzyaaab65)=xyzzyaaae65*xyzzyaaah65/xyzzyaaag65
xyzzyaaao65(xyzzyaaab65)=(-xyzzyaaah65+xyzzyaaae65*xyzzyaaai65)/xyzzya&
&aag65
xyzzyaaap65(xyzzyaaab65)=(-2.d0*xyzzyaaai65+xyzzyaaae65*xyzzyaaaj65)/x&
&yzzyaaag65
xyzzyaaaq65(xyzzyaaab65)=xyzzyaaaf65*xyzzyaaah65/xyzzyaaag65
xyzzyaaar65(xyzzyaaab65)=(-xyzzyaaah65+xyzzyaaaf65*xyzzyaaai65)/xyzzya&
&aag65
xyzzyaaas65(xyzzyaaab65)=(-2.d0*xyzzyaaai65+xyzzyaaaf65*xyzzyaaaj65)/x&
&yzzyaaag65
enddo
if(xyzzyaaac65*2<npoints-xyzzyaaaa65)then
xyzzyaaak65=xyzzyaaan65(xyzzyaaac65+1)
xyzzyaaal65=xyzzyaaao65(xyzzyaaac65+1)
xyzzyaaam65=xyzzyaaap65(xyzzyaaac65+1)
else
xyzzyaaak65=xyzzyaaaq65(xyzzyaaac65)
xyzzyaaal65=xyzzyaaar65(xyzzyaaac65)
xyzzyaaam65=xyzzyaaas65(xyzzyaaac65)
xyzzyaaac65=xyzzyaaac65-1
endif
val_interp=val_interp+xyzzyaaak65
fd_interp=fd_interp+xyzzyaaal65
sd_interp=sd_interp+xyzzyaaam65
enddo
end subroutine interp_nev_with_derivs
subroutine lookup(xarr,n,x,jl)
implicit none
integer,intent(in) :: n
integer,intent(out) :: jl
real(dp),intent(in) :: xarr(*),x
integer xyzzyaaaa66,xyzzyaaab66
jl=0
xyzzyaaab66=n+1
if(xarr(n)>xarr(1))then
do while(xyzzyaaab66-jl>1)
xyzzyaaaa66=(xyzzyaaab66+jl)/2
if(x>xarr(xyzzyaaaa66))then
jl=xyzzyaaaa66
else
xyzzyaaab66=xyzzyaaaa66
endif
enddo
else
do while(xyzzyaaab66-jl>1)
xyzzyaaaa66=(xyzzyaaab66+jl)/2
if(x<=xarr(xyzzyaaaa66))then
jl=xyzzyaaaa66
else
xyzzyaaab66=xyzzyaaaa66
endif
enddo
endif
end subroutine lookup
subroutine interp_nev_2d_with_derivs(x1a,x2a,ya,m,n,x1,x2,y,y_d1,y_d2,&
&y_dd1,y_dd2,y_d1d2)
implicit none
integer,intent(in) :: m,n
real(dp),intent(in) :: x1,x2,x1a(m),x2a(n),ya(m,n)
real(dp),intent(out) :: y,y_d1,y_d2,y_dd1,y_dd2,y_d1d2
integer xyzzyaaaa67,xyzzyaaab67
integer,parameter :: xyzzyaaac67=10,xyzzyaaad67=10
real(dp) xyzzyaaae67(xyzzyaaad67),xyzzyaaaf67(xyzzyaaac67),xyzzyaaag67&
&(xyzzyaaad67),xyzzyaaah67(xyzzyaaad67),xyzzyaaai67,xyzzyaaaj67,xyzzya&
&aak67
do xyzzyaaaa67=1,m
do xyzzyaaab67=1,n
xyzzyaaaf67(xyzzyaaab67)=ya(xyzzyaaaa67,xyzzyaaab67)
enddo
call interp_nev_with_derivs(x2a,xyzzyaaaf67,n,x2,xyzzyaaae67(xyzzyaaaa&
&67),xyzzyaaag67(xyzzyaaaa67),xyzzyaaah67(xyzzyaaaa67))
enddo
call interp_nev_with_derivs(x1a,xyzzyaaae67,n,x1,xyzzyaaai67,xyzzyaaaj&
&67,xyzzyaaak67)
y=xyzzyaaai67
y_d1=xyzzyaaaj67
y_dd1=xyzzyaaak67
call interp_nev_with_derivs(x1a,xyzzyaaag67,n,x1,xyzzyaaai67,xyzzyaaaj&
&67,xyzzyaaak67)
y_d2=xyzzyaaai67
y_d1d2=xyzzyaaaj67
call interp_nev_with_derivs(x1a,xyzzyaaah67,n,x1,xyzzyaaai67,xyzzyaaaj&
&67,xyzzyaaak67)
y_dd2=xyzzyaaai67
end subroutine interp_nev_2d_with_derivs
subroutine dadd(n,x,incx,y,incy)
implicit none
integer,intent(in) :: incx,incy,n
real(dp),intent(in) :: x(*)
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68
if(n<=0)return
if(incx==1.and.incy==1)then
y(1:n)=y(1:n)+x(1:n)
else
if(incx>=0)then
xyzzyaaab68=1
else
xyzzyaaab68=(-n+1)*incx+1
endif
if(incy>=0)then
xyzzyaaac68=1
else
xyzzyaaac68=(-n+1)*incy+1
endif
do xyzzyaaaa68=1,n
y(xyzzyaaac68)=y(xyzzyaaac68)+x(xyzzyaaab68)
xyzzyaaab68=xyzzyaaab68+incx
xyzzyaaac68=xyzzyaaac68+incy
enddo
endif
end subroutine dadd
subroutine zadd(n,x,incx,y,incy)
implicit none
integer,intent(in) :: incx,incy,n
complex(dp),intent(in) :: x(*)
complex(dp),intent(inout) :: y(*)
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69
if(n<=0)return
if(incx==1.and.incy==1)then
y(1:n)=y(1:n)+x(1:n)
else
if(incx>=0)then
xyzzyaaab69=1
else
xyzzyaaab69=(-n+1)*incx+1
endif
if(incy>=0)then
xyzzyaaac69=1
else
xyzzyaaac69=(-n+1)*incy+1
endif
do xyzzyaaaa69=1,n
y(xyzzyaaac69)=y(xyzzyaaac69)+x(xyzzyaaab69)
xyzzyaaab69=xyzzyaaab69+incx
xyzzyaaac69=xyzzyaaac69+incy
enddo
endif
end subroutine zadd
subroutine dcopy3(n,n3,x,incx,y,incy)
implicit none
integer,intent(in) :: n,n3,incx,incy
real(dp),intent(in) :: x(*)
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70
if(n<1)return
if(incx==1.and.incy==1)then
y(1:n3)=x(1:n3)
else
xyzzyaaab70=1
xyzzyaaac70=1
do xyzzyaaaa70=1,n
y(xyzzyaaac70)=x(xyzzyaaab70)
xyzzyaaab70=xyzzyaaab70+1
xyzzyaaac70=xyzzyaaac70+1
y(xyzzyaaac70)=x(xyzzyaaab70)
xyzzyaaab70=xyzzyaaab70+1
xyzzyaaac70=xyzzyaaac70+1
y(xyzzyaaac70)=x(xyzzyaaab70)
xyzzyaaab70=xyzzyaaab70+incx
xyzzyaaac70=xyzzyaaac70+incy
enddo
endif
end subroutine dcopy3
subroutine dcopy6(n,n6,x,incx,y,incy)
implicit none
integer,intent(in) :: n,n6,incx,incy
real(dp),intent(in) :: x(*)
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71
if(n<1)return
if(incx==1.and.incy==1)then
y(1:n6)=x(1:n6)
else
xyzzyaaab71=1
xyzzyaaac71=1
do xyzzyaaaa71=1,n
y(xyzzyaaac71)=x(xyzzyaaab71)
xyzzyaaab71=xyzzyaaab71+1
xyzzyaaac71=xyzzyaaac71+1
y(xyzzyaaac71)=x(xyzzyaaab71)
xyzzyaaab71=xyzzyaaab71+1
xyzzyaaac71=xyzzyaaac71+1
y(xyzzyaaac71)=x(xyzzyaaab71)
xyzzyaaab71=xyzzyaaab71+1
xyzzyaaac71=xyzzyaaac71+1
y(xyzzyaaac71)=x(xyzzyaaab71)
xyzzyaaab71=xyzzyaaab71+1
xyzzyaaac71=xyzzyaaac71+1
y(xyzzyaaac71)=x(xyzzyaaab71)
xyzzyaaab71=xyzzyaaab71+1
xyzzyaaac71=xyzzyaaac71+1
y(xyzzyaaac71)=x(xyzzyaaab71)
xyzzyaaab71=xyzzyaaab71+incx
xyzzyaaac71=xyzzyaaac71+incy
enddo
endif
end subroutine dcopy6
subroutine dcopy_ee(n,ii,x,y)
implicit none
integer,intent(in) :: n,ii
real(dp),intent(in) :: x(*)
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa72,xyzzyaaab72,xyzzyaaac72,xyzzyaaad72,xyzzyaaae72
xyzzyaaad72=4*n
xyzzyaaae72=4*ii-4
xyzzyaaac72=1+xyzzyaaae72
xyzzyaaab72=1
do xyzzyaaaa72=1,ii-1
y(xyzzyaaac72)=-x(xyzzyaaab72)
y(xyzzyaaac72+1)=-x(xyzzyaaab72+1)
y(xyzzyaaac72+2)=-x(xyzzyaaab72+2)
y(xyzzyaaac72+3)=x(xyzzyaaab72+3)
xyzzyaaac72=xyzzyaaac72+xyzzyaaad72
xyzzyaaab72=xyzzyaaab72+4
enddo
y(xyzzyaaac72-xyzzyaaae72:xyzzyaaac72+xyzzyaaad72-xyzzyaaae72-1)=x(1:x&
&yzzyaaad72)
xyzzyaaab72=xyzzyaaab72+4
xyzzyaaac72=xyzzyaaac72+xyzzyaaad72
do xyzzyaaaa72=ii+1,n
y(xyzzyaaac72)=-x(xyzzyaaab72)
y(xyzzyaaac72+1)=-x(xyzzyaaab72+1)
y(xyzzyaaac72+2)=-x(xyzzyaaab72+2)
y(xyzzyaaac72+3)=x(xyzzyaaab72+3)
xyzzyaaac72=xyzzyaaac72+xyzzyaaad72
xyzzyaaab72=xyzzyaaab72+4
enddo
end subroutine dcopy_ee
subroutine dcopy_ee_flip(n,ii,x,y)
implicit none
integer,intent(in) :: n,ii
real(dp),intent(in) :: x(*)
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa73,xyzzyaaab73,xyzzyaaac73,xyzzyaaad73
xyzzyaaad73=4*n
xyzzyaaac73=1+4*ii-4
xyzzyaaab73=1
do xyzzyaaaa73=1,ii-1
y(xyzzyaaac73)=-x(xyzzyaaab73)
y(xyzzyaaac73+1)=-x(xyzzyaaab73+1)
y(xyzzyaaac73+2)=-x(xyzzyaaab73+2)
y(xyzzyaaac73+3)=x(xyzzyaaab73+3)
xyzzyaaac73=xyzzyaaac73+xyzzyaaad73
xyzzyaaab73=xyzzyaaab73+4
enddo
xyzzyaaab73=xyzzyaaab73+4
xyzzyaaac73=xyzzyaaac73+xyzzyaaad73
do xyzzyaaaa73=ii+1,n
y(xyzzyaaac73)=-x(xyzzyaaab73)
y(xyzzyaaac73+1)=-x(xyzzyaaab73+1)
y(xyzzyaaac73+2)=-x(xyzzyaaab73+2)
y(xyzzyaaac73+3)=x(xyzzyaaab73+3)
xyzzyaaac73=xyzzyaaac73+xyzzyaaad73
xyzzyaaab73=xyzzyaaab73+4
enddo
end subroutine dcopy_ee_flip
subroutine xyzzyaaad1(n,xr,xi,incx,y,incy)
implicit none
integer,intent(in) :: n,incx,incy
real(dp),intent(in) :: xr(*),xi(*)
complex(dp),intent(inout) :: y(*)
integer xyzzyaaaa74,xyzzyaaab74,xyzzyaaac74
if(n<1)return
if(incx==1.and.incy==1)then
y(1:n)=cmplx(xr(1:n),xi(1:n),dp)
else
xyzzyaaab74=1
xyzzyaaac74=1
do xyzzyaaaa74=1,n
y(xyzzyaaac74)=cmplx(xr(xyzzyaaab74),xi(xyzzyaaab74),dp)
xyzzyaaab74=xyzzyaaab74+incx
xyzzyaaac74=xyzzyaaac74+incy
enddo
endif
end subroutine xyzzyaaad1
subroutine xyzzyaaae1(n,x,incx,yr,yi,incy)
implicit none
integer,intent(in) :: n,incx,incy
real(dp),intent(inout) :: yr(*),yi(*)
complex(dp),intent(in) :: x(*)
integer xyzzyaaaa75,xyzzyaaab75,xyzzyaaac75
if(n<1)return
if(incx==1.and.incy==1)then
yr(1:n)=dble(x(1:n))
yi(1:n)=aimag(x(1:n))
else
xyzzyaaab75=1
xyzzyaaac75=1
do xyzzyaaaa75=1,n
yr(xyzzyaaac75)=dble(x(xyzzyaaab75))
yi(xyzzyaaac75)=aimag(x(xyzzyaaab75))
xyzzyaaab75=xyzzyaaab75+incx
xyzzyaaac75=xyzzyaaac75+incy
enddo
endif
end subroutine xyzzyaaae1
subroutine xyzzyaaaf1(x,y)
implicit none
real(dp),intent(inout) :: x,y
real(dp) xyzzyaaaa76
xyzzyaaaa76=x
x=y
y=xyzzyaaaa76
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(x,y)
implicit none
complex(dp),intent(inout) :: x,y
complex(dp) xyzzyaaaa77
xyzzyaaaa77=x
x=y
y=xyzzyaaaa77
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(i,j)
implicit none
integer,intent(inout) :: i,j
integer xyzzyaaaa78
xyzzyaaaa78=i
i=j
j=xyzzyaaaa78
end subroutine xyzzyaaah1
subroutine xyzzyaaai1(a,b)
implicit none
character(*),intent(inout) :: a,b
character(max(len_trim(a),len_trim(b))) c
c=trim(a)
a=trim(b)
b=trim(c)
end subroutine xyzzyaaai1
subroutine iswap(n,x,incx,y,incy)
implicit none
integer,intent(in) :: n,incx,incy
integer,intent(inout) :: x(*),y(*)
integer xyzzyaaaa80,xyzzyaaab80,xyzzyaaac80,xyzzyaaad80
if(n<1)return
xyzzyaaab80=1
xyzzyaaac80=1
do xyzzyaaaa80=1,n
xyzzyaaad80=y(xyzzyaaac80)
y(xyzzyaaac80)=x(xyzzyaaab80)
x(xyzzyaaab80)=xyzzyaaad80
xyzzyaaab80=xyzzyaaab80+incx
xyzzyaaac80=xyzzyaaac80+incy
enddo
end subroutine iswap
function dsum(n,x,incx) result(s)
implicit none
integer,intent(in) :: n,incx
real(dp),intent(in) :: x(*)
integer xyzzyaaaa81,xyzzyaaab81
real(dp) s
s=0.d0
if(n<1)return
if(incx==1)then
s=sum(x(1:n))
else
xyzzyaaab81=1
do xyzzyaaaa81=1,n
s=s+x(xyzzyaaab81)
xyzzyaaab81=xyzzyaaab81+incx
enddo
endif
end function dsum
function dsum3(n,x,incx) result(s)
implicit none
integer,intent(in) :: n,incx
real(dp),intent(in) :: x(*)
integer xyzzyaaaa82,xyzzyaaab82
real(dp) s(3),xyzzyaaac82,xyzzyaaad82,xyzzyaaae82
s=0.d0
if(n<1)return
xyzzyaaac82=0.d0
xyzzyaaad82=0.d0
xyzzyaaae82=0.d0
xyzzyaaab82=1
do xyzzyaaaa82=1,n
xyzzyaaac82=xyzzyaaac82+x(xyzzyaaab82)
xyzzyaaab82=xyzzyaaab82+1
xyzzyaaad82=xyzzyaaad82+x(xyzzyaaab82)
xyzzyaaab82=xyzzyaaab82+1
xyzzyaaae82=xyzzyaaae82+x(xyzzyaaab82)
xyzzyaaab82=xyzzyaaab82+incx
enddo
s(1)=xyzzyaaac82
s(2)=xyzzyaaad82
s(3)=xyzzyaaae82
end function dsum3
real(dp) function ddot_s(m,rmap,x,incx,y,incy)
implicit none
integer,intent(in) :: m,rmap(*),incx,incy
real(dp),intent(in) :: x(*),y(*)
integer xyzzyaaaa83,xyzzyaaab83,xyzzyaaac83,xyzzyaaad83
ddot_s=0.d0
if(incx==1.and.incy==1)then
do xyzzyaaaa83=1,m
xyzzyaaab83=rmap(xyzzyaaaa83)
ddot_s=ddot_s+x(xyzzyaaab83)*y(xyzzyaaab83)
enddo
elseif(incy==1)then
do xyzzyaaaa83=1,m
xyzzyaaab83=rmap(xyzzyaaaa83)
xyzzyaaac83=1+(xyzzyaaab83-1)*incx
ddot_s=ddot_s+x(xyzzyaaac83)*y(xyzzyaaab83)
enddo
elseif(incx==1)then
do xyzzyaaaa83=1,m
xyzzyaaab83=rmap(xyzzyaaaa83)
xyzzyaaad83=1+(xyzzyaaab83-1)*incy
ddot_s=ddot_s+x(xyzzyaaab83)*y(xyzzyaaad83)
enddo
else
do xyzzyaaaa83=1,m
xyzzyaaab83=rmap(xyzzyaaaa83)
xyzzyaaac83=1+(xyzzyaaab83-1)*incx
xyzzyaaad83=1+(xyzzyaaab83-1)*incy
ddot_s=ddot_s+x(xyzzyaaac83)*y(xyzzyaaad83)
enddo
endif
end function ddot_s
real(dp) function ddot_alt(n,dx,size_dx,start_dx,incx,dy)
integer,intent(in) :: n,size_dx,start_dx,incx
real(dp) dx(size_dx),dy(*)
integer xyzzyaaaa84,xyzzyaaab84
ddot_alt=0.d0
xyzzyaaab84=start_dx
do xyzzyaaaa84=1,n
ddot_alt=ddot_alt+dx(xyzzyaaab84)*dy(xyzzyaaaa84)
xyzzyaaab84=xyzzyaaab84+incx
enddo
end function ddot_alt
complex(dp) function zdotu_alt(n,zx,size_zx,start_zx,incx,zy)
integer n,size_zx,start_zx,incx
complex(dp) zx(size_zx),zy(*)
integer xyzzyaaaa85,xyzzyaaab85
zdotu_alt=(0.0d0,0.0d0)
xyzzyaaab85=start_zx
do xyzzyaaaa85=1,n
zdotu_alt=zdotu_alt+zx(xyzzyaaab85)*zy(xyzzyaaaa85)
xyzzyaaab85=xyzzyaaab85+incx
enddo
end function zdotu_alt
function multi_ddot_s(nx,ny,rmap,mask,x,size_x,ldx,y,s) result(p)
implicit none
integer,intent(in) :: nx,ny,rmap(*),ldx,size_x
real(dp),intent(in) :: x(size_x),y(*),s
logical,intent(in) :: mask(nx)
integer xyzzyaaaa86,xyzzyaaab86,xyzzyaaac86,xyzzyaaad86,xyzzyaaae86,xy&
&zzyaaaf86,xyzzyaaag86,xyzzyaaah86,xyzzyaaai86,xyzzyaaaj86
integer,parameter :: xyzzyaaak86=6
real(dp) xyzzyaaal86,xyzzyaaam86,xyzzyaaan86,xyzzyaaao86,xyzzyaaap86,x&
&yzzyaaaq86,p(nx)
p=0.d0
xyzzyaaad86=ny-xyzzyaaak86+1
do xyzzyaaab86=1,xyzzyaaad86,xyzzyaaak86
xyzzyaaac86=rmap(xyzzyaaab86)
xyzzyaaae86=(xyzzyaaac86-1)*ldx
xyzzyaaal86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+1)
xyzzyaaaf86=(xyzzyaaac86-1)*ldx
xyzzyaaam86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+2)
xyzzyaaag86=(xyzzyaaac86-1)*ldx
xyzzyaaan86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+3)
xyzzyaaah86=(xyzzyaaac86-1)*ldx
xyzzyaaao86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+4)
xyzzyaaai86=(xyzzyaaac86-1)*ldx
xyzzyaaap86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+5)
xyzzyaaaj86=(xyzzyaaac86-1)*ldx
xyzzyaaaq86=y(xyzzyaaac86)
do xyzzyaaaa86=1,nx
if(mask(xyzzyaaaa86))p(xyzzyaaaa86)=p(xyzzyaaaa86)+xyzzyaaal86*x(xyzzy&
&aaae86+xyzzyaaaa86)+xyzzyaaam86*x(xyzzyaaaf86+xyzzyaaaa86)+xyzzyaaan8&
&6*x(xyzzyaaag86+xyzzyaaaa86)+xyzzyaaao86*x(xyzzyaaah86+xyzzyaaaa86)+x&
&yzzyaaap86*x(xyzzyaaai86+xyzzyaaaa86)+xyzzyaaaq86*x(xyzzyaaaj86+xyzzy&
&aaaa86)
enddo
enddo
select case(ny-xyzzyaaab86)
case(4)
xyzzyaaac86=rmap(xyzzyaaab86)
xyzzyaaae86=(xyzzyaaac86-1)*ldx
xyzzyaaal86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+1)
xyzzyaaaf86=(xyzzyaaac86-1)*ldx
xyzzyaaam86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+2)
xyzzyaaag86=(xyzzyaaac86-1)*ldx
xyzzyaaan86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+3)
xyzzyaaah86=(xyzzyaaac86-1)*ldx
xyzzyaaao86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+4)
xyzzyaaai86=(xyzzyaaac86-1)*ldx
xyzzyaaap86=y(xyzzyaaac86)
do xyzzyaaaa86=1,nx
if(mask(xyzzyaaaa86))p(xyzzyaaaa86)=p(xyzzyaaaa86)+xyzzyaaal86*x(xyzzy&
&aaae86+xyzzyaaaa86)+xyzzyaaam86*x(xyzzyaaaf86+xyzzyaaaa86)+xyzzyaaan8&
&6*x(xyzzyaaag86+xyzzyaaaa86)+xyzzyaaao86*x(xyzzyaaah86+xyzzyaaaa86)+x&
&yzzyaaap86*x(xyzzyaaai86+xyzzyaaaa86)
enddo
case(3)
xyzzyaaac86=rmap(xyzzyaaab86)
xyzzyaaae86=(xyzzyaaac86-1)*ldx
xyzzyaaal86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+1)
xyzzyaaaf86=(xyzzyaaac86-1)*ldx
xyzzyaaam86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+2)
xyzzyaaag86=(xyzzyaaac86-1)*ldx
xyzzyaaan86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+3)
xyzzyaaah86=(xyzzyaaac86-1)*ldx
xyzzyaaao86=y(xyzzyaaac86)
do xyzzyaaaa86=1,nx
if(mask(xyzzyaaaa86))p(xyzzyaaaa86)=p(xyzzyaaaa86)+xyzzyaaal86*x(xyzzy&
&aaae86+xyzzyaaaa86)+xyzzyaaam86*x(xyzzyaaaf86+xyzzyaaaa86)+xyzzyaaan8&
&6*x(xyzzyaaag86+xyzzyaaaa86)+xyzzyaaao86*x(xyzzyaaah86+xyzzyaaaa86)
enddo
case(2)
xyzzyaaac86=rmap(xyzzyaaab86)
xyzzyaaae86=(xyzzyaaac86-1)*ldx
xyzzyaaal86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+1)
xyzzyaaaf86=(xyzzyaaac86-1)*ldx
xyzzyaaam86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+2)
xyzzyaaag86=(xyzzyaaac86-1)*ldx
xyzzyaaan86=y(xyzzyaaac86)
do xyzzyaaaa86=1,nx
if(mask(xyzzyaaaa86))p(xyzzyaaaa86)=p(xyzzyaaaa86)+xyzzyaaal86*x(xyzzy&
&aaae86+xyzzyaaaa86)+xyzzyaaam86*x(xyzzyaaaf86+xyzzyaaaa86)+xyzzyaaan8&
&6*x(xyzzyaaag86+xyzzyaaaa86)
enddo
case(1)
xyzzyaaac86=rmap(xyzzyaaab86)
xyzzyaaae86=(xyzzyaaac86-1)*ldx
xyzzyaaal86=y(xyzzyaaac86)
xyzzyaaac86=rmap(xyzzyaaab86+1)
xyzzyaaaf86=(xyzzyaaac86-1)*ldx
xyzzyaaam86=y(xyzzyaaac86)
do xyzzyaaaa86=1,nx
if(mask(xyzzyaaaa86))p(xyzzyaaaa86)=p(xyzzyaaaa86)+xyzzyaaal86*x(xyzzy&
&aaae86+xyzzyaaaa86)+xyzzyaaam86*x(xyzzyaaaf86+xyzzyaaaa86)
enddo
case(0)
xyzzyaaac86=rmap(xyzzyaaab86)
xyzzyaaae86=(xyzzyaaac86-1)*ldx
xyzzyaaal86=y(xyzzyaaac86)
do xyzzyaaaa86=1,nx
if(mask(xyzzyaaaa86))p(xyzzyaaaa86)=p(xyzzyaaaa86)+xyzzyaaal86*x(xyzzy&
&aaae86+xyzzyaaaa86)
enddo
end select
call dscal(nx,s,p,1)
end function multi_ddot_s
function multi_zdotu_s(nx,ny,rmap,mask,x,size_x,ldx,y,s) result(p)
implicit none
integer,intent(in) :: nx,ny,rmap(*),size_x,ldx
real(dp),intent(in) :: s
complex(dp),intent(in) :: x(size_x),y(*)
logical,intent(in) :: mask(nx)
integer xyzzyaaaa87,xyzzyaaab87,xyzzyaaac87,xyzzyaaad87,xyzzyaaae87,xy&
&zzyaaaf87,xyzzyaaag87,xyzzyaaah87,xyzzyaaai87,xyzzyaaaj87
integer,parameter :: xyzzyaaak87=6
complex(dp) xyzzyaaal87,xyzzyaaam87,xyzzyaaan87,xyzzyaaao87,xyzzyaaap8&
&7,xyzzyaaaq87,p(nx)
p=czero
xyzzyaaad87=ny-xyzzyaaak87+1
do xyzzyaaab87=1,xyzzyaaad87,xyzzyaaak87
xyzzyaaac87=rmap(xyzzyaaab87)
xyzzyaaae87=(xyzzyaaac87-1)*ldx
xyzzyaaal87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+1)
xyzzyaaaf87=(xyzzyaaac87-1)*ldx
xyzzyaaam87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+2)
xyzzyaaag87=(xyzzyaaac87-1)*ldx
xyzzyaaan87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+3)
xyzzyaaah87=(xyzzyaaac87-1)*ldx
xyzzyaaao87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+4)
xyzzyaaai87=(xyzzyaaac87-1)*ldx
xyzzyaaap87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+5)
xyzzyaaaj87=(xyzzyaaac87-1)*ldx
xyzzyaaaq87=y(xyzzyaaac87)
do xyzzyaaaa87=1,nx
if(mask(xyzzyaaaa87))p(xyzzyaaaa87)=p(xyzzyaaaa87)+xyzzyaaal87*x(xyzzy&
&aaae87+xyzzyaaaa87)+xyzzyaaam87*x(xyzzyaaaf87+xyzzyaaaa87)+xyzzyaaan8&
&7*x(xyzzyaaag87+xyzzyaaaa87)+xyzzyaaao87*x(xyzzyaaah87+xyzzyaaaa87)+x&
&yzzyaaap87*x(xyzzyaaai87+xyzzyaaaa87)+xyzzyaaaq87*x(xyzzyaaaj87+xyzzy&
&aaaa87)
enddo
enddo
select case(ny-xyzzyaaab87)
case(4)
xyzzyaaac87=rmap(xyzzyaaab87)
xyzzyaaae87=(xyzzyaaac87-1)*ldx
xyzzyaaal87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+1)
xyzzyaaaf87=(xyzzyaaac87-1)*ldx
xyzzyaaam87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+2)
xyzzyaaag87=(xyzzyaaac87-1)*ldx
xyzzyaaan87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+3)
xyzzyaaah87=(xyzzyaaac87-1)*ldx
xyzzyaaao87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+4)
xyzzyaaai87=(xyzzyaaac87-1)*ldx
xyzzyaaap87=y(xyzzyaaac87)
do xyzzyaaaa87=1,nx
if(mask(xyzzyaaaa87))p(xyzzyaaaa87)=p(xyzzyaaaa87)+xyzzyaaal87*x(xyzzy&
&aaae87+xyzzyaaaa87)+xyzzyaaam87*x(xyzzyaaaf87+xyzzyaaaa87)+xyzzyaaan8&
&7*x(xyzzyaaag87+xyzzyaaaa87)+xyzzyaaao87*x(xyzzyaaah87+xyzzyaaaa87)+x&
&yzzyaaap87*x(xyzzyaaai87+xyzzyaaaa87)
enddo
case(3)
xyzzyaaac87=rmap(xyzzyaaab87)
xyzzyaaae87=(xyzzyaaac87-1)*ldx
xyzzyaaal87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+1)
xyzzyaaaf87=(xyzzyaaac87-1)*ldx
xyzzyaaam87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+2)
xyzzyaaag87=(xyzzyaaac87-1)*ldx
xyzzyaaan87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+3)
xyzzyaaah87=(xyzzyaaac87-1)*ldx
xyzzyaaao87=y(xyzzyaaac87)
do xyzzyaaaa87=1,nx
if(mask(xyzzyaaaa87))p(xyzzyaaaa87)=p(xyzzyaaaa87)+xyzzyaaal87*x(xyzzy&
&aaae87+xyzzyaaaa87)+xyzzyaaam87*x(xyzzyaaaf87+xyzzyaaaa87)+xyzzyaaan8&
&7*x(xyzzyaaag87+xyzzyaaaa87)+xyzzyaaao87*x(xyzzyaaah87+xyzzyaaaa87)
enddo
case(2)
xyzzyaaac87=rmap(xyzzyaaab87)
xyzzyaaae87=(xyzzyaaac87-1)*ldx
xyzzyaaal87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+1)
xyzzyaaaf87=(xyzzyaaac87-1)*ldx
xyzzyaaam87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+2)
xyzzyaaag87=(xyzzyaaac87-1)*ldx
xyzzyaaan87=y(xyzzyaaac87)
do xyzzyaaaa87=1,nx
if(mask(xyzzyaaaa87))p(xyzzyaaaa87)=p(xyzzyaaaa87)+xyzzyaaal87*x(xyzzy&
&aaae87+xyzzyaaaa87)+xyzzyaaam87*x(xyzzyaaaf87+xyzzyaaaa87)+xyzzyaaan8&
&7*x(xyzzyaaag87+xyzzyaaaa87)
enddo
case(1)
xyzzyaaac87=rmap(xyzzyaaab87)
xyzzyaaae87=(xyzzyaaac87-1)*ldx
xyzzyaaal87=y(xyzzyaaac87)
xyzzyaaac87=rmap(xyzzyaaab87+1)
xyzzyaaaf87=(xyzzyaaac87-1)*ldx
xyzzyaaam87=y(xyzzyaaac87)
do xyzzyaaaa87=1,nx
if(mask(xyzzyaaaa87))p(xyzzyaaaa87)=p(xyzzyaaaa87)+xyzzyaaal87*x(xyzzy&
&aaae87+xyzzyaaaa87)+xyzzyaaam87*x(xyzzyaaaf87+xyzzyaaaa87)
enddo
case(0)
xyzzyaaac87=rmap(xyzzyaaab87)
xyzzyaaae87=(xyzzyaaac87-1)*ldx
xyzzyaaal87=y(xyzzyaaac87)
do xyzzyaaaa87=1,nx
if(mask(xyzzyaaaa87))p(xyzzyaaaa87)=p(xyzzyaaaa87)+xyzzyaaal87*x(xyzzy&
&aaae87+xyzzyaaaa87)
enddo
end select
call dzscal(nx,s,p,1)
end function multi_zdotu_s
subroutine d_sparse_matrix_vector(m,alpha,a,lda,ias,ldias,inums,x,incx&
&,y,incy)
implicit none
integer,intent(in) :: m,lda,ldias,ias(ldias,*),inums(lda),incx,incy
real(dp),intent(in) :: a(lda,*),x(*),alpha
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa88,xyzzyaaab88,xyzzyaaac88
real(dp) xyzzyaaad88
if(incx==1.and.incy==1)then
do xyzzyaaaa88=1,m
xyzzyaaad88=0.d0
do xyzzyaaab88=1,inums(xyzzyaaaa88)
xyzzyaaac88=ias(xyzzyaaab88,xyzzyaaaa88)
xyzzyaaad88=xyzzyaaad88+a(xyzzyaaaa88,xyzzyaaac88)*x(xyzzyaaac88)
enddo
y(xyzzyaaaa88)=alpha*xyzzyaaad88
enddo
elseif(incy==1)then
do xyzzyaaaa88=1,m
xyzzyaaad88=0.d0
do xyzzyaaab88=1,inums(xyzzyaaaa88)
xyzzyaaac88=ias(xyzzyaaab88,xyzzyaaaa88)
xyzzyaaad88=xyzzyaaad88+a(xyzzyaaaa88,xyzzyaaac88)*x(1+(xyzzyaaac88-1)&
&*incx)
enddo
y(xyzzyaaaa88)=alpha*xyzzyaaad88
enddo
elseif(incx==1)then
do xyzzyaaaa88=1,m
xyzzyaaad88=0.d0
do xyzzyaaab88=1,inums(xyzzyaaaa88)
xyzzyaaac88=ias(xyzzyaaab88,xyzzyaaaa88)
xyzzyaaad88=xyzzyaaad88+a(xyzzyaaaa88,xyzzyaaac88)*x(xyzzyaaac88)
enddo
y(1+(xyzzyaaaa88-1)*incy)=alpha*xyzzyaaad88
enddo
else
do xyzzyaaaa88=1,m
xyzzyaaad88=0.d0
do xyzzyaaab88=1,inums(xyzzyaaaa88)
xyzzyaaac88=ias(xyzzyaaab88,xyzzyaaaa88)
xyzzyaaad88=xyzzyaaad88+a(xyzzyaaaa88,xyzzyaaac88)*x(1+(xyzzyaaac88-1)&
&*incx)
enddo
y(1+(xyzzyaaaa88-1)*incy)=alpha*xyzzyaaad88
enddo
endif
end subroutine d_sparse_matrix_vector
subroutine d2_sparse_matrix_vector(m,alpha,a,lda,ias,ldias,inums,xr,xi&
&,incx,yr,yi,incy)
implicit none
integer,intent(in) :: m,lda,ldias,ias(ldias,*),inums(lda),incx,incy
real(dp),intent(in) :: a(lda,*),xr(*),xi(*),alpha
real(dp),intent(inout) :: yr(*),yi(*)
integer xyzzyaaaa89,xyzzyaaab89,xyzzyaaac89
real(dp) xyzzyaaad89,xyzzyaaae89
if(incx==1.and.incy==1)then
do xyzzyaaaa89=1,m
xyzzyaaad89=0.d0
xyzzyaaae89=0.d0
do xyzzyaaab89=1,inums(xyzzyaaaa89)
xyzzyaaac89=ias(xyzzyaaab89,xyzzyaaaa89)
xyzzyaaad89=xyzzyaaad89+a(xyzzyaaaa89,xyzzyaaac89)*xr(xyzzyaaac89)
xyzzyaaae89=xyzzyaaae89+a(xyzzyaaaa89,xyzzyaaac89)*xi(xyzzyaaac89)
enddo
yr(xyzzyaaaa89)=alpha*xyzzyaaad89
yi(xyzzyaaaa89)=alpha*xyzzyaaae89
enddo
elseif(incy==1)then
do xyzzyaaaa89=1,m
xyzzyaaad89=0.d0
xyzzyaaae89=0.d0
do xyzzyaaab89=1,inums(xyzzyaaaa89)
xyzzyaaac89=ias(xyzzyaaab89,xyzzyaaaa89)
xyzzyaaad89=xyzzyaaad89+a(xyzzyaaaa89,xyzzyaaac89)*xr(1+(xyzzyaaac89-1&
&)*incx)
xyzzyaaae89=xyzzyaaae89+a(xyzzyaaaa89,xyzzyaaac89)*xi(1+(xyzzyaaac89-1&
&)*incx)
enddo
yr(xyzzyaaaa89)=alpha*xyzzyaaad89
yi(xyzzyaaaa89)=alpha*xyzzyaaae89
enddo
elseif(incx==1)then
do xyzzyaaaa89=1,m
xyzzyaaad89=0.d0
xyzzyaaae89=0.d0
do xyzzyaaab89=1,inums(xyzzyaaaa89)
xyzzyaaac89=ias(xyzzyaaab89,xyzzyaaaa89)
xyzzyaaad89=xyzzyaaad89+a(xyzzyaaaa89,xyzzyaaac89)*xr(xyzzyaaac89)
xyzzyaaae89=xyzzyaaae89+a(xyzzyaaaa89,xyzzyaaac89)*xi(xyzzyaaac89)
enddo
yr(1+(xyzzyaaaa89-1)*incy)=alpha*xyzzyaaad89
yi(1+(xyzzyaaaa89-1)*incy)=alpha*xyzzyaaae89
enddo
else
do xyzzyaaaa89=1,m
xyzzyaaad89=0.d0
xyzzyaaae89=0.d0
do xyzzyaaab89=1,inums(xyzzyaaaa89)
xyzzyaaac89=ias(xyzzyaaab89,xyzzyaaaa89)
xyzzyaaad89=xyzzyaaad89+a(xyzzyaaaa89,xyzzyaaac89)*xr(1+(xyzzyaaac89-1&
&)*incx)
xyzzyaaae89=xyzzyaaae89+a(xyzzyaaaa89,xyzzyaaac89)*xi(1+(xyzzyaaac89-1&
&)*incx)
enddo
yr(1+(xyzzyaaaa89-1)*incy)=alpha*xyzzyaaad89
yi(1+(xyzzyaaaa89-1)*incy)=alpha*xyzzyaaae89
enddo
endif
end subroutine d2_sparse_matrix_vector
subroutine d_sparse_vector_matrix(m,alpha,a,lda,ias,ldias,inums,x,incx&
&,y,incy)
implicit none
integer,intent(in) :: m,lda,ldias,ias(ldias,*),inums(lda),incx,incy
real(dp),intent(in) :: a(lda,*),x(*),alpha
real(dp),intent(inout) :: y(*)
integer xyzzyaaaa90,xyzzyaaab90,xyzzyaaac90
real(dp) xyzzyaaad90
if(incx==1.and.incy==1)then
do xyzzyaaab90=1,m
xyzzyaaad90=0.d0
do xyzzyaaaa90=1,inums(xyzzyaaab90)
xyzzyaaac90=ias(xyzzyaaaa90,xyzzyaaab90)
xyzzyaaad90=xyzzyaaad90+x(xyzzyaaac90)*a(xyzzyaaac90,xyzzyaaab90)
enddo
y(xyzzyaaab90)=alpha*xyzzyaaad90
enddo
elseif(incy==1)then
do xyzzyaaab90=1,m
xyzzyaaad90=0.d0
do xyzzyaaaa90=1,inums(xyzzyaaab90)
xyzzyaaac90=ias(xyzzyaaaa90,xyzzyaaab90)
xyzzyaaad90=xyzzyaaad90+x(1+(xyzzyaaac90-1)*incx)*a(xyzzyaaac90,xyzzya&
&aab90)
enddo
y(xyzzyaaab90)=alpha*xyzzyaaad90
enddo
elseif(incx==1)then
do xyzzyaaab90=1,m
xyzzyaaad90=0.d0
do xyzzyaaaa90=1,inums(xyzzyaaab90)
xyzzyaaac90=ias(xyzzyaaaa90,xyzzyaaab90)
xyzzyaaad90=xyzzyaaad90+x(xyzzyaaac90)*a(xyzzyaaac90,xyzzyaaab90)
enddo
y(1+(xyzzyaaab90-1)*incy)=alpha*xyzzyaaad90
enddo
else
do xyzzyaaab90=1,m
xyzzyaaad90=0.d0
do xyzzyaaaa90=1,inums(xyzzyaaab90)
xyzzyaaac90=ias(xyzzyaaaa90,xyzzyaaab90)
xyzzyaaad90=xyzzyaaad90+x(1+(xyzzyaaac90-1)*incx)*a(xyzzyaaac90,xyzzya&
&aab90)
enddo
y(1+(xyzzyaaab90-1)*incy)=alpha*xyzzyaaad90
enddo
endif
end subroutine d_sparse_vector_matrix
subroutine d2_sparse_vector_matrix(m,alpha,a,lda,ias,ldias,inums,xr,xi&
&,incx,yr,yi,incy)
implicit none
integer,intent(in) :: m,lda,ldias,ias(ldias,*),inums(lda),incx,incy
real(dp),intent(in) :: a(lda,*),xr(*),xi(*),alpha
real(dp),intent(inout) :: yr(*),yi(*)
integer xyzzyaaaa91,xyzzyaaab91,xyzzyaaac91
real(dp) xyzzyaaad91,xyzzyaaae91
if(incx==1.and.incy==1)then
do xyzzyaaab91=1,m
xyzzyaaad91=0.d0
xyzzyaaae91=0.d0
do xyzzyaaaa91=1,inums(xyzzyaaab91)
xyzzyaaac91=ias(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaad91=xyzzyaaad91+xr(xyzzyaaac91)*a(xyzzyaaac91,xyzzyaaab91)
xyzzyaaae91=xyzzyaaae91+xi(xyzzyaaac91)*a(xyzzyaaac91,xyzzyaaab91)
enddo
yr(xyzzyaaab91)=alpha*xyzzyaaad91
yi(xyzzyaaab91)=alpha*xyzzyaaae91
enddo
elseif(incy==1)then
do xyzzyaaab91=1,m
xyzzyaaad91=0.d0
xyzzyaaae91=0.d0
do xyzzyaaaa91=1,inums(xyzzyaaab91)
xyzzyaaac91=ias(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaad91=xyzzyaaad91+xr(1+(xyzzyaaac91-1)*incx)*a(xyzzyaaac91,xyzzy&
&aaab91)
xyzzyaaae91=xyzzyaaae91+xi(1+(xyzzyaaac91-1)*incx)*a(xyzzyaaac91,xyzzy&
&aaab91)
enddo
yr(xyzzyaaab91)=alpha*xyzzyaaad91
yi(xyzzyaaab91)=alpha*xyzzyaaae91
enddo
elseif(incx==1)then
do xyzzyaaab91=1,m
xyzzyaaad91=0.d0
xyzzyaaae91=0.d0
do xyzzyaaaa91=1,inums(xyzzyaaab91)
xyzzyaaac91=ias(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaad91=xyzzyaaad91+xr(xyzzyaaac91)*a(xyzzyaaac91,xyzzyaaab91)
xyzzyaaae91=xyzzyaaae91+xi(xyzzyaaac91)*a(xyzzyaaac91,xyzzyaaab91)
enddo
yr(1+(xyzzyaaab91-1)*incy)=alpha*xyzzyaaad91
yi(1+(xyzzyaaab91-1)*incy)=alpha*xyzzyaaae91
enddo
else
do xyzzyaaab91=1,m
xyzzyaaad91=0.d0
xyzzyaaae91=0.d0
do xyzzyaaaa91=1,inums(xyzzyaaab91)
xyzzyaaac91=ias(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaad91=xyzzyaaad91+xr(1+(xyzzyaaac91-1)*incx)*a(xyzzyaaac91,xyzzy&
&aaab91)
xyzzyaaae91=xyzzyaaae91+xi(1+(xyzzyaaac91-1)*incx)*a(xyzzyaaac91,xyzzy&
&aaab91)
enddo
yr(1+(xyzzyaaab91-1)*incy)=alpha*xyzzyaaad91
yi(1+(xyzzyaaab91-1)*incy)=alpha*xyzzyaaae91
enddo
endif
end subroutine d2_sparse_vector_matrix
complex(dp) function zdot_rc(n,x,incx,y,incy)
implicit none
integer,intent(in) :: n,incx,incy
real(dp),intent(in) :: x(*)
complex(dp),intent(in) :: y(*)
integer xyzzyaaaa92,xyzzyaaab92,xyzzyaaac92
complex(dp) xyzzyaaad92
if(n<=0)then
zdot_rc=czero
elseif(incx==1.and.incy==1)then
zdot_rc=sum(x(1:n)*y(1:n))
else
if(incx>=0)then
xyzzyaaaa92=1
else
xyzzyaaaa92=(-n+1)*incx+1
endif
if(incy>=0)then
xyzzyaaab92=1
else
xyzzyaaab92=(-n+1)*incy+1
endif
xyzzyaaad92=czero
do xyzzyaaac92=1,n
xyzzyaaad92=xyzzyaaad92+x(xyzzyaaaa92)*y(xyzzyaaab92)
xyzzyaaaa92=xyzzyaaaa92+incx
xyzzyaaab92=xyzzyaaab92+incy
enddo
zdot_rc=xyzzyaaad92
endif
end function zdot_rc
subroutine dzscal(n,sa,x,incx)
implicit none
integer,intent(in) :: incx,n
real(dp),intent(in) :: sa
complex(dp),intent(inout) :: x(*)
integer xyzzyaaaa93,xyzzyaaab93,xyzzyaaac93
if(n<=0)then
return
elseif(incx==1)then
xyzzyaaac93=mod(n,5)
do xyzzyaaaa93=1,xyzzyaaac93
x(xyzzyaaaa93)=sa*x(xyzzyaaaa93)
enddo
do xyzzyaaaa93=xyzzyaaac93+1,n,5
x(xyzzyaaaa93)=sa*x(xyzzyaaaa93)
x(xyzzyaaaa93+1)=sa*x(xyzzyaaaa93+1)
x(xyzzyaaaa93+2)=sa*x(xyzzyaaaa93+2)
x(xyzzyaaaa93+3)=sa*x(xyzzyaaaa93+3)
x(xyzzyaaaa93+4)=sa*x(xyzzyaaaa93+4)
enddo
else
if(incx>=0)then
xyzzyaaab93=1
else
xyzzyaaab93=(-n+1)*incx+1
endif
do xyzzyaaaa93=1,n
x(xyzzyaaab93)=sa*x(xyzzyaaab93)
xyzzyaaab93=xyzzyaaab93+incx
enddo
endif
end subroutine dzscal
complex(dp) function zdotu_cc(n,x_re,x_im,incx,y_re,y_im,incy)
implicit none
integer,intent(in) :: n,incx,incy
real(dp),intent(in) :: x_re(*),y_re(*),x_im(*),y_im(*)
real(dp) xyzzyaaaa94,xyzzyaaab94
xyzzyaaaa94=ddot(n,x_re(1),incx,y_re(1),incy)-ddot(n,x_im(1),incx,y_im&
&(1),incy)
xyzzyaaab94=ddot(n,x_re(1),incx,y_im(1),incy)+ddot(n,x_im(1),incx,y_re&
&(1),incy)
zdotu_cc=cmplx(xyzzyaaaa94,xyzzyaaab94,dp)
end function
subroutine zaxpy_cc(n,a_r,a_i,x_r,x_i,incx,y_r,y_i,incy)
implicit none
integer,intent(in) :: n,incx,incy
real(dp),intent(in) :: a_r,a_i,x_r(*),x_i(*)
real(dp),intent(inout) :: y_r(*),y_i(*)
call daxpy(n,a_r,x_r(1),incx,y_r(1),incy)
call daxpy(n,a_i,x_r(1),incx,y_i(1),incy)
call daxpy(n,a_r,x_i(1),incx,y_i(1),incy)
call daxpy(n,-a_i,x_i(1),incx,y_r(1),incy)
end subroutine
subroutine zscal_cc(n,a_r,a_i,zx_r,zx_i,incx)
implicit none
integer,intent(in) :: n,incx
real(dp),intent(in) :: a_r,a_i
real(dp),intent(inout) :: zx_r(*),zx_i(*)
real(dp) xyzzyaaaa96,xyzzyaaab96
integer xyzzyaaac96
do xyzzyaaac96=1,n*incx,incx
xyzzyaaaa96=zx_r(xyzzyaaac96)
xyzzyaaab96=zx_i(xyzzyaaac96)
zx_r(xyzzyaaac96)=a_r*xyzzyaaaa96-a_i*xyzzyaaab96
zx_i(xyzzyaaac96)=a_r*xyzzyaaab96+a_i*xyzzyaaaa96
enddo
end subroutine
subroutine zgemm_cc(transa,transb,m,n,k,alpha,a_re,a_im,lda,b_re,b_im,&
&ldb,beta,c_re,c_im,ldc)
implicit none
character,intent(in) :: transa,transb
integer,intent(in) :: m,n,k,lda,ldb,ldc
real(dp),intent(in) :: alpha,a_re(lda,*),a_im(lda,*),b_re(ldb,*),b_im(&
&ldb,*),beta
real(dp),intent(inout) :: c_re(ldc,*),c_im(ldc,*)
call dgemm(transa,transb,m,n,k,alpha,a_re(1,1),lda,b_re(1,1),ldb,beta,&
&c_re(1,1),ldc)
call dgemm(transa,transb,m,n,k,alpha,a_im(1,1),lda,b_re(1,1),ldb,beta,&
&c_im(1,1),ldc)
call dgemm(transa,transb,m,n,k,alpha,a_re(1,1),lda,b_im(1,1),ldb,1.d0,&
&c_im(1,1),ldc)
call dgemm(transa,transb,m,n,k,-alpha,a_im(1,1),lda,b_im(1,1),ldb,1.d0&
&,c_re(1,1),ldc)
end subroutine
subroutine zgemv_cc(trans,m,n,alpha,a_r,a_i,lda,x_r,x_i,incx,beta,y_r,&
&y_i,incy)
implicit none
character,intent(in) :: trans
integer,intent(in) :: m,n,lda,incx,incy
real(dp),intent(in) :: alpha,beta,a_r(lda,*),a_i(lda,*),x_r(*),x_i(*)
real(dp),intent(inout) :: y_r(*),y_i(*)
call dgemv(trans,m,n,alpha,a_r(1,1),lda,x_r(1),incx,beta,y_r(1),incy)
call dgemv(trans,m,n,-alpha,a_i(1,1),lda,x_i(1),incx,1d0,y_r(1),incy)
call dgemv(trans,m,n,alpha,a_r(1,1),lda,x_i(1),incx,beta,y_i(1),incy)
call dgemv(trans,m,n,alpha,a_i(1,1),lda,x_r(1),incx,1d0,y_i(1),incy)
end subroutine
subroutine zger_cc(m,n,alpha,x_re,x_im,incx,y_re,y_im,incy,a_re,a_im,l&
&da)
implicit none
integer,intent(in) :: m,n,incx,incy,lda
real(dp),intent(in) :: x_re(*),x_im(*),y_re(*),y_im(*),alpha
real(dp),intent(inout) :: a_re(lda,*),a_im(lda,*)
call dger(m,n,alpha,x_re,incx,y_re,incy,a_re,lda)
call dger(m,n,alpha,x_im,incx,y_re,incy,a_im,lda)
call dger(m,n,alpha,x_re,incx,y_im,incy,a_im,lda)
call dger(m,n,-alpha,x_im,incx,y_im,incy,a_re,lda)
end subroutine
subroutine d_or_z_gemm(is_complex,transa,transb,m,n,k,alpha,a_re,a_im,&
&lda,b_re,b_im,ldb,beta,c_re,c_im,ldc)
implicit none
logical,intent(in) :: is_complex
character,intent(in) :: transa,transb
integer,intent(in) :: m,n,k,lda,ldb,ldc
real(dp),intent(in) :: alpha,a_re(lda,*),a_im(lda,*),b_re(ldb,*),b_im(&
&ldb,*),beta
real(dp),intent(inout) :: c_re(ldc,*),c_im(ldc,*)
if(is_complex)then
call zgemm_cc(transa,transb,m,n,k,alpha,a_re,a_im,lda,b_re,b_im,ldb,be&
&ta,c_re,c_im,ldc)
else
call dgemm(transa,transb,m,n,k,alpha,a_re,lda,b_re,ldb,beta,c_re,ldc)
endif
end subroutine
subroutine d_or_z_gemv(is_complex,trans,m,n,alpha,a_r,a_i,lda,x_r,x_i,&
&incx,beta,y_r,y_i,incy)
implicit none
logical,intent(in) :: is_complex
character,intent(in) :: trans
integer,intent(in) :: m,n,lda,incx,incy
real(dp),intent(in) :: alpha,beta,a_r(lda,*),a_i(lda,*),x_r(*),x_i(*)
real(dp),intent(inout) :: y_r(*),y_i(*)
if(is_complex)then
call zgemv_cc(trans,m,n,alpha,a_r,a_i,lda,x_r,x_i,incx,beta,y_r,y_i,in&
&cy)
else
call dgemv(trans,m,n,alpha,a_r,lda,x_r,incx,beta,y_r,incy)
endif
end subroutine
subroutine d_or_z_ger(is_complex,m,n,alpha,x_re,x_im,incx,y_re,y_im,in&
&cy,a_re,a_im,lda)
implicit none
logical,intent(in) :: is_complex
integer,intent(in) :: m,n,incx,incy,lda
real(dp),intent(in) :: x_re(*),x_im(*),y_re(*),y_im(*),alpha
real(dp),intent(inout) :: a_re(lda,*),a_im(lda,*)
if(is_complex)then
call zger_cc(m,n,alpha,x_re,x_im,incx,y_re,y_im,incy,a_re,a_im,lda)
else
call dger(m,n,alpha,x_re,incx,y_re,incy,a_re,lda)
endif
end subroutine
subroutine d_or_z_scal(is_complex,n,a_r,a_i,zx_r,zx_i,incx)
implicit none
integer,intent(in) :: n,incx
logical, intent(in) :: is_complex
real(dp),intent(in) :: a_r,a_i
real(dp),intent(inout) :: zx_r(*),zx_i(*)
if(is_complex)then
call zscal_cc(n,a_r,a_i,zx_r,zx_i,incx)
else
call dscal(n,a_r,zx_r,incx)
endif
end subroutine
real(dp) function param_limiting(smooth,param,h,hilim,l,lolim)
implicit none
real(dp),intent(in) :: param,hilim,lolim
logical,intent(in) :: smooth,h,l
param_limiting=param
if(.not.(h.or.l))return
if(smooth)then
param_limiting=xyzzyaaaj1(param,h,hilim,l,lolim)
else
if(h.and.param>hilim)param_limiting=hilim
if(l.and.param<lolim)param_limiting=lolim
endif
end function param_limiting
real(dp) function xyzzyaaaj1(param,h,hilim,l,lolim)
implicit none
real(dp),intent(in) :: param,hilim,lolim
logical,intent(in) :: h,l
real(dp) xyzzyaaaa105,xyzzyaaab105,xyzzyaaac105
if(h.and.l.and.hilim<=lolim)call errstop('LIMIT_PARAB','Called with lo&
&wer limit >= upper limit.')
if(h.and.l)then
xyzzyaaac105=param-lolim
xyzzyaaaa105=hilim+lolim
xyzzyaaab105=hilim-lolim
xyzzyaaaj1=.5d0*(xyzzyaaaa105-xyzzyaaab105*cos(pi*xyzzyaaac105/xyzzyaa&
&ab105))
elseif(h)then
xyzzyaaac105=abs(hilim-param)
xyzzyaaaj1=hilim-xyzzyaaac105*xyzzyaaac105/(xyzzyaaac105+1)
elseif(l)then
xyzzyaaac105=abs(param-lolim)
xyzzyaaaj1=lolim+xyzzyaaac105*xyzzyaaac105/(xyzzyaaac105+1)
else
xyzzyaaaj1=param
endif
end function xyzzyaaaj1
real(dp) function param_unlimiting(smooth,param,h,hilim,l,lolim)
implicit none
real(dp),intent(in) :: param,hilim,lolim
logical,intent(in) :: smooth,h,l
param_unlimiting=param
if(.not.(h.or.l))return
if(smooth)param_unlimiting=xyzzyaaak1(param,h,hilim,l,lolim)
end function param_unlimiting
real(dp) function xyzzyaaak1(param,h,hilim,l,lolim)
implicit none
real(dp),intent(in) :: param,hilim,lolim
logical,intent(in) :: h,l
real(dp) xyzzyaaaa107,xyzzyaaab107,xyzzyaaac107,xyzzyaaad107,xyzzyaaae&
&107,xyzzyaaaf107
if(h.and.l.and.hilim<=lolim)call errstop('UNLIMIT_PARAB','Called with &
&lower limit >= upper limit.')
if(h.and.l)then
xyzzyaaaa107=param
if(xyzzyaaaa107>hilim)xyzzyaaaa107=hilim
if(xyzzyaaaa107<lolim)xyzzyaaaa107=lolim
xyzzyaaac107=lolim+hilim
xyzzyaaad107=hilim-lolim
xyzzyaaae107=(xyzzyaaac107-xyzzyaaaa107-xyzzyaaaa107)/xyzzyaaad107
xyzzyaaae107=min(1.d0,xyzzyaaae107)
xyzzyaaae107=max(-1.d0,xyzzyaaae107)
xyzzyaaak1=lolim+xyzzyaaad107*one_over_pi*acos(xyzzyaaae107)
elseif(h)then
xyzzyaaaa107=param
if(xyzzyaaaa107>hilim)xyzzyaaaa107=hilim
xyzzyaaab107=abs(hilim-xyzzyaaaa107)
xyzzyaaaf107=xyzzyaaab107*(xyzzyaaab107+4.d0)
xyzzyaaaf107=max(0.d0,xyzzyaaaf107)
xyzzyaaak1=hilim-.5d0*(xyzzyaaab107+sqrt(xyzzyaaaf107))
elseif(l)then
xyzzyaaaa107=param
if(xyzzyaaaa107<lolim)xyzzyaaaa107=lolim
xyzzyaaab107=abs(xyzzyaaaa107-lolim)
xyzzyaaaf107=xyzzyaaab107*(xyzzyaaab107+4.d0)
xyzzyaaaf107=max(0.d0,xyzzyaaaf107)
xyzzyaaak1=lolim+.5d0*(xyzzyaaab107+sqrt(xyzzyaaaf107))
else
xyzzyaaak1=param
endif
end function xyzzyaaak1
integer function choose(a,b)
implicit none
integer,intent(in) :: a,b
integer xyzzyaaaa108
integer(i64) xyzzyaaab108,xyzzyaaac108
if(a<0)call errstop('CHOOSE','a<0')
if(b<0)call errstop('CHOOSE','b<0.')
if(b>a)call errstop('CHOOSE','b>a.')
xyzzyaaab108=1_i64
xyzzyaaac108=1_i64
do xyzzyaaaa108=1,b
xyzzyaaac108=xyzzyaaac108*int(a+1-xyzzyaaaa108,i64)
xyzzyaaab108=xyzzyaaab108*int(xyzzyaaaa108,i64)
enddo
choose=int(xyzzyaaac108/xyzzyaaab108)
end function choose
subroutine next_permutation(n,iv,even,finished)
implicit none
integer,intent(in) :: n
integer,intent(inout) :: iv(n)
logical,intent(out) :: finished
logical,intent(inout) :: even
integer xyzzyaaaa109,xyzzyaaab109,xyzzyaaac109,xyzzyaaad109
finished=n<2
do xyzzyaaaa109=n-1,1,-1
if(iv(xyzzyaaaa109)<iv(xyzzyaaaa109+1))exit
if(xyzzyaaaa109==1)finished=.true.
enddo
if(.not.finished)then
do xyzzyaaab109=n,xyzzyaaaa109+1,-1
if(iv(xyzzyaaab109)>iv(xyzzyaaaa109))exit
enddo
call swap1(iv(xyzzyaaab109),iv(xyzzyaaaa109))
even=.not.even
else
xyzzyaaaa109=0
endif
xyzzyaaad109=n
do xyzzyaaac109=xyzzyaaaa109+1,(xyzzyaaaa109+n)/2
call swap1(iv(xyzzyaaac109),iv(xyzzyaaad109))
even=.not.even
xyzzyaaad109=xyzzyaaad109-1
enddo
end subroutine next_permutation
subroutine correlation_time(odata,otau,otau_err)
implicit none
real(dp),intent(in) :: odata(:)
real(dp),intent(out) :: otau,otau_err
integer xyzzyaaaa110,xyzzyaaab110,xyzzyaaac110,xyzzyaaad110,xyzzyaaae1&
&10
real(dp) xyzzyaaaf110,xyzzyaaag110,xyzzyaaah110,xyzzyaaai110,xyzzyaaaj&
&110,xyzzyaaak110,xyzzyaaal110,xyzzyaaam110,xyzzyaaan110
real(dp),parameter :: xyzzyaaao110=1.d-100
otau=-1.d0
otau_err=-1.d0
xyzzyaaaa110=size(odata,1)
if(xyzzyaaaa110<10)return
xyzzyaaan110=1.d0/real(xyzzyaaaa110,dp)
xyzzyaaaf110=sum(odata)*xyzzyaaan110
xyzzyaaag110=xyzzyaaaf110**2
xyzzyaaah110=sum(odata**2)*xyzzyaaan110
xyzzyaaai110=xyzzyaaah110-xyzzyaaag110
if(xyzzyaaai110<xyzzyaaao110)return
xyzzyaaaj110=1.d0/xyzzyaaai110
otau=1.d0
xyzzyaaad110=xyzzyaaaa110
do xyzzyaaab110=1,xyzzyaaaa110-1
xyzzyaaad110=xyzzyaaad110-1
xyzzyaaal110=0.d0
xyzzyaaae110=xyzzyaaab110
do xyzzyaaac110=1,xyzzyaaad110
xyzzyaaae110=xyzzyaaae110+1
xyzzyaaal110=xyzzyaaal110+(odata(xyzzyaaac110)-xyzzyaaaf110)*(odata(xy&
&zzyaaae110)-xyzzyaaaf110)
enddo
xyzzyaaak110=xyzzyaaal110*xyzzyaaaj110/real(xyzzyaaad110,dp)
otau=otau+2.d0*xyzzyaaak110
if(xyzzyaaab110>=nint(3.d0*otau).or.xyzzyaaab110==xyzzyaaaa110-1)then
xyzzyaaam110=real(xyzzyaaab110,dp)
exit
endif
enddo
otau_err=otau*sqrt((4.d0*xyzzyaaam110+2.d0)*xyzzyaaan110)
end subroutine correlation_time
subroutine correlation_time_alt(odata1,odata2,otau,otau_err)
implicit none
real(dp),intent(in) :: odata1(:),odata2(:)
real(dp),intent(out) :: otau,otau_err
integer xyzzyaaaa111,xyzzyaaab111,xyzzyaaac111,xyzzyaaad111,xyzzyaaae1&
&11
real(dp) xyzzyaaaf111,xyzzyaaag111,xyzzyaaah111,xyzzyaaai111,xyzzyaaaj&
&111,xyzzyaaak111,xyzzyaaal111
otau=-1.d0
otau_err=-1.d0
xyzzyaaaa111=size(odata1,1)
if(xyzzyaaaa111<10)return
xyzzyaaal111=1.d0/real(xyzzyaaaa111,dp)
xyzzyaaaf111=sum(odata1)/sum(odata2)
xyzzyaaag111=sum((odata1-odata2*xyzzyaaaf111)**2)*xyzzyaaal111
xyzzyaaah111=1.d0/xyzzyaaag111
otau=1.d0
xyzzyaaad111=xyzzyaaaa111
do xyzzyaaab111=1,xyzzyaaaa111-1
xyzzyaaad111=xyzzyaaad111-1
xyzzyaaaj111=0.d0
xyzzyaaae111=xyzzyaaab111
do xyzzyaaac111=1,xyzzyaaad111
xyzzyaaae111=xyzzyaaae111+1
xyzzyaaaj111=xyzzyaaaj111+(odata1(xyzzyaaac111)-odata2(xyzzyaaac111)*x&
&yzzyaaaf111)*(odata1(xyzzyaaae111)-odata2(xyzzyaaae111)*xyzzyaaaf111)
enddo
xyzzyaaai111=xyzzyaaaj111*xyzzyaaah111/real(xyzzyaaad111,dp)
otau=otau+2.d0*xyzzyaaai111
if(xyzzyaaab111>=nint(3.d0*otau).or.xyzzyaaab111==xyzzyaaaa111-1)then
xyzzyaaak111=real(xyzzyaaab111,dp)
exit
endif
enddo
otau_err=otau*sqrt((4.d0*xyzzyaaak111+2.d0)*xyzzyaaal111)
end subroutine correlation_time_alt
subroutine reaverage(n,nmove,vector,dvector,prefactor,mean,std_dev)
implicit none
integer,intent(in) :: n,nmove
real(dp),intent(in) :: vector(:),dvector(:),prefactor
real(dp),intent(out) :: mean,std_dev
real(dp) xyzzyaaaa112,xyzzyaaab112,xyzzyaaac112,xyzzyaaad112,xyzzyaaae&
&112,xyzzyaaaf112
mean=0.d0
std_dev=0.d0
if(n<1)then
return
elseif(n==1)then
mean=vector(1)
std_dev=dvector(1)*prefactor
else
xyzzyaaad112=1.d0/dble(n)
xyzzyaaae112=1.d0/dble(n*nmove-1)
xyzzyaaaa112=sum(vector(1:n))
xyzzyaaab112=sum(vector(1:n)**2)
xyzzyaaac112=sum(dvector(1:n)**2)
mean=xyzzyaaaa112*xyzzyaaad112
xyzzyaaaf112=xyzzyaaab112+dble(nmove-1)*xyzzyaaac112
xyzzyaaaf112=xyzzyaaaf112*xyzzyaaad112-mean*mean
xyzzyaaaf112=xyzzyaaaf112*xyzzyaaae112
if(xyzzyaaaf112<0.d0)xyzzyaaaf112=0.d0
std_dev=prefactor*sqrt(xyzzyaaaf112)
endif
end subroutine reaverage
real(dp) function xyzzyaaal1(p,n,x)
implicit none
integer,intent(in) :: p,n
real(dp),intent(in) :: x(n)
integer xyzzyaaaa113,xyzzyaaab113,xyzzyaaac113,xyzzyaaad113,xyzzyaaae1&
&13,xyzzyaaaf113,xyzzyaaag113(n)
real(dp) xyzzyaaah113
do xyzzyaaaa113=1,n
xyzzyaaag113(xyzzyaaaa113)=xyzzyaaaa113
enddo
xyzzyaaae113=1
xyzzyaaad113=n
do
if(xyzzyaaad113<=xyzzyaaae113+1)then
if(xyzzyaaad113==xyzzyaaae113+1)then
if(x(xyzzyaaag113(xyzzyaaad113))<x(xyzzyaaag113(xyzzyaaae113)))call sw&
&ap1(xyzzyaaag113(xyzzyaaad113),xyzzyaaag113(xyzzyaaae113))
endif
xyzzyaaal1=x(xyzzyaaag113(p))
return
endif
xyzzyaaaf113=(xyzzyaaae113+xyzzyaaad113)/2
call swap1(xyzzyaaag113(xyzzyaaaf113),xyzzyaaag113(xyzzyaaae113+1))
if(x(xyzzyaaag113(xyzzyaaae113))>x(xyzzyaaag113(xyzzyaaad113)))call sw&
&ap1(xyzzyaaag113(xyzzyaaae113),xyzzyaaag113(xyzzyaaad113))
if(x(xyzzyaaag113(xyzzyaaae113+1))>x(xyzzyaaag113(xyzzyaaad113)))call &
&swap1(xyzzyaaag113(xyzzyaaae113+1),xyzzyaaag113(xyzzyaaad113))
if(x(xyzzyaaag113(xyzzyaaae113))>x(xyzzyaaag113(xyzzyaaae113+1)))call &
&swap1(xyzzyaaag113(xyzzyaaae113),xyzzyaaag113(xyzzyaaae113+1))
xyzzyaaab113=xyzzyaaae113+1
xyzzyaaac113=xyzzyaaad113
xyzzyaaah113=x(xyzzyaaag113(xyzzyaaae113+1))
do
do while(xyzzyaaab113<n)
xyzzyaaab113=xyzzyaaab113+1
if(x(xyzzyaaag113(xyzzyaaab113))>=xyzzyaaah113)exit
enddo
do while(xyzzyaaac113>1)
xyzzyaaac113=xyzzyaaac113-1
if(x(xyzzyaaag113(xyzzyaaac113))<=xyzzyaaah113)exit
enddo
if(xyzzyaaac113<xyzzyaaab113)exit
call swap1(xyzzyaaag113(xyzzyaaab113),xyzzyaaag113(xyzzyaaac113))
enddo
call swap1(xyzzyaaag113(xyzzyaaae113+1),xyzzyaaag113(xyzzyaaac113))
if(xyzzyaaac113>=p)xyzzyaaad113=xyzzyaaac113-1
if(xyzzyaaac113<=p)xyzzyaaae113=xyzzyaaab113
enddo
end function xyzzyaaal1
real(dp) function median(n,x)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: x(n)
integer xyzzyaaaa114
real(dp) xyzzyaaab114,xyzzyaaac114
if(mod(n,2)==0)then
xyzzyaaaa114=n/2
xyzzyaaab114=xyzzyaaal1(xyzzyaaaa114,n,x)
xyzzyaaac114=xyzzyaaal1(xyzzyaaaa114+1,n,x)
median=0.5d0*(xyzzyaaab114+xyzzyaaac114)
else
xyzzyaaaa114=(n+1)/2
median=xyzzyaaal1(xyzzyaaaa114,n,x)
endif
end function median
subroutine reduced_echelon(m,n,a)
implicit none
integer,intent(in) :: m,n
real(dp),intent(inout) :: a(m,n)
integer xyzzyaaaa115,xyzzyaaab115,xyzzyaaac115,xyzzyaaad115,xyzzyaaae1&
&15
real(dp) xyzzyaaaf115
real(dp),parameter :: xyzzyaaag115=1.d-10
xyzzyaaaa115=1
do xyzzyaaab115=1,n
if(m-xyzzyaaaa115>=1)then
xyzzyaaac115=xyzzyaaaa115
xyzzyaaaf115=abs(a(xyzzyaaaa115,xyzzyaaab115))
do xyzzyaaad115=xyzzyaaaa115+1,m
if(abs(a(xyzzyaaad115,xyzzyaaab115))>xyzzyaaaf115)then
xyzzyaaaf115=abs(a(xyzzyaaad115,xyzzyaaab115))
xyzzyaaac115=xyzzyaaad115
endif
enddo
if(xyzzyaaac115/=xyzzyaaaa115)then
do xyzzyaaae115=1,n
xyzzyaaaf115=a(xyzzyaaaa115,xyzzyaaae115)
a(xyzzyaaaa115,xyzzyaaae115)=a(xyzzyaaac115,xyzzyaaae115)
a(xyzzyaaac115,xyzzyaaae115)=xyzzyaaaf115
enddo
endif
endif
if(abs(a(xyzzyaaaa115,xyzzyaaab115))>xyzzyaaag115)then
call dscal(n,1.d0/a(xyzzyaaaa115,xyzzyaaab115),a(xyzzyaaaa115,1),m)
do xyzzyaaad115=1,xyzzyaaaa115-1
call daxpy(n,-a(xyzzyaaad115,xyzzyaaab115),a(xyzzyaaaa115,1),m,a(xyzzy&
&aaad115,1),m)
enddo
do xyzzyaaad115=xyzzyaaaa115+1,m
call daxpy(n,-a(xyzzyaaad115,xyzzyaaab115),a(xyzzyaaaa115,1),m,a(xyzzy&
&aaad115,1),m)
enddo
a(1:m,xyzzyaaab115)=0.d0
a(xyzzyaaaa115,xyzzyaaab115)=1.d0
xyzzyaaaa115=xyzzyaaaa115+1
else
a(xyzzyaaaa115:m,xyzzyaaab115)=0.d0
endif
if(xyzzyaaaa115==m+1)exit
enddo
end subroutine reduced_echelon
subroutine gauss_quadrature(x1,x2,n,func,integral)
integer,intent(in) :: n
real(dp),intent(in) :: x1,x2
real(dp),intent(out) :: integral
integer xyzzyaaaa116,xyzzyaaab116
real(dp),allocatable :: x(:),xyzzyaaac116(:)
interface
real(kind(1.d0)) function func(x)
real(kind(1.d0)),intent(in) :: x
end function func
end interface
allocate(x(n),xyzzyaaac116(n),stat=xyzzyaaab116)
call check_alloc(xyzzyaaab116,'GAUSS_QUADRATURE','')
call xyzzyaaad116
integral=0.d0
do xyzzyaaaa116=1,n
integral=integral+xyzzyaaac116(xyzzyaaaa116)*func(x(xyzzyaaaa116))
enddo
contains
subroutine xyzzyaaad116
implicit none
integer xyzzyaaaa118,xyzzyaaab118,xyzzyaaac118
integer,parameter :: xyzzyaaad118=10
real(dp) xyzzyaaae118,xyzzyaaaf118
real(dp),parameter :: xyzzyaaag118=1.d-14
real(dp),dimension((size(x)+1)/2) :: xyzzyaaah118,xyzzyaaai118,xyzzyaa&
&aj118,xyzzyaaak118,xyzzyaaal118,xyzzyaaam118
logical,dimension((size(x)+1)/2) :: xyzzyaaan118
if(size(x)/=size(xyzzyaaac116))call errstop('GAUSS_LEGENDRE_XW','Numbe&
&r of abscissa not equal to number of weights.')
xyzzyaaac118=(n+1)/2
xyzzyaaaf118=0.5d0*(x2+x1)
xyzzyaaae118=0.5d0*(x2-x1)
xyzzyaaal118=cos(pi*(xyzzyaaae116(1,1,xyzzyaaac118)-0.25d0)/(n+0.5d0))
xyzzyaaan118=.true.
do xyzzyaaaa118=1,xyzzyaaad118
where(xyzzyaaan118)
xyzzyaaah118=1.d0
xyzzyaaai118=0.d0
end where
do xyzzyaaab118=1,n
where(xyzzyaaan118)
xyzzyaaaj118=xyzzyaaai118
xyzzyaaai118=xyzzyaaah118
xyzzyaaah118=((2.d0*xyzzyaaab118-1.d0)*xyzzyaaal118*xyzzyaaai118-(xyzz&
&yaaab118-1.d0)*xyzzyaaaj118)/xyzzyaaab118
end where
enddo
where(xyzzyaaan118)
xyzzyaaak118=n*(xyzzyaaal118*xyzzyaaah118-xyzzyaaai118)/(xyzzyaaal118*&
&xyzzyaaal118-1.d0)
xyzzyaaam118=xyzzyaaal118
xyzzyaaal118=xyzzyaaam118-xyzzyaaah118/xyzzyaaak118
xyzzyaaan118=(abs(xyzzyaaal118-xyzzyaaam118)>xyzzyaaag118)
end where
if(.not.any(xyzzyaaan118))exit
enddo
if(xyzzyaaaa118==xyzzyaaad118+1)call errstop('GAUSS_LEGENDRE_XW','Too &
&many iterations.')
x(1:xyzzyaaac118)=xyzzyaaaf118-xyzzyaaae118*xyzzyaaal118
x(n:n-xyzzyaaac118+1:-1)=xyzzyaaaf118+xyzzyaaae118*xyzzyaaal118
xyzzyaaac116(1:xyzzyaaac118)=2.d0*xyzzyaaae118/((1.d0-xyzzyaaal118**2)&
&*xyzzyaaak118**2)
xyzzyaaac116(n:n-xyzzyaaac118+1:-1)=xyzzyaaac116(1:xyzzyaaac118)
end subroutine xyzzyaaad116
function xyzzyaaae116(first,increment,n)
integer,intent(in) :: first,increment,n
integer xyzzyaaaa119,xyzzyaaab119,xyzzyaaac119,xyzzyaaae116(n)
integer,parameter :: xyzzyaaad119=16,xyzzyaaae119=8
if(n>0)xyzzyaaae116(1)=first
if(n<=xyzzyaaad119)then
do xyzzyaaaa119=2,n
xyzzyaaae116(xyzzyaaaa119)=xyzzyaaae116(xyzzyaaaa119-1)+increment
enddo
else
do xyzzyaaaa119=2,xyzzyaaae119
xyzzyaaae116(xyzzyaaaa119)=xyzzyaaae116(xyzzyaaaa119-1)+increment
enddo
xyzzyaaac119=increment*xyzzyaaae119
xyzzyaaaa119=xyzzyaaae119
do
if(xyzzyaaaa119>=n)exit
xyzzyaaab119=xyzzyaaaa119+xyzzyaaaa119
xyzzyaaae116(xyzzyaaaa119+1:min(xyzzyaaab119,n))=xyzzyaaac119+xyzzyaaa&
&e116(1:min(xyzzyaaaa119,n-xyzzyaaaa119))
xyzzyaaac119=xyzzyaaac119+xyzzyaaac119
xyzzyaaaa119=xyzzyaaab119
enddo
endif
end function xyzzyaaae116
end subroutine gauss_quadrature
real(dp) function bessel_j0(x)
implicit none
real(dp),intent(in) :: x
real(dp) xyzzyaaaa120,xyzzyaaab120,xyzzyaaac120,xyzzyaaad120
real(dp),parameter :: xyzzyaaae120=1.d0,xyzzyaaaf120=-0.1098628627d-2,&
&xyzzyaaag120=0.2734510407d-4,xyzzyaaah120=-0.2073370639d-5,xyzzyaaai1&
&20=0.2093887211d-6
real(dp),parameter :: xyzzyaaaj120=-0.1562499995d-1,xyzzyaaak120=0.143&
&0488765d-3,xyzzyaaal120=-0.6911147651d-5,xyzzyaaam120=0.7621095161d-6&
&,xyzzyaaan120=-0.934945152d-7
real(dp),parameter :: xyzzyaaao120=57568490574.d0,xyzzyaaap120=-133625&
&90354.d0,xyzzyaaaq120=651619640.7d0,xyzzyaaar120=-11214424.18d0,xyzzy&
&aaas120=77392.33017d0,xyzzyaaat120=-184.9052456d0
real(dp),parameter :: xyzzyaaau120=57568490411.d0,xyzzyaaav120=1029532&
&985.d0,xyzzyaaaw120=9494680.718d0,xyzzyaaax120=59272.64853d0,xyzzyaaa&
&y120=267.8532712d0,xyzzyaaaz120=1.d0
xyzzyaaaa120=abs(x)
if(xyzzyaaaa120<8.d0)then
xyzzyaaad120=x*x
bessel_j0=(xyzzyaaao120+xyzzyaaad120*(xyzzyaaap120+xyzzyaaad120*(xyzzy&
&aaaq120+xyzzyaaad120*(xyzzyaaar120+xyzzyaaad120*(xyzzyaaas120+xyzzyaa&
&ad120*xyzzyaaat120)))))/(xyzzyaaau120+xyzzyaaad120*(xyzzyaaav120+xyzz&
&yaaad120*(xyzzyaaaw120+xyzzyaaad120*(xyzzyaaax120+xyzzyaaad120*(xyzzy&
&aaay120+xyzzyaaad120*xyzzyaaaz120)))))
else
xyzzyaaac120=8.d0/xyzzyaaaa120
xyzzyaaad120=xyzzyaaac120*xyzzyaaac120
xyzzyaaab120=xyzzyaaaa120-0.78539816339744830961d0
bessel_j0=sqrt(0.0795774715d0*xyzzyaaac120)*(cos(xyzzyaaab120)*(xyzzya&
&aae120+xyzzyaaad120*(xyzzyaaaf120+xyzzyaaad120*(xyzzyaaag120+xyzzyaaa&
&d120*(xyzzyaaah120+xyzzyaaad120*xyzzyaaai120))))-xyzzyaaac120*sin(xyz&
&zyaaab120)*(xyzzyaaaj120+xyzzyaaad120*(xyzzyaaak120+xyzzyaaad120*(xyz&
&zyaaal120+xyzzyaaad120*(xyzzyaaam120+xyzzyaaad120*xyzzyaaan120)))))
endif
end function bessel_j0
real(dp) function bessel_j1(x)
implicit none
real(dp),intent(in) :: x
real(dp) xyzzyaaaa121,xyzzyaaab121,xyzzyaaac121,xyzzyaaad121
real(dp),parameter :: xyzzyaaae121=72362614232.d0,xyzzyaaaf121=-789505&
&9235.d0,xyzzyaaag121=242396853.1d0,xyzzyaaah121=-2972611.439d0,xyzzya&
&aai121=15704.48260d0,xyzzyaaaj121=-30.16036606d0
real(dp),parameter :: xyzzyaaak121=144725228442.d0,xyzzyaaal121=230053&
&5178.d0,xyzzyaaam121=18583304.74d0,xyzzyaaan121=99447.43394d0,xyzzyaa&
&ao121=376.9991397d0
real(dp),parameter :: xyzzyaaap121=0.183105d-2,xyzzyaaaq121=-0.3516396&
&496d-4,xyzzyaaar121=0.2457520174d-5,xyzzyaaas121=-0.240337019d-6
real(dp),parameter :: xyzzyaaat121=0.04687499995d0,xyzzyaaau121=-0.200&
&2690873d-3,xyzzyaaav121=0.8449199096d-5,xyzzyaaaw121=-0.88228987d-6,x&
&yzzyaaax121=0.105787412d-6
xyzzyaaad121=abs(x)
if(xyzzyaaad121<8.d0)then
xyzzyaaaa121=x*x
bessel_j1=x*(xyzzyaaae121+xyzzyaaaa121*(xyzzyaaaf121+xyzzyaaaa121*(xyz&
&zyaaag121+xyzzyaaaa121*(xyzzyaaah121+xyzzyaaaa121*(xyzzyaaai121+xyzzy&
&aaaa121*xyzzyaaaj121))))) /(xyzzyaaak121+xyzzyaaaa121*(xyzzyaaal121+x&
&yzzyaaaa121*(xyzzyaaam121+xyzzyaaaa121*(xyzzyaaan121+xyzzyaaaa121*(xy&
&zzyaaao121+xyzzyaaaa121)))))
else
xyzzyaaab121=8.d0/xyzzyaaad121
xyzzyaaaa121=xyzzyaaab121*xyzzyaaab121
xyzzyaaac121=xyzzyaaad121-2.35619449019234492883d0
if(x>=0.d0)then
bessel_j1=sqrt(0.0795774715d0*xyzzyaaab121) *(cos(xyzzyaaac121)*(1.d0+&
&xyzzyaaaa121*(xyzzyaaap121+xyzzyaaaa121*(xyzzyaaaq121+xyzzyaaaa121*(x&
&yzzyaaar121+xyzzyaaaa121*xyzzyaaas121)))) -xyzzyaaab121*sin(xyzzyaaac&
&121)*(xyzzyaaat121+xyzzyaaaa121*(xyzzyaaau121+xyzzyaaaa121*(xyzzyaaav&
&121+xyzzyaaaa121*(xyzzyaaaw121+xyzzyaaaa121*xyzzyaaax121)))))
else
bessel_j1=-sqrt(0.0795774715d0*xyzzyaaab121) *(cos(xyzzyaaac121)*(1.d0&
&+xyzzyaaaa121*(xyzzyaaap121+xyzzyaaaa121*(xyzzyaaaq121+xyzzyaaaa121*(&
&xyzzyaaar121+xyzzyaaaa121*xyzzyaaas121)))) -xyzzyaaab121*sin(xyzzyaaa&
&c121)*(xyzzyaaat121+xyzzyaaaa121*(xyzzyaaau121+xyzzyaaaa121*(xyzzyaaa&
&v121+xyzzyaaaa121*(xyzzyaaaw121+xyzzyaaaa121*xyzzyaaax121)))))
endif
endif
end function bessel_j1
subroutine make_representable(x)
implicit none
real(dp),intent(inout) :: x
x=x
end subroutine make_representable
recursive subroutine bracket(xyzzyaaab123,xyzzyaaac123,f,xyzzyaaaa123,&
&xyzzyaaad123,xyzzyaaae123,xyzzyaaaf123,xyzzyaaag123)
implicit none
integer,intent(in) :: xyzzyaaaa123
real(dp),intent(in) :: xyzzyaaab123,xyzzyaaac123
real(dp),intent(out) :: xyzzyaaad123,xyzzyaaae123,xyzzyaaaf123
logical,intent(out) :: xyzzyaaag123
integer xyzzyaaah123
real(dp) xyzzyaaai123,x,xyzzyaaaj123,xyzzyaaak123,xyzzyaaal123
interface
function f(x)
use dsp,only : dp
real(dp),intent(in) :: x
real(dp) f
end function f
end interface
xyzzyaaai123=(xyzzyaaac123-xyzzyaaab123)/dble(xyzzyaaaa123-1)
xyzzyaaak123=huge(1.d0)
do xyzzyaaah123=2,xyzzyaaaa123-1
x=xyzzyaaab123+xyzzyaaai123*dble(xyzzyaaah123-1)
xyzzyaaaj123=f(x)
if(xyzzyaaaj123<xyzzyaaak123)then
xyzzyaaal123=x
xyzzyaaak123=xyzzyaaaj123
endif
enddo
xyzzyaaad123=xyzzyaaal123-xyzzyaaai123
xyzzyaaae123=xyzzyaaal123
xyzzyaaaf123=xyzzyaaal123+xyzzyaaai123
xyzzyaaag123=(f(xyzzyaaad123)>xyzzyaaak123.and.f(xyzzyaaaf123)>xyzzyaa&
&ak123)
end subroutine bracket
recursive subroutine gsbracket(xyzzyaaab125,xyzzyaaac125,xyzzyaaaa125,&
&f)
implicit none
interface
function f(x)
use dsp,only : dp
real(dp),intent(in) :: x
real(dp) f
end function f
end interface
real(dp),intent(out) :: xyzzyaaaa125
real(dp),intent(inout) :: xyzzyaaab125,xyzzyaaac125
real(dp) xyzzyaaad125,xyzzyaaae125,xyzzyaaaf125,xyzzyaaag125,xyzzyaaah&
&125,xyzzyaaai125,xyzzyaaaj125,x,xyzzyaaak125,xyzzyaaal125,xyzzyaaam12&
&5
real(dp),parameter :: xyzzyaaan125=1.6180339887498948482d0,xyzzyaaao12&
&5=100.d0,xyzzyaaap125=1000.d0*tiny(1.d0),xyzzyaaaq125=1.d-10
if(abs(xyzzyaaac125-xyzzyaaab125)<xyzzyaaaq125*abs(xyzzyaaab125))call &
&errstop('GSBRACKET','Initial guesses should be further apart.')
xyzzyaaad125=f(xyzzyaaab125)
xyzzyaaae125=f(xyzzyaaac125)
if(xyzzyaaae125>xyzzyaaad125)then
xyzzyaaag125=xyzzyaaab125
xyzzyaaab125=xyzzyaaac125
xyzzyaaac125=xyzzyaaag125
xyzzyaaag125=xyzzyaaad125
xyzzyaaad125=xyzzyaaae125
xyzzyaaae125=xyzzyaaag125
endif
xyzzyaaaa125=xyzzyaaac125+xyzzyaaan125*(xyzzyaaac125-xyzzyaaab125)
xyzzyaaaf125=f(xyzzyaaaa125)
do while(xyzzyaaae125>=xyzzyaaaf125)
xyzzyaaam125=xyzzyaaao125*(xyzzyaaaa125-xyzzyaaac125)
xyzzyaaaj125=(xyzzyaaac125-xyzzyaaab125)*(xyzzyaaae125-xyzzyaaaf125)
xyzzyaaai125=(xyzzyaaac125-xyzzyaaaa125)*(xyzzyaaae125-xyzzyaaad125)
xyzzyaaak125=xyzzyaaai125-xyzzyaaaj125
if(abs(xyzzyaaak125)>xyzzyaaap125)then
xyzzyaaal125=-0.5d0*((xyzzyaaac125-xyzzyaaaa125)*xyzzyaaai125-(xyzzyaa&
&ac125-xyzzyaaab125)*xyzzyaaaj125)/xyzzyaaak125
if(abs(xyzzyaaal125)>abs(xyzzyaaam125).or.(xyzzyaaal125>0.d0.neqv.xyzz&
&yaaam125>0.d0))xyzzyaaal125=xyzzyaaam125
else
xyzzyaaal125=xyzzyaaam125
endif
x=xyzzyaaac125+xyzzyaaal125
if(xyzzyaaac125>x.eqv.x>xyzzyaaaa125)then
xyzzyaaah125=f(x)
if(xyzzyaaah125<xyzzyaaaf125)then
xyzzyaaab125=xyzzyaaac125
xyzzyaaac125=x
exit
elseif(xyzzyaaah125>xyzzyaaae125)then
xyzzyaaaa125=x
exit
endif
elseif(x>=xyzzyaaaa125.eqv.xyzzyaaaa125>=xyzzyaaab125)then
xyzzyaaah125=f(x)
if(xyzzyaaaf125<xyzzyaaah125)then
xyzzyaaab125=xyzzyaaac125
xyzzyaaac125=xyzzyaaaa125
xyzzyaaaa125=x
exit
endif
xyzzyaaac125=xyzzyaaaa125
xyzzyaaae125=xyzzyaaaf125
xyzzyaaaa125=x
xyzzyaaaf125=xyzzyaaah125
endif
x=xyzzyaaaa125+xyzzyaaan125*(xyzzyaaaa125-xyzzyaaac125)
xyzzyaaab125=xyzzyaaac125
xyzzyaaac125=xyzzyaaaa125
xyzzyaaaa125=x
xyzzyaaad125=xyzzyaaae125
xyzzyaaae125=xyzzyaaaf125
xyzzyaaaf125=f(xyzzyaaaa125)
enddo
if(xyzzyaaaa125<xyzzyaaab125)then
xyzzyaaag125=xyzzyaaaa125
xyzzyaaaa125=xyzzyaaab125
xyzzyaaab125=xyzzyaaag125
xyzzyaaag125=xyzzyaaaf125
xyzzyaaaf125=xyzzyaaad125
xyzzyaaad125=xyzzyaaag125
endif
end subroutine gsbracket
recursive subroutine brentmin(xyzzyaaaa127,xyzzyaaab127,xyzzyaaac127,f&
&,x,xyzzyaaae127,xyzzyaaad127)
implicit none
real(dp),intent(in) :: xyzzyaaaa127,xyzzyaaab127,xyzzyaaac127
real(dp),intent(in),optional :: xyzzyaaad127
real(dp),intent(out) :: x,xyzzyaaae127
integer xyzzyaaaf127
integer,parameter :: xyzzyaaag127=100
real(dp) xyzzyaaah127,xyzzyaaai127,xyzzyaaaj127,xyzzyaaak127,xyzzyaaal&
&127,xyzzyaaam127,xyzzyaaan127,xyzzyaaao127,xyzzyaaap127,xyzzyaaaq127,&
&xyzzyaaar127,xyzzyaaas127,xyzzyaaat127,xyzzyaaau127,xyzzyaaav127,xyzz&
&yaaaw127,xyzzyaaax127
real(dp),save :: xyzzyaaay127
real(dp),parameter :: xyzzyaaaz127=0.3819660112501051518d0
logical,save :: xyzzyaaba127=.true.
interface
function f(x)
use dsp,only : dp
real(dp),intent(in) :: x
real(dp) f
end function f
end interface
if(xyzzyaaba127)then
xyzzyaaay127=sqrt(epsilon(1.d0))
xyzzyaaba127=.false.
endif
if(present(xyzzyaaad127))then
xyzzyaaap127=xyzzyaaad127
else
xyzzyaaap127=xyzzyaaay127
endif
xyzzyaaah127=xyzzyaaaa127
x=xyzzyaaab127
xyzzyaaai127=xyzzyaaac127
xyzzyaaae127=f(x)
xyzzyaaau127=x
xyzzyaaat127=x
xyzzyaaaw127=xyzzyaaae127
xyzzyaaax127=xyzzyaaae127
xyzzyaaak127=0.d0
do xyzzyaaaf127=1,xyzzyaaag127
xyzzyaaal127=0.5d0*(xyzzyaaah127+xyzzyaaai127)
xyzzyaaaq127=xyzzyaaay127*abs(x)+xyzzyaaap127*third
xyzzyaaar127=xyzzyaaaq127+xyzzyaaaq127
if(abs(x-xyzzyaaal127)<=xyzzyaaar127-0.5d0*(xyzzyaaai127-xyzzyaaah127)&
&)return
if(abs(xyzzyaaak127)>xyzzyaaaq127)then
xyzzyaaao127=(x-xyzzyaaau127)*(xyzzyaaae127-xyzzyaaaw127)
xyzzyaaan127=(x-xyzzyaaat127)*(xyzzyaaae127-xyzzyaaax127)
xyzzyaaam127=(x-xyzzyaaat127)*xyzzyaaan127-(x-xyzzyaaau127)*xyzzyaaao1&
&27
xyzzyaaan127=2.d0*(xyzzyaaan127-xyzzyaaao127)
if(xyzzyaaan127>0.d0)xyzzyaaam127=-xyzzyaaam127
xyzzyaaan127=abs(xyzzyaaan127)
if(abs(xyzzyaaam127)<abs(0.5d0*xyzzyaaan127*xyzzyaaak127).and.xyzzyaaa&
&m127>xyzzyaaan127*(xyzzyaaah127-x).and.xyzzyaaam127<xyzzyaaan127*(xyz&
&zyaaai127-x))then
xyzzyaaak127=xyzzyaaaj127
if(xyzzyaaas127-xyzzyaaah127<xyzzyaaar127.or.xyzzyaaai127-xyzzyaaas127&
&<xyzzyaaar127)then
xyzzyaaaj127=sign(xyzzyaaaq127,xyzzyaaal127-x)
else
xyzzyaaaj127=xyzzyaaam127/xyzzyaaan127
endif
else
if(x>=xyzzyaaal127)then
xyzzyaaak127=xyzzyaaah127-x
else
xyzzyaaak127=xyzzyaaai127-x
endif
xyzzyaaaj127=xyzzyaaaz127*xyzzyaaak127
endif
else
if(x>=xyzzyaaal127)then
xyzzyaaak127=xyzzyaaah127-x
else
xyzzyaaak127=xyzzyaaai127-x
endif
xyzzyaaaj127=xyzzyaaaz127*xyzzyaaak127
endif
if(abs(xyzzyaaaj127)>=xyzzyaaaq127)then
xyzzyaaas127=x+xyzzyaaaj127
else
xyzzyaaas127=x+sign(xyzzyaaaq127,xyzzyaaaj127)
endif
xyzzyaaav127=f(xyzzyaaas127)
if(xyzzyaaav127>xyzzyaaae127)then
if(xyzzyaaas127<x)then
xyzzyaaah127=xyzzyaaas127
else
xyzzyaaai127=xyzzyaaas127
endif
if(xyzzyaaav127<=xyzzyaaax127.or.xyzzyaaau127==x)then
xyzzyaaat127=xyzzyaaau127
xyzzyaaaw127=xyzzyaaax127
xyzzyaaau127=xyzzyaaas127
xyzzyaaax127=xyzzyaaav127
elseif(xyzzyaaav127<=xyzzyaaaw127.or.xyzzyaaat127==x.or.xyzzyaaat127==&
&xyzzyaaau127)then
xyzzyaaat127=xyzzyaaas127
xyzzyaaaw127=xyzzyaaav127
endif
else
if(xyzzyaaas127>=x)then
xyzzyaaah127=x
else
xyzzyaaai127=x
endif
xyzzyaaat127=xyzzyaaau127
xyzzyaaaw127=xyzzyaaax127
xyzzyaaau127=x
xyzzyaaax127=xyzzyaaae127
x=xyzzyaaas127
xyzzyaaae127=xyzzyaaav127
endif
enddo
call errwarn('BRENTMIN','Exceded maximum number of iterations.')
end subroutine brentmin
subroutine parabolic_min(x1,x2,x3,y1,y2,y3,x0,y0,rejected)
implicit none
real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
real(dp),intent(out) :: x0,y0
logical,intent(out) :: rejected
real(dp) xyzzyaaaa129,xyzzyaaab129,xyzzyaaac129,xyzzyaaad129,xyzzyaaae&
&129,xyzzyaaaf129,xyzzyaaag129,xyzzyaaah129,xyzzyaaai129,xyzzyaaaj129,&
&xyzzyaaak129,xyzzyaaal129,xyzzyaaam129,xyzzyaaan129,xyzzyaaao129,xyzz&
&yaaap129,xyzzyaaaq129
rejected=.false.
if(x1>=x2.or.x2>=x3.or.(y2>=y1.eqv.y3>=y2))then
rejected=.true.
return
endif
xyzzyaaah129=x1*x1
xyzzyaaai129=x2*x2
xyzzyaaaj129=x3*x3
xyzzyaaam129=x3-x1
xyzzyaaak129=(x2-x1)/xyzzyaaam129
xyzzyaaal129=(x3-x2)/xyzzyaaam129
xyzzyaaan129=y1*xyzzyaaal129
xyzzyaaao129=y2
xyzzyaaap129=y3*xyzzyaaak129
xyzzyaaag129=-xyzzyaaah129*xyzzyaaal129+xyzzyaaai129-xyzzyaaaj129*xyzz&
&yaaak129
xyzzyaaad129=-xyzzyaaan129+xyzzyaaao129-xyzzyaaap129
xyzzyaaae129=xyzzyaaan129*(x2+x3)-xyzzyaaao129*(x1+x3)+xyzzyaaap129*(x&
&1+x2)
xyzzyaaaf129=-xyzzyaaan129*x2*x3+xyzzyaaao129*x3*x1-xyzzyaaap129*x1*x2
xyzzyaaaq129=1.d0/xyzzyaaag129
xyzzyaaaa129=xyzzyaaad129*xyzzyaaaq129
xyzzyaaab129=xyzzyaaae129*xyzzyaaaq129
xyzzyaaac129=xyzzyaaaf129*xyzzyaaaq129
x0=-.5d0*xyzzyaaae129/xyzzyaaad129
y0=(xyzzyaaaa129*x0+xyzzyaaab129)*x0+xyzzyaaac129
if(x1==x0.or.x2==x0.or.x3==x0.or.y1==y0.or.y2==y0.or.y3==y0)then
rejected=.true.
return
endif
end subroutine parabolic_min
subroutine sph_bessel(r,n,r_c,bessel,dbessel,d2bessel)
use slaarnaag, only : pi,twopi
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,r_c
real(dp),intent(inout) :: bessel(n)
real(dp),intent(inout),optional :: dbessel(n),d2bessel(n)
integer xyzzyaaaa130
real(dp) xyzzyaaab130,xyzzyaaac130,xyzzyaaad130,xyzzyaaae130,xyzzyaaaf&
&130,xyzzyaaag130,xyzzyaaah130,xyzzyaaai130,xyzzyaaaj130
logical xyzzyaaak130,xyzzyaaal130
xyzzyaaak130=present(dbessel)
xyzzyaaal130=xyzzyaaak130.and.present(d2bessel)
xyzzyaaab130=1.d0/sqrt(twopi*r_c)
xyzzyaaac130=pi/r_c
xyzzyaaad130=0.d0
if(r/=0.d0)then
xyzzyaaah130=1.d0/r
if(xyzzyaaak130)xyzzyaaai130=xyzzyaaah130*xyzzyaaah130
if(xyzzyaaal130)xyzzyaaaj130=xyzzyaaai130*xyzzyaaah130
endif
do xyzzyaaaa130=1,n
xyzzyaaad130=xyzzyaaad130+xyzzyaaac130
xyzzyaaae130=xyzzyaaad130*r
xyzzyaaag130=sin(xyzzyaaae130)
if(r==0.d0)then
bessel(xyzzyaaaa130)=xyzzyaaab130*xyzzyaaad130
if(xyzzyaaak130)dbessel(xyzzyaaaa130)=0.d0
if(xyzzyaaal130)d2bessel(xyzzyaaaa130)=-2.d0*bessel(xyzzyaaaa130)*xyzz&
&yaaad130*xyzzyaaad130
else
bessel(xyzzyaaaa130)=xyzzyaaab130*xyzzyaaah130*xyzzyaaag130
if(xyzzyaaak130)then
xyzzyaaaf130=cos(xyzzyaaae130)
dbessel(xyzzyaaaa130)=xyzzyaaab130*xyzzyaaai130*(xyzzyaaae130*xyzzyaaa&
&f130-xyzzyaaag130)
if(xyzzyaaal130)d2bessel(xyzzyaaaa130)=xyzzyaaab130*xyzzyaaaj130*((2.d&
&0-xyzzyaaae130*xyzzyaaae130)*xyzzyaaag130-2.d0*xyzzyaaae130*xyzzyaaaf&
&130)
endif
endif
enddo
end subroutine sph_bessel
subroutine chebyshev(r,n,r_c,cheb,dcheb,d2cheb)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,r_c
real(dp),intent(inout) :: cheb(n)
real(dp),intent(inout),optional :: dcheb(n),d2cheb(n)
integer xyzzyaaaa131
real(dp) xyzzyaaab131
logical xyzzyaaac131,xyzzyaaad131
if(n==0)return
xyzzyaaab131=2.d0*r/r_c-1.d0
xyzzyaaac131=present(dcheb)
xyzzyaaad131=xyzzyaaac131.and.present(d2cheb)
cheb(1)=1.d0
if(xyzzyaaac131)dcheb(1)=0.d0
if(xyzzyaaad131)d2cheb(1)=0.d0
if(n>1)then
cheb(2)=xyzzyaaab131
if(xyzzyaaac131)dcheb(2)=1.d0
if(xyzzyaaad131)d2cheb(2)=0.d0
endif
do xyzzyaaaa131=3,n
cheb(xyzzyaaaa131)=2*xyzzyaaab131*cheb(xyzzyaaaa131-1)-cheb(xyzzyaaaa1&
&31-2)
if(xyzzyaaac131)dcheb(xyzzyaaaa131)=2.d0*(xyzzyaaab131*dcheb(xyzzyaaaa&
&131-1)+cheb(xyzzyaaaa131-1))-dcheb(xyzzyaaaa131-2)
if(xyzzyaaad131)d2cheb(xyzzyaaaa131)=2.d0*(xyzzyaaab131*d2cheb(xyzzyaa&
&aa131-1)+2.d0*dcheb(xyzzyaaaa131-1))-d2cheb(xyzzyaaaa131-2)
enddo
do xyzzyaaaa131=1,n
if(xyzzyaaac131)dcheb(xyzzyaaaa131)=2.d0/r_c*dcheb(xyzzyaaaa131)
if(xyzzyaaad131)d2cheb(xyzzyaaaa131)=4.d0/r_c**2*d2cheb(xyzzyaaaa131)
enddo
end subroutine chebyshev
subroutine exponential(r,n,k_0,exp_step,expn,dexpn,d2expn)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,k_0,exp_step
real(dp),intent(inout) :: expn(n)
real(dp),intent(inout),optional :: dexpn(n),d2expn(n)
integer xyzzyaaaa132,xyzzyaaab132
real(dp) xyzzyaaac132
logical xyzzyaaad132,xyzzyaaae132
xyzzyaaad132=present(dexpn)
xyzzyaaae132=xyzzyaaad132.and.present(d2expn)
xyzzyaaab132=-((n-1)/2)
xyzzyaaac132=-k_0*exp_step**xyzzyaaab132
do xyzzyaaaa132=1,n
expn(xyzzyaaaa132)=exp(xyzzyaaac132*r)
if(xyzzyaaad132)dexpn(xyzzyaaaa132)=xyzzyaaac132*expn(xyzzyaaaa132)
if(xyzzyaaae132)d2expn(xyzzyaaaa132)=xyzzyaaac132*dexpn(xyzzyaaaa132)
xyzzyaaac132=exp_step*xyzzyaaac132
enddo
end subroutine exponential
subroutine partial_rank_int(xdont,irngt,nord)
integer,intent(in)  :: nord,xdont(:)
integer,intent(out) :: irngt(:)
integer xyzzyaaaa133,xyzzyaaab133,xyzzyaaac133,xyzzyaaad133,xyzzyaaae1&
&33,xyzzyaaaf133,xyzzyaaag133,xyzzyaaah133,xyzzyaaai133,xyzzyaaaj133,x&
&yzzyaaak133,xyzzyaaal133,xyzzyaaam133,xyzzyaaan133,xyzzyaaao133,xyzzy&
&aaap133,xyzzyaaaq133,xyzzyaaar133,xyzzyaaas133,xyzzyaaat133,xyzzyaaau&
&133,xyzzyaaav133,xyzzyaaaw133,xyzzyaaax133,xyzzyaaay133,xyzzyaaaz133
integer,dimension(size(xdont)) :: xyzzyaaba133,xyzzyaabb133
xyzzyaaag133=size(xdont)
if(xyzzyaaag133<2)then
if(nord>=1)irngt(1)=1
return
endif
if(xdont(2)<xdont(1))then
xyzzyaaba133(1)=2
xyzzyaabb133(1)=1
else
xyzzyaaba133(1)=1
xyzzyaabb133(1)=2
endif
if(xyzzyaaag133<3)then
if(nord>=1)irngt(1)=xyzzyaabb133(1)
if(nord>=2)irngt(2)=xyzzyaaba133(1)
return
endif
if(xdont(3)>xdont(xyzzyaaba133(1)))then
xyzzyaaba133(2)=xyzzyaaba133(1)
if(xdont(3)>xdont(xyzzyaabb133(1)))then
xyzzyaaba133(1)=xyzzyaabb133(1)
xyzzyaabb133(1)=3
else
xyzzyaaba133(1)=3
endif
else
xyzzyaaba133(2)=3
endif
if(xyzzyaaag133<4)then
if(nord>=1)irngt(1)=xyzzyaabb133(1)
if(nord>=2)irngt(2)=xyzzyaaba133(1)
if(nord>=3)irngt(3)=xyzzyaaba133(2)
return
endif
if(xdont(xyzzyaaag133)>xdont(xyzzyaaba133(1)))then
xyzzyaaba133(3)=xyzzyaaba133(2)
xyzzyaaba133(2)=xyzzyaaba133(1)
if(xdont(xyzzyaaag133)>xdont(xyzzyaabb133(1)))then
xyzzyaaba133(1)=xyzzyaabb133(1)
xyzzyaabb133(1)=xyzzyaaag133
else
xyzzyaaba133(1)=xyzzyaaag133
endif
else
if(xdont(xyzzyaaag133)>xdont(xyzzyaaba133(2)))then
xyzzyaaba133(3)=xyzzyaaba133(2)
xyzzyaaba133(2)=xyzzyaaag133
else
xyzzyaaba133(3)=xyzzyaaag133
endif
endif
if(xyzzyaaag133<5)then
if(nord>=1)irngt(1)=xyzzyaabb133(1)
if(nord>=2)irngt(2)=xyzzyaaba133(1)
if(nord>=3)irngt(3)=xyzzyaaba133(2)
if(nord>=4)irngt(4)=xyzzyaaba133(3)
return
endif
xyzzyaaap133=0
xyzzyaaao133=xyzzyaaap133+1
xyzzyaaah133=xyzzyaaao133
xyzzyaaai133=3
xyzzyaaaa133=xdont(xyzzyaabb133(xyzzyaaao133))+int(real(2*nord,dp)/rea&
&l(xyzzyaaag133+nord,dp)*real((xdont(xyzzyaaba133(3))-xdont(xyzzyaabb1&
&33(xyzzyaaao133))),dp))
if(xyzzyaaaa133>=xdont(xyzzyaaba133(1)))then
xyzzyaaaa133=xdont(xyzzyaabb133(xyzzyaaao133))+int(real(2*nord,dp)/rea&
&l(xyzzyaaag133+nord,dp)* real((xdont(xyzzyaaba133(2))-xdont(xyzzyaabb&
&133(xyzzyaaao133))),dp))
if(xyzzyaaaa133>=xdont(xyzzyaaba133(1)))xyzzyaaaa133=xdont(xyzzyaabb13&
&3(xyzzyaaao133))+int(real(2*nord,dp)/real(xyzzyaaag133+nord,dp)*real(&
&(xdont(xyzzyaaba133(1))-xdont(xyzzyaabb133(xyzzyaaao133))),dp))
endif
xyzzyaaab133=xyzzyaaaa133
if(xdont(xyzzyaaag133)<xyzzyaaaa133)then
xyzzyaaat133=3
do
xyzzyaaat133=xyzzyaaat133+1
if(xdont(xyzzyaaat133)<xyzzyaaaa133)then
if(xyzzyaaat133>=xyzzyaaag133)exit
xyzzyaaai133=xyzzyaaai133+1
xyzzyaaba133(xyzzyaaai133)=xyzzyaaat133
else
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaat133
if(xyzzyaaah133>=nord)exit
endif
enddo
if(xyzzyaaat133<xyzzyaaag133-1)then
do
xyzzyaaat133=xyzzyaaat133+1
if(xdont(xyzzyaaat133)>=xyzzyaaaa133)then
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaat133
elseif(xyzzyaaat133>=xyzzyaaag133)then
exit
endif
enddo
endif
else
do xyzzyaaat133=4,xyzzyaaag133-1
if(xdont(xyzzyaaat133)<xyzzyaaaa133)then
xyzzyaaai133=xyzzyaaai133+1
xyzzyaaba133(xyzzyaaai133)=xyzzyaaat133
else
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaat133
if(xyzzyaaah133>=nord)exit
endif
enddo
if(xyzzyaaat133<xyzzyaaag133-1)then
do
xyzzyaaat133=xyzzyaaat133+1
if(xdont(xyzzyaaat133)>=xyzzyaaaa133)then
if(xyzzyaaat133>=xyzzyaaag133)exit
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaat133
endif
enddo
endif
endif
xyzzyaaaw133=0
xyzzyaaax133=0
xyzzyaaay133=0
xyzzyaaaz133=0
do
if(xyzzyaaah133==nord)exit
if(xyzzyaaay133==xyzzyaaah133.and.xyzzyaaaw133==xyzzyaaai133)then
if(nord>xyzzyaaah133)then
xyzzyaaaf133=xdont(xyzzyaaba133(1))
xyzzyaaav133=1
do xyzzyaaat133=2,xyzzyaaai133
if(xdont(xyzzyaaba133(xyzzyaaat133))>xyzzyaaaf133)then
xyzzyaaaf133=xdont(xyzzyaaba133(xyzzyaaat133))
xyzzyaaav133=xyzzyaaat133
endif
enddo
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(xyzzyaaav133)
xyzzyaaba133(xyzzyaaav133)=xyzzyaaba133(xyzzyaaai133)
xyzzyaaai133=xyzzyaaai133-1
else
xyzzyaaaj133=xyzzyaabb133(xyzzyaaah133)
xyzzyaaae133=xdont(xyzzyaaaj133)
do xyzzyaaat133=1,xyzzyaaah133
if(xdont(xyzzyaabb133(xyzzyaaat133))<xyzzyaaae133)then
xyzzyaaak133=xyzzyaabb133(xyzzyaaat133)
xyzzyaaae133=xdont(xyzzyaaak133)
xyzzyaabb133(xyzzyaaat133)=xyzzyaaaj133
xyzzyaaaj133=xyzzyaaak133
endif
enddo
xyzzyaaah133=xyzzyaaah133-1
endif
endif
xyzzyaaaw133=xyzzyaaax133
xyzzyaaax133=xyzzyaaai133
xyzzyaaay133=xyzzyaaaz133
xyzzyaaaz133=xyzzyaaah133
select case (nord-xyzzyaaah133)
case (2:)
select case (xyzzyaaai133)
case(2)
if(xdont(xyzzyaaba133(1))>=xdont(xyzzyaaba133(2)))then
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(1)
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(2)
else
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(2)
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133 (1)
endif
exit
case (3)
xyzzyaaal133=xyzzyaaba133(1)
xyzzyaaam133=xyzzyaaba133(2)
xyzzyaaan133=xyzzyaaba133(3)
if(xdont(xyzzyaaam133)>xdont(xyzzyaaal133))then
xyzzyaaba133(1)=xyzzyaaam133
xyzzyaaba133(2)=xyzzyaaal133
xyzzyaaam133=xyzzyaaal133
endif
if(xdont(xyzzyaaam133)<xdont(xyzzyaaan133))then
xyzzyaaba133(3)=xyzzyaaam133
xyzzyaaba133(2)=xyzzyaaan133
xyzzyaaam133=xyzzyaaan133
if(xdont(xyzzyaaam133)>xdont(xyzzyaaba133(1)))then
xyzzyaaba133(2)=xyzzyaaba133(1)
xyzzyaaba133(1)=xyzzyaaam133
endif
endif
xyzzyaaai133=0
do xyzzyaaat133=xyzzyaaah133+1,nord
xyzzyaaai133=xyzzyaaai133+1
xyzzyaabb133(xyzzyaaat133)=xyzzyaaba133(xyzzyaaai133)
enddo
xyzzyaaah133=nord
exit
case(4:)
xyzzyaaab133=xyzzyaaaa133
xyzzyaaar133=xyzzyaaai133
xyzzyaaal133=xyzzyaaba133(1)
xyzzyaaam133=xyzzyaaba133(2)
xyzzyaaan133=xyzzyaaba133(xyzzyaaar133)
if(xdont(xyzzyaaam133)>xdont(xyzzyaaal133))then
xyzzyaaba133(1)=xyzzyaaam133
xyzzyaaba133(2)=xyzzyaaal133
xyzzyaaam133=xyzzyaaal133
endif
if(xdont(xyzzyaaam133)<xdont(xyzzyaaan133))then
xyzzyaaba133(xyzzyaaar133)=xyzzyaaam133
xyzzyaaba133(2)=xyzzyaaan133
xyzzyaaam133=xyzzyaaan133
if(xdont(xyzzyaaam133)>xdont(xyzzyaabb133(1)))then
xyzzyaaba133(2)=xyzzyaaba133(1)
xyzzyaaba133(1)=xyzzyaaam133
endif
endif
xyzzyaaap133=xyzzyaaah133
xyzzyaaas133=nord-xyzzyaaah133
xyzzyaaal133=xyzzyaaba133(1)
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaal133
xyzzyaaaa133=xdont(xyzzyaaal133)+int(real(xyzzyaaas133,dp)/real(nord+x&
&yzzyaaas133,dp)* real((xdont(xyzzyaaba133(xyzzyaaar133))-xdont(xyzzya&
&aal133)),dp))
xyzzyaaai133=0
do xyzzyaaat133=2,xyzzyaaar133
if(xdont(xyzzyaaba133(xyzzyaaat133))>=xyzzyaaaa133)then
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(xyzzyaaat133)
if(xyzzyaaah133>=nord)exit
else
xyzzyaaai133=xyzzyaaai133+1
xyzzyaaba133(xyzzyaaai133)=xyzzyaaba133(xyzzyaaat133)
endif
enddo
do xyzzyaaat133=xyzzyaaat133+1,xyzzyaaar133
if(xdont(xyzzyaaba133(xyzzyaaat133))>=xyzzyaaaa133)then
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(xyzzyaaat133)
endif
enddo
end select
case (1)
xyzzyaaaf133=xdont(xyzzyaaba133(1))
xyzzyaaav133=1
do xyzzyaaat133=2,xyzzyaaai133
if(xdont(xyzzyaaba133(xyzzyaaat133))>xyzzyaaaf133)then
xyzzyaaaf133=xdont(xyzzyaaba133(xyzzyaaat133))
xyzzyaaav133=xyzzyaaat133
endif
enddo
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaaba133(xyzzyaaav133)
exit
case (0)
exit
case (-5:-1)
irngt(1)=xyzzyaabb133(1)
do xyzzyaaat133=2,nord
xyzzyaaak133=xyzzyaabb133(xyzzyaaat133)
xyzzyaaac133=xdont(xyzzyaaak133)
do xyzzyaaau133=xyzzyaaat133-1,1,-1
if(xyzzyaaac133>xdont(irngt(xyzzyaaau133)))then
irngt(xyzzyaaau133+1)=irngt(xyzzyaaau133)
else
exit
endif
enddo
irngt(xyzzyaaau133+1)=xyzzyaaak133
enddo
xyzzyaaad133=xdont(irngt(nord))
do xyzzyaaat133=nord+1,xyzzyaaah133
if(xdont(xyzzyaabb133(xyzzyaaat133))>xyzzyaaad133)then
xyzzyaaac133=xdont(xyzzyaabb133(xyzzyaaat133))
do xyzzyaaau133=nord-1,1,-1
if (xyzzyaaac133<=xdont(irngt(xyzzyaaau133)))exit
irngt(xyzzyaaau133+1)=irngt(xyzzyaaau133)
enddo
irngt(xyzzyaaau133+1)=xyzzyaabb133(xyzzyaaat133)
xyzzyaaad133=xdont(irngt(nord))
endif
enddo
return
case(:-6)
xyzzyaaao133=xyzzyaaap133+1
xyzzyaaaq133=(xyzzyaaah133+xyzzyaaao133)/2
xyzzyaaar133=xyzzyaaah133
if(xdont(xyzzyaabb133(xyzzyaaaq133))>xdont(xyzzyaabb133(xyzzyaaao133))&
&)then
xyzzyaaak133=xyzzyaabb133(xyzzyaaao133)
xyzzyaabb133(xyzzyaaao133)=xyzzyaabb133(xyzzyaaaq133)
xyzzyaabb133(xyzzyaaaq133)=xyzzyaaak133
endif
if(xdont(xyzzyaabb133(xyzzyaaaq133))<xdont(xyzzyaabb133(xyzzyaaar133))&
&)then
xyzzyaaak133=xyzzyaabb133(xyzzyaaar133)
xyzzyaabb133(xyzzyaaar133)=xyzzyaabb133(xyzzyaaaq133)
xyzzyaabb133(xyzzyaaaq133)=xyzzyaaak133
if(xdont(xyzzyaabb133(xyzzyaaaq133))>xdont(xyzzyaabb133(xyzzyaaao133))&
&)then
xyzzyaaak133=xyzzyaabb133(xyzzyaaao133)
xyzzyaabb133(xyzzyaaao133)=xyzzyaabb133(xyzzyaaaq133)
xyzzyaabb133(xyzzyaaaq133)=xyzzyaaak133
endif
endif
if(xyzzyaaar133<=3)exit
xyzzyaaaa133=xdont(xyzzyaabb133(1))+int(real(nord)/real(xyzzyaaah133+n&
&ord)* real((xdont(xyzzyaabb133(xyzzyaaar133))-xdont(xyzzyaabb133(1)))&
&,dp))
if(xyzzyaaap133>0)then
if(xyzzyaaaa133<=xyzzyaaab133)xyzzyaaaa133=xyzzyaaab133+int(real(2*nor&
&d-xyzzyaaap133,dp)/real(xyzzyaaah133+nord,dp)* real((xdont(xyzzyaabb1&
&33(xyzzyaaar133))-xyzzyaaab133),dp))
else
xyzzyaaao133=1
endif
xyzzyaaai133=0
xyzzyaaah133=xyzzyaaap133
if(xdont(xyzzyaabb133(xyzzyaaar133))<xyzzyaaaa133)then
xyzzyaaat133=xyzzyaaap133
do
xyzzyaaat133=xyzzyaaat133+1
if(xdont(xyzzyaabb133(xyzzyaaat133))<xyzzyaaaa133)then
xyzzyaaai133=xyzzyaaai133+1
xyzzyaaba133(xyzzyaaai133)=xyzzyaabb133(xyzzyaaat133)
if(xyzzyaaat133>=xyzzyaaar133)exit
else
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaabb133(xyzzyaaat133)
if(xyzzyaaah133>=nord)exit
endif
enddo
if(xyzzyaaat133<xyzzyaaar133)then
do
xyzzyaaat133=xyzzyaaat133+1
if(xdont(xyzzyaabb133(xyzzyaaat133))>=xyzzyaaaa133)then
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaabb133(xyzzyaaat133)
else
if(xyzzyaaat133>=xyzzyaaar133)exit
endif
enddo
endif
else
do xyzzyaaat133=xyzzyaaao133,xyzzyaaar133
if(xdont(xyzzyaabb133(xyzzyaaat133))<xyzzyaaaa133)then
xyzzyaaai133=xyzzyaaai133+1
xyzzyaaba133(xyzzyaaai133)=xyzzyaabb133(xyzzyaaat133)
else
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaabb133(xyzzyaaat133)
if(xyzzyaaah133>=nord)exit
endif
enddo
do xyzzyaaat133=xyzzyaaat133+1,xyzzyaaar133
if(xdont(xyzzyaabb133(xyzzyaaat133))>=xyzzyaaaa133)then
xyzzyaaah133=xyzzyaaah133+1
xyzzyaabb133(xyzzyaaah133)=xyzzyaabb133(xyzzyaaat133)
endif
enddo
endif
end select
enddo
irngt(1)=xyzzyaabb133(1)
do xyzzyaaat133=2,nord
xyzzyaaak133=xyzzyaabb133(xyzzyaaat133)
xyzzyaaac133=xdont(xyzzyaaak133)
do xyzzyaaau133=xyzzyaaat133-1,1,-1
if(xyzzyaaac133>xdont(irngt(xyzzyaaau133)))then
irngt(xyzzyaaau133+1)=irngt(xyzzyaaau133)
else
exit
endif
enddo
irngt(xyzzyaaau133+1)=xyzzyaaak133
enddo
end subroutine partial_rank_int
subroutine resize_pointer_dble1(dims,pt,init)
implicit none
integer,intent(in) :: dims(1)
real(dp),intent(in),optional :: init
integer xyzzyaaaa134(1),xyzzyaaab134
real(dp),parameter :: xyzzyaaac134=0.d0
real(dp),pointer :: pt(:),xyzzyaaad134(:)
if(.not.associated(pt))then
allocate(pt(dims(1)),stat=xyzzyaaab134)
call check_alloc(xyzzyaaab134,'RESIZE_POINTER_DBLE1','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac134
endif
return
endif
xyzzyaaaa134=shape(pt)
if(all(xyzzyaaaa134==dims))return
allocate(xyzzyaaad134(dims(1)),stat=xyzzyaaab134)
if(any(xyzzyaaaa134<dims))then
if(present(init))then
xyzzyaaad134=init
else
xyzzyaaad134=xyzzyaaac134
endif
endif
xyzzyaaad134(1:min(xyzzyaaaa134(1),dims(1)))=pt(1:min(xyzzyaaaa134(1),&
&dims(1)))
deallocate(pt)
pt=>xyzzyaaad134
nullify(xyzzyaaad134)
end subroutine resize_pointer_dble1
subroutine resize_pointer_dble2(dims,pt,init)
implicit none
integer,intent(in) :: dims(2)
real(dp),intent(in),optional :: init
integer xyzzyaaaa135(2),xyzzyaaab135
real(dp),parameter :: xyzzyaaac135=0.d0
real(dp),pointer :: pt(:,:),xyzzyaaad135(:,:)
if(.not.associated(pt))then
allocate(pt(dims(1),dims(2)),stat=xyzzyaaab135)
call check_alloc(xyzzyaaab135,'RESIZE_POINTER_DBLE2','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac135
endif
return
endif
xyzzyaaaa135=shape(pt)
if(all(xyzzyaaaa135==dims))return
allocate(xyzzyaaad135(dims(1),dims(2)),stat=xyzzyaaab135)
if(any(xyzzyaaaa135<dims))then
if(present(init))then
xyzzyaaad135=init
else
xyzzyaaad135=xyzzyaaac135
endif
endif
xyzzyaaad135(1:min(xyzzyaaaa135(1),dims(1)),1:min(xyzzyaaaa135(2),dims&
&(2)))=pt(1:min(xyzzyaaaa135(1),dims(1)),1:min(xyzzyaaaa135(2),dims(2)&
&))
deallocate(pt)
pt=>xyzzyaaad135
nullify(xyzzyaaad135)
end subroutine resize_pointer_dble2
subroutine resize_pointer_int1(dims,pt,init)
implicit none
integer,intent(in) :: dims(1)
integer,intent(in),optional :: init
integer xyzzyaaaa136(1),xyzzyaaab136
integer,parameter :: xyzzyaaac136=0
integer,pointer :: pt(:),xyzzyaaad136(:)
if(.not.associated(pt))then
allocate(pt(dims(1)),stat=xyzzyaaab136)
call check_alloc(xyzzyaaab136,'RESIZE_POINTER_INT1','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac136
endif
return
endif
xyzzyaaaa136=shape(pt)
if(all(xyzzyaaaa136==dims))return
allocate(xyzzyaaad136(dims(1)),stat=xyzzyaaab136)
if(any(xyzzyaaaa136<dims))then
if(present(init))then
xyzzyaaad136=init
else
xyzzyaaad136=xyzzyaaac136
endif
endif
xyzzyaaad136(1:min(xyzzyaaaa136(1),dims(1)))=pt(1:min(xyzzyaaaa136(1),&
&dims(1)))
deallocate(pt)
pt=>xyzzyaaad136
nullify(xyzzyaaad136)
end subroutine resize_pointer_int1
subroutine resize_pointer_int2(dims,pt,init)
implicit none
integer,intent(in) :: dims(2)
integer,intent(in),optional :: init
integer xyzzyaaaa137(2),xyzzyaaab137
integer,parameter :: xyzzyaaac137=0
integer,pointer :: pt(:,:),xyzzyaaad137(:,:)
if(.not.associated(pt))then
allocate(pt(dims(1),dims(2)),stat=xyzzyaaab137)
call check_alloc(xyzzyaaab137,'RESIZE_POINTER_INT2','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac137
endif
return
endif
xyzzyaaaa137=shape(pt)
if(all(xyzzyaaaa137==dims))return
allocate(xyzzyaaad137(dims(1),dims(2)),stat=xyzzyaaab137)
if(any(xyzzyaaaa137<dims))then
if(present(init))then
xyzzyaaad137=init
else
xyzzyaaad137=xyzzyaaac137
endif
endif
xyzzyaaad137(1:min(xyzzyaaaa137(1),dims(1)),1:min(xyzzyaaaa137(2),dims&
&(2)))=pt(1:min(xyzzyaaaa137(1),dims(1)),1:min(xyzzyaaaa137(2),dims(2)&
&))
deallocate(pt)
pt=>xyzzyaaad137
nullify(xyzzyaaad137)
end subroutine resize_pointer_int2
subroutine resize_pointer_int3(dims,pt,init)
implicit none
integer,intent(in) :: dims(3)
integer,intent(in),optional :: init
integer xyzzyaaaa138(3),xyzzyaaab138
integer,parameter :: xyzzyaaac138=0
integer,pointer :: pt(:,:,:),xyzzyaaad138(:,:,:)
if(.not.associated(pt))then
allocate(pt(dims(1),dims(2),dims(3)),stat=xyzzyaaab138)
call check_alloc(xyzzyaaab138,'RESIZE_POINTER_INT3','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac138
endif
return
endif
xyzzyaaaa138=shape(pt)
if(all(xyzzyaaaa138==dims))return
allocate(xyzzyaaad138(dims(1),dims(2),dims(3)),stat=xyzzyaaab138)
if(any(xyzzyaaaa138<dims))then
if(present(init))then
xyzzyaaad138=init
else
xyzzyaaad138=xyzzyaaac138
endif
endif
xyzzyaaad138(1:min(xyzzyaaaa138(1),dims(1)),1:min(xyzzyaaaa138(2),dims&
&(2)),1:min(xyzzyaaaa138(3),dims(3)))=pt(1:min(xyzzyaaaa138(1),dims(1)&
&),1:min(xyzzyaaaa138(2),dims(2)),1:min(xyzzyaaaa138(3),dims(3)))
deallocate(pt)
pt=>xyzzyaaad138
nullify(xyzzyaaad138)
end subroutine resize_pointer_int3
subroutine resize_pointer_bool1(dims,pt,init)
implicit none
integer,intent(in) :: dims(1)
logical,intent(in),optional :: init
integer xyzzyaaaa139(1),xyzzyaaab139
logical,parameter :: xyzzyaaac139=.false.
logical,pointer :: pt(:),xyzzyaaad139(:)
if(.not.associated(pt))then
allocate(pt(dims(1)),stat=xyzzyaaab139)
call check_alloc(xyzzyaaab139,'RESIZE_POINTER_BOOL1','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac139
endif
return
endif
xyzzyaaaa139=shape(pt)
if(all(xyzzyaaaa139==dims))return
allocate(xyzzyaaad139(dims(1)),stat=xyzzyaaab139)
if(any(xyzzyaaaa139<dims))then
if(present(init))then
xyzzyaaad139=init
else
xyzzyaaad139=xyzzyaaac139
endif
endif
xyzzyaaad139(1:min(xyzzyaaaa139(1),dims(1)))=pt(1:min(xyzzyaaaa139(1),&
&dims(1)))
deallocate(pt)
pt=>xyzzyaaad139
nullify(xyzzyaaad139)
end subroutine resize_pointer_bool1
subroutine resize_pointer_bool2(dims,pt,init)
implicit none
integer,intent(in) :: dims(2)
logical,intent(in),optional :: init
integer xyzzyaaaa140(2),xyzzyaaab140
logical,parameter :: xyzzyaaac140=.false.
logical,pointer :: pt(:,:),xyzzyaaad140(:,:)
if(.not.associated(pt))then
allocate(pt(dims(1),dims(2)),stat=xyzzyaaab140)
call check_alloc(xyzzyaaab140,'RESIZE_POINTER_BOOL2','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac140
endif
return
endif
xyzzyaaaa140=shape(pt)
if(all(xyzzyaaaa140==dims))return
allocate(xyzzyaaad140(dims(1),dims(2)),stat=xyzzyaaab140)
if(any(xyzzyaaaa140<dims))then
if(present(init))then
xyzzyaaad140=init
else
xyzzyaaad140=xyzzyaaac140
endif
endif
xyzzyaaad140(1:min(xyzzyaaaa140(1),dims(1)),1:min(xyzzyaaaa140(2),dims&
&(2)))=pt(1:min(xyzzyaaaa140(1),dims(1)),1:min(xyzzyaaaa140(2),dims(2)&
&))
deallocate(pt)
pt=>xyzzyaaad140
nullify(xyzzyaaad140)
end subroutine resize_pointer_bool2
subroutine resize_pointer_bool3(dims,pt,init)
implicit none
integer,intent(in) :: dims(3)
logical,intent(in),optional :: init
integer xyzzyaaaa141(3),xyzzyaaab141
logical,parameter :: xyzzyaaac141=.false.
logical,pointer :: pt(:,:,:),xyzzyaaad141(:,:,:)
if(.not.associated(pt))then
allocate(pt(dims(1),dims(2),dims(3)),stat=xyzzyaaab141)
call check_alloc(xyzzyaaab141,'RESIZE_POINTER_BOOL3','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac141
endif
return
endif
xyzzyaaaa141=shape(pt)
if(all(xyzzyaaaa141==dims))return
allocate(xyzzyaaad141(dims(1),dims(2),dims(3)),stat=xyzzyaaab141)
if(any(xyzzyaaaa141<dims))then
if(present(init))then
xyzzyaaad141=init
else
xyzzyaaad141=xyzzyaaac141
endif
endif
xyzzyaaad141(1:min(xyzzyaaaa141(1),dims(1)),1:min(xyzzyaaaa141(2),dims&
&(2)),1:min(xyzzyaaaa141(3),dims(3)))=pt(1:min(xyzzyaaaa141(1),dims(1)&
&),1:min(xyzzyaaaa141(2),dims(2)),1:min(xyzzyaaaa141(3),dims(3)))
deallocate(pt)
pt=>xyzzyaaad141
nullify(xyzzyaaad141)
end subroutine resize_pointer_bool3
subroutine resize_pointer_char1(sz,dims,pt,init)
implicit none
integer,intent(in) :: sz,dims(1)
character(sz),intent(in),optional :: init
integer xyzzyaaaa142(1),xyzzyaaab142
character(0),parameter :: xyzzyaaac142=''
character(sz),pointer :: pt(:),xyzzyaaad142(:)
if(.not.associated(pt))then
allocate(pt(dims(1)),stat=xyzzyaaab142)
call check_alloc(xyzzyaaab142,'RESIZE_POINTER_CHAR1','pt')
if(present(init))then
pt=init
else
pt=xyzzyaaac142
endif
return
endif
xyzzyaaaa142=shape(pt)
if(all(xyzzyaaaa142==dims))return
allocate(xyzzyaaad142(dims(1)),stat=xyzzyaaab142)
if(any(xyzzyaaaa142<dims))then
if(present(init))then
xyzzyaaad142=init
else
xyzzyaaad142=xyzzyaaac142
endif
endif
xyzzyaaad142(1:min(xyzzyaaaa142(1),dims(1)))=pt(1:min(xyzzyaaaa142(1),&
&dims(1)))
deallocate(pt)
pt=>xyzzyaaad142
nullify(xyzzyaaad142)
end subroutine resize_pointer_char1
subroutine get_numerical_orbmask(norb,nemax,nspin,ndet,orbmap,orbmask)
implicit none
integer,intent(in) :: norb,nemax,nspin,ndet,orbmap(nemax,nspin,ndet)
logical,intent(out) :: orbmask(norb,nspin)
integer xyzzyaaaa143,xyzzyaaab143,xyzzyaaac143,xyzzyaaad143
orbmask=.false.
do xyzzyaaab143=1,nspin
do xyzzyaaaa143=1,ndet
do xyzzyaaac143=1,nemax
xyzzyaaad143=orbmap(xyzzyaaac143,xyzzyaaab143,xyzzyaaaa143)
if(xyzzyaaad143==0)cycle
orbmask(xyzzyaaad143,xyzzyaaab143)=.true.
enddo
enddo
enddo
end subroutine get_numerical_orbmask
subroutine get_numerical_orbrange(norb,nspin,orbmask,norbrange,orbrang&
&e)
implicit none
integer,intent(in) :: norb,nspin
integer,intent(out) :: norbrange(nspin),orbrange(2,norb,nspin)
logical,intent(in) :: orbmask(norb,nspin)
integer xyzzyaaaa144,xyzzyaaab144,xyzzyaaac144,xyzzyaaad144
norbrange=0
orbrange=0
do xyzzyaaaa144=1,nspin
xyzzyaaad144=0
do
xyzzyaaac144=0
do xyzzyaaab144=xyzzyaaad144+1,norb
if(orbmask(xyzzyaaab144,xyzzyaaaa144))then
xyzzyaaac144=xyzzyaaab144
exit
endif
enddo
if(xyzzyaaac144==0)exit
xyzzyaaad144=0
do xyzzyaaab144=xyzzyaaac144+1,norb
if(.not.orbmask(xyzzyaaab144,xyzzyaaaa144))then
xyzzyaaad144=xyzzyaaab144-1
exit
endif
enddo
norbrange(xyzzyaaaa144)=norbrange(xyzzyaaaa144)+1
if(xyzzyaaad144>0)then
orbrange(1:2,norbrange(xyzzyaaaa144),xyzzyaaaa144)=(/xyzzyaaac144,xyzz&
&yaaad144/)
else
orbrange(1:2,norbrange(xyzzyaaaa144),xyzzyaaaa144)=(/xyzzyaaac144,norb&
&/)
endif
if(xyzzyaaad144==0)exit
enddo
enddo
end subroutine get_numerical_orbrange
logical function isdiag_int(n,a)
implicit none
integer,intent(in) :: n,a(n,n)
integer xyzzyaaaa145,xyzzyaaab145
isdiag_int=.true.
do xyzzyaaaa145=2,n
do xyzzyaaab145=1,xyzzyaaaa145-1
if(a(xyzzyaaaa145,xyzzyaaab145)/=0.or.a(xyzzyaaab145,xyzzyaaaa145)/=0)&
&then
isdiag_int=.false.
return
endif
enddo
enddo
end function isdiag_int
function approx_equal_sp(x,y,epsilon) result(approx_equal)
implicit none
real(sp),intent(in) :: x,y,epsilon
real(sp) xyzzyaaaa146
logical approx_equal
approx_equal=.false.
xyzzyaaaa146=abs(x-y)
if(xyzzyaaaa146<=epsilon*max(1._sp,abs(x),abs(y)))approx_equal=.true.
end function approx_equal_sp
function approx_equal_dp(x,y,epsilon) result(approx_equal)
implicit none
real(dp),intent(in) :: x,y,epsilon
real(dp) xyzzyaaaa147
logical approx_equal
approx_equal=.false.
xyzzyaaaa147=abs(x-y)
if(xyzzyaaaa147<=epsilon*max(1.d0,abs(x),abs(y)))approx_equal=.true.
end function approx_equal_dp
end module slaarnabt
