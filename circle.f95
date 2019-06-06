         subroutine ccenter(p1,p2,p3,xc,yc,t1)
!-----------------------------------------------------------
! input ...
! p1,p2,p3 - three points to fit a circle
!if p1 and p2 are the edges of a cord the subroutine will compute the tangent angle as well
! output ...
! xc,yc - circle center coordinates
!t1=tangent
! comments ...
! the D matrix can be used to compute the radius
! during the calculation
!===========================================================
                implicit none
                                Double precision :: p1(1,2),p2(1,2),p3(1,2)
                                Double precision :: A(3,3),B(3,3),C(3,3),D(3,3)
                                Double precision :: detA,detB,detC,detD 
                                Double precision :: xc,yc,r,m1,t1,r2d                   

                A(1,1)=p1(1,1)
                A(1,2)=p1(1,2)
                A(1,3)=1

                A(2,1)=p2(1,1)
                A(2,2)=p2(1,2)
                A(2,3)=1

                A(3,1)=p3(1,1)
                A(3,2)=p3(1,2)
                A(3,3)=1

                call determinant(A,detA)
                !print *, 'determinante di A', detA


                B(1,1)=p1(1,1)**2+p1(1,2)**2
                B(1,2)=p1(1,2)
                B(1,3)=1

                B(2,1)=p2(1,1)**2+p2(1,2)**2
                B(2,2)=p2(1,2)
                B(2,3)=1

                B(3,1)=p3(1,1)**2+p3(1,2)**2
                B(3,2)=p3(1,2)
                B(3,3)=1

                B=-B
                call determinant(B,detB)

                !print *, 'determinante di B', detB

                C(1,1)=p1(1,1)**2+p1(1,2)**2
                C(1,2)=p1(1,1)
                C(1,3)=1

                C(2,1)=p2(1,1)**2+p2(1,2)**2
                C(2,2)=p2(1,1)
                C(2,3)=1

                C(3,1)=p3(1,1)**2+p3(1,2)**2
                C(3,2)=p3(1,1)
                C(3,3)=1

                call determinant(C,detC)

                !print *, 'determinante di C', detC

!                D(1,1)=p1(1,1)**2+p1(1,2)**2
!                D(1,2)=p1(1,1)
!                D(1,3)=p1(1,2)

!                D(2,1)=p2(1,1)**2+p2(1,2)**2
!                D(2,2)=p2(1,1)
!                D(2,3)=p2(1,2)

!                D(3,1)=p3(1,1)**2+p3(1,2)**2
!                D(3,2)=p3(1,1)
!                D(3,3)=p3(1,2)

!                D=-D


!                call determinant(D,detD)

                !print *, 'determinante di D', detD   

                !KNOWING THESE DETERMINANTS IT IS POSSIBLE TO EXTRACT THE CENTER COORDINATES OF THE FITTING CIRCLE
                xc=-(detB)/(2*detA)
                yc=-(detC)/(2*detA)
!                r=sqrt((detB**2+detC**2-4*detA*detD)/(4*detA**2))

                m1=(p1(1,1)-xc)/(p1(1,2)-yc);

                t1=-1/((m1))
                if (abs(m1).lt.0.1) then
                t1=0
                !print *, t1,m1, abs(m1)
                end if

                !CONVERT TANGENT ANGLE TO RADIANTS
                r2d= 57.295779513;
                t1=datan(t1)*r2d
                                end subroutine ccenter
