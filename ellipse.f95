!-----------------------------------------------------------
! input ...
! P(r1,c1) - ellipse data points 
! san - number of points used to draw the fitted ellipse 
! output ...
! xx - coordinates of the fitted ellipse
!===========================================================



        subroutine ellipse(P,counter,r1,c1,san,xx)
            implicit none
            integer:: counter,i,j,r1,c1,k
            integer:: sh(2),sh2(2),shAT(2),san
            Double Precision :: ris,pi,phi, up, down1, down2,res1,res2,R,step
            Double precision,intent(in) :: P(r1,c1),zav
            Double precision ::S(6,6),Sinv(6,6),C(6,6),SC(6,6),evec(6),xx(san,2)
            Double precision, Allocatable :: A(:,:),At(:,:)
            character (len=1):: jobvl, jobvr
            Double precision :: WR(6),WI(6),VL(6,6),VR(6,6)
            Double precision, ALLOCATABLE :: WORK(:)
            Integer :: LWORK,  INFO, meval
            Double precision :: coeff1, coeff2, coeff3,coeff4,coeff5,coeff6,den,x0,y0
            integer :: splt,splt2
           !WR and WI contain the real and imaginary parts, respectively, of the computed eigenvalues.
           !VL the left eigenvectors u(j) are stored one after another in the columns of VL, in te same order as their eigenvalues
           !VR the right eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues




           pi=4*atan(1.0d0)
           C=0.0d0
           C(1,3)=2
           C(2,2)=-1
           C(3,1)=2


           ALLOCATE(A(counter,6),At(6,counter))
                      do i=1,counter
                        ! print *, P(i,1),P(i,2),P(i,3)

                         A(i,1)=P(i,1)*P(i,1)
                         A(i,2)=P(i,1)*P(i,2)
                         A(i,3)=P(i,2)*P(i,2)
                         A(i,4)=P(i,1)
                         A(i,5)=P(i,2)
                         A(i,6)=1
                       end do
           sh=shape(A)
           At=transpose(A)
         !    print *, sh, 'A'
                do i=1,sh(2)
                do j=1,sh(1)
                At(i,j)=A(j,i)
                end do
                end do
          shAT=shape(At)

                S=0
                j=1
                !transposed matrix x matrix
                do i=1,shAT(1)
                do j=1,sh(2)
                ris=0
                do k=1,sh(1)
                ris = ris+At(i,k)*A(k,j)
                S(i,j)=ris
                end do
                end do
                end do

           sh=shape(S)
           call inv(S,Sinv,6)
           SC=matmul(Sinv,C)
           sh2=shape(SC)
           LWORK=4*sh2(1)
           ALLOCATE(Work(MAX(1,LWORK)))
           JOBVL='V'
           JOBVR='V'

           call dgeev(JOBVL,JOBVR,sh2(1),SC,sh2(1),WR,WI,VL,sh2(1),VR,sh2(1),WORK,LWORK,INFO)

           meval=int(Maxloc(WR,dim=1))

           evec=VR(:,meval)
           coeff1=evec(1)
           coeff2=evec(2)/2
           coeff3=evec(3)
           coeff4=evec(4)/2
           coeff5=evec(5)/2
           coeff6=evec(6)
!!!!!!!!!!!FIND CENTER OF THE ELLIPSE

           den=coeff2*coeff2-coeff1*coeff3;
           x0=(coeff3*coeff4-coeff2*coeff5)/den;
           y0=(coeff1*coeff5-coeff2*coeff4)/den;
         ! print *, 'center of the ellipse', x0,y0,P(1,3)
!!!!!!!!!!!FIND ROTATION ANGLE OF THE ELLIPSE
                if (coeff2.eq. 0)then

                    if (coeff1.gt.coeff3) then
                        phi=0
                    else
                        phi=pi/2
                     endif
                else
                    if (coeff1.gt.coeff3) then
                      phi=atan(2*coeff2/(coeff1-coeff3))/2
                    else
                      phi=pi/2+atan(2*coeff2/(coeff1-coeff3))/2
                    endif

                 endif
        !    print *, 'rotation angle'
        !    print *, phi       
!!!!!!!!!!!!FIND AXES OF THE ELLIPSE
        up = 2*(coeff1*coeff5*coeff5+coeff3*coeff4*coeff4+coeff6*coeff2*coeff2-2*coeff2*coeff4*coeff5-coeff1*coeff3*coeff6)
        down1=(coeff2*coeff2-coeff1*coeff3)*((coeff3-coeff1)*sqrt(1+4*coeff2*coeff2&
        /((coeff1-coeff3)*(coeff1-coeff3)))-(coeff3+coeff1))
        down2=(coeff2*coeff2-coeff1*coeff3)*((coeff1-coeff3)*sqrt(1+4*coeff2*coeff2&
        /((coeff1-coeff3)*(coeff1-coeff3)))-(coeff3+coeff1))
        res2=sqrt(up/down1)
        res1=sqrt(up/down2)
        !print *, 'semiassi' 
        !print *, res1,res2
        !print *, 'area',pi*res1*res2
!!!!!!!!!!! DRAW ELLIPSE
          step=(2*pi/dble(san))
          !print *, 'STEP',step,san*step,'SAN',san
          splt=0

          !print *, 'splt',splt,2*pi,step, int(2*pi/step) 
          !ALLOCATE (XX(san), YY(san))
          R=0.0d0+step
          do while (R.lt.2*pi)
          splt=splt+1
          xx(splt,1) = x0 + res1*cos(R)*cos(phi) - res2*sin(R)*sin(phi)
          xx(splt,2) = y0 + res1*cos(R)*sin(phi) + res2*sin(R)*cos(phi)
          !write(185,*) xx(splt,1),xx(splt,2)
          R=R+step
          end do
          !print *, splt 
           DEALLOCATE(Work,A,At)
           end subroutine ellipse
