c * * * * * * * * * * * * * * * * * * * * * * * * *
c    Driver for ROCK4 (or ROCK2) at Brusselator-2dim problem
c * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    This driver shows how to use ROCK4 (or ROCK2). It solves a
c    system of ODEs resulting from the 2-dimensional space 
c    discretization of the Brusselator equations (u=u(x,y,t),v=v(x,y,t)):
c     
c    u_t=1+u^2*v-4.4*u+0.1*(u_{xx}+u_{yy})+f(x,y,t)
c    v_t=3.4*u-u^2*v+0.1*(v_{xx}+v_{yy})     for t>=0, 0<= x <= 1, 0<= y <= 1
c
c    with initial conditions 
c
c    u(x,y,0)=22*y*(1-y)^{3/2}  v(x,y,0)=27*x*(1-x)^{3/2}
c
c    and periodic boundary conditions
c
c    u(x+1,y,t)=u(x,y,t),  v(x,y+1,t)=v(x,y,t).
c
c    The function f is defined by (inhomogeneity)
c
c    f(x,y,t)=5 if (x-0.3)^2+(y-0.6)^2<= 0.1^2 and t>=1.1
c            =0 else
c
c    We discretize the space variables with
c    x_i=i/(N+1), y_i=i/(N+1) for i=0,1,...,N,
c    with N=128. We obtain a system of 32768
c    equations. The spectral radius of the Jacobian 
c    can be estimated with the Gershgorin theorem   
c    (13200 is an estimation for it). Thus we 
c    provide an external function RHO, giving 
c    the spectral radius of the Jacobian matrix.
c    As output point we choose t_out=1.5 
c
c--------------------------------------------------------
c ----- to integrate with rock2.f ----- 
c     include 'rock2.f'
c  
      include 'rock4.f'     
      implicit double precision (a-h,o-z)
      parameter (nsd=128,neqn=nsd*nsd*2)
c--------------------------------------------------------
c      Work is of length 7*neqn because the spectral radius 
c      is computed externally (otherwise should be 8*neqn).
c      If integrating with ROCK2 define work(4*neqn).
c-------------------------------------------------------- 
c ----- to integrate with rock2.f 
c     dimension y(neqn),work(4*neqn) 
c 
      dimension y(neqn),work(7*neqn)
      integer iwork(12),idid
      external fbrus
c --- common parameters for the problem -----
      common/trans/alf,ns,nssq,nsnsm1,nsm1sq
c ----- file for solution -----
       open(8,file='sol.out')
       rewind 8
c ----- dimensions -----
      ns=nsd
      nssq=ns*ns
      nsnsm1=ns*(ns-1)
      n=2*nssq 
      alf=1.0d-1
      alph=alf
c--------------------------------------------------------
c     Initialize iwork: 
c      iwork(1)=1  RHO returns an upper bound for the spectral radius.
c      iwork(2)=1  The Jacobian is constant (RHO is called once).
c      iwork(3)=0  Return and solution at tend.
c      iwork(4)=0  Atol and rtol are scalars.
c--------------------------------------------------------
      iwork(1)=1
      iwork(2)=1
      iwork(3)=0
      iwork(4)=0
c ----- initial and end point of integration -----
      t=0.0d0
      tend=1.5d0
c ----- initial values -----
      ans=ns
      do j=1,ns
        yy=(j-1)/ans
        do i=1,ns
          y((j-1)*ns+i)=22.d0*yy*(1.d0-yy)**(1.5d0)
        end do
      end do
      do i=1,ns
        xx=(i-1)/ans
        do j=1,ns
           y((j-1)*ns+i+nssq)=27.d0*xx*(1.d0-xx)**(1.5d0)
        end do
      end do
c ----- required tolerance -----
      rtol=0.1d0**2
      atol=rtol
c ----- initial step size -----
      h=1.0d-4
c ----- integration -----
      write(6,*) 'Integration of the 2-dim Brusselator problem' 
c ----- to integrate with rock2.f     
c      call rock2(neqn,t,tend,h,y,fbrus,atol,rtol,work,
c     &           iwork,idid) 
c
c ----- call of the subroutine rock4 -----
      call rock4(neqn,t,tend,h,y,fbrus,atol,rtol,work,
     &           iwork,idid)    
c ----- print solution -----
      do j=1,neqn,267
        write (8,*) y(j)
      end do
c ----- print statistics -----
      write(6,*) 'Solution is tabulated in file sol.out'
      write(6,*) 'The value of IDID is',idid
      write(6,*) 'Max estimation of the spectral radius=',iwork(11)
      write(6,*) 'Min estimation of the spectral radius=',iwork(12)
      write(6,*) 'Max number of stages used=',iwork(10)
      write(6,*) 'Number of f eval. for the spectr. radius=',iwork(9)
      write (6,91) iwork(5),iwork(6),iwork(7),iwork(8)
 91   format(' Number of f evaluations=',i5,' steps=',i4,
     &        ' accpt=',i4,' rejct=',i3) 
c--------------------------------------------------------
c     End of main program
c--------------------------------------------------------
      end      
c--------------------------------------------------------
c     The subroutine RHO gives an estimation of the spectral 
c     radius of the Jacobian matrix of the problem. This
c     is a bound for the whole interval and thus RHO is called
c     once.
c--------------------------------------------------------
      double precision function rho(neqn,t,y)
      implicit double precision (a-h,o-z)
      common/trans/alf,ns,nssq,nsnsm1,nsm1sq
        rho = 8.0d0*nssq*alf + 2.d0 
      return
      end 
c--------------------------------------------------------
c     The subroutine FBRUS compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine fbrus(neqn,x,y,f)
c ----- brusselator with diffusion in 2 dim. space -----
      implicit double precision (a-h,o-z)
      dimension y(neqn),f(neqn)
      common/trans/alf,ns,nssq,nsnsm1,nsm1sq
c ----- constants for inhomogenity -----
      ans=ns
      radsq=0.1d0**2
      if(x.ge.1.1d0)then
        bet=5.00d0
      else
        bet=0.00d0
      end if
c ----- big loop -----
      do i=1,nssq
c ----- left neighbour -----
         if(mod(i,ns).eq.1)then
            uleft=y(i+ns-1)
            vleft=y(nssq+i+ns-1)
         else
            uleft=y(i-1)
            vleft=y(nssq+i-1)
         end if
c ----- right neighbour -----
         if(mod(i,ns).eq.0)then
            uright=y(i-ns+1)
            vright=y(nssq+i-ns+1)
         else
            uright=y(i+1)
            vright=y(nssq+i+1)
         end if
c ----- lower neighbour -----
         if(i.le.ns)then
            ulow=y(i+nsnsm1)
            vlow=y(nssq+i+nsnsm1)
         else
            ulow=y(i-ns)
            vlow=y(nssq+i-ns)
         end if
c ----- upper neighbour -----
         if(i.gt.nsnsm1)then
            uup=y(i-nsnsm1)
            vup=y(nssq+i-nsnsm1)
         else
            uup=y(i+ns)
            vup=y(nssq+i+ns)
         end if
c ----- the derivative -----
         uij=y(i)
         vij=y(i+nssq)
         f(i)=1.d0+uij*uij*vij-4.4d0*uij
     &        +alf*nssq*(uleft+uright+ulow+uup-4.d0*uij)
         f(i+nssq)=3.4d0*uij - uij*uij*vij
     &        +alf*nssq*(vleft+vright+vlow+vup-4.d0*vij)
c ----- inhomogenity -----
         iy=(i-1)/ns+1
         ix=i-(iy-1)*ns
         yy=iy/ans
         xx=ix/ans
         if(((xx-0.3d0)**2+(yy-0.6d0)**2).le.radsq)then
           f(i)=f(i)+bet
         end if
      end do
      return
      end  
      
      


