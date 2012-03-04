c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Driver for ROCK4 (or ROCK2) at FitzHugh and Nagumo problem
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    This driver shows how to use ROCK4 (or ROCK2). It solves a
c    system of ODEs resulting from the space discretization of 
c    the FitzHugh and Nagumo equations (u=u(x,t),v=v(x,t)):
c     
c    u_t=u_{xx}-f(u)-v  for 0<= t <=400, 0<= x <= 100
c    v_t=n(u-b*v)       where a=0.139, n=0.008, b=2.54
c
c    with initial condition 
c
c    u(x,0)=v(x,0)=0
c
c    and boundary conditions
c
c    u_t(0,t)=-0.3,  u_t(100,t)=0
c
c    The function f is defined by 
c
c    f(u)=u*(u-a)*(u-1)
c
c    We discretize the space variables with
c    x_i=(2*i+1)/4, for i=0,1,...,199.
c    We obtain a system of 400 equations.
c    The spectral radius of the Jacobian can
c    be estimated with the Gerhgorin theorem   
c    Thus we provide an external function RHO,
c    giving the spectral radius of the Jacobian matrix.
c    As output point we choose t_out=400 
c
c--------------------------------------------------------
c ----- to integrate with rock2.f ----- 
c     include 'rock2.f'  
c
      include 'rock4.f'       
      implicit double precision (a-h,o-z)
c --- parameters for the problem -----      
       parameter (neqn=400)
c--------------------------------------------------------
c      Work is of lenght 7*neqn because the spectral radius 
c      is computed externally (otherwise should by 8*neqn).
c      If integrating with rock2 define work(4*neqn).
c-------------------------------------------------------- 
c ----- to integrate with rock2.f -----
c     dimension y(neqn),work(4*neqn) 
c  
      dimension y(neqn),work(7*neqn)
      integer iwork(12),idid  
      external finag    
c --- common parameters for the problem -----
      common /par/ alpha,eta,tau,rhodn,dis
      data alpha/ 0.139d0 / 
      data beta / 2.54d0  /
      data eta  / 0.008d0 /
      data rhop  / 0.3d0   /
      data dif  / 0.010d0 /
c ----- file for solution -----
      open(8,file='sol.out')
      rewind 8
c -------- dimensions and initialisations --------
      n=neqn
      hd=200
      dis=(hd*dif)**2
      tau=-eta*beta
      rhodn=rhop/(hd*dif)
c--------------------------------------------------------
c     Initialise iwork: 
c     iwork(1)=1  RHO returns an upper bound for the spectral radius
c     iwork(2)=1  The Jacobian is constant (RHO is called once)
c     iwork(3)=0  Return and solution at tend
c     iwork(4)=0  Atol and rtol are scalar
c--------------------------------------------------------
      iwork(1)=1
      iwork(2)=1
      iwork(3)=0
      iwork(4)=0
c ----- initial values -----
      do i=1,n
        y(i)=0.d0
      end do 
c ----- initial and end point of integration -----
      t=0.0d0
      tend=400.d0
c ----- required tolerance -----
      rtol=0.01d0
      atol=rtol
c ----- initial step size -----
      h=1d-4
c ----- integration -----
      write(6,*) 'Integration of the Finag problem'     
c ----- call of the subroutine rock2 -----
c      call rock2(neqn,t,tend,h,y,finag,atol,rtol,work,
c     &           iwork,idid) 
c
c ----- call of the subroutine rock4 -----
      call rock4(neqn,t,tend,h,y,finag,atol,rtol,work,
     &           iwork,idid)
c ----- print solution -----
      do j=1,400,7
        write (8,*) y(j)
      end do
c ----- print statistics -----
      write(6,*) 'Solution is tabuled in file sol.out'
      write(6,*) 'The value of IDID is',idid
      write(6,*) 'Max estimation of the spectral radius=',iwork(11)
      write(6,*) 'Min estimation of the spectral radius=',iwork(12)
      write(6,*) 'Max number of stages used=',iwork(10)
      write(6,*) 'Number of f eval. for the spectr. radius=',iwork(9)
      write (6,91) iwork(5),iwork(6),iwork(7),iwork(8)
 91   format(' Number of f evaluation=',i5,' steps=',i4,
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
      rho=17.d0
      return
      end
c--------------------------------------------------------
c     The subroutine FINAG compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine finag(neqn,t,y,dy)
      implicit double precision (a-h,o-z)
      dimension y(neqn),dy(neqn)
      common /par/ alpha,eta,tau,rhodn,dis
c ---------- discretisation  -------
      fv=((y(1)-alpha-1.0d0)*y(1)+alpha)*y(1)
      dy(1)=dis*(rhodn-y(1)+y(3))-fv-y(2)
      do i=3,neqn-3,2
        fv=((y(i)-alpha-1.0d0)*y(i)+alpha)*y(i)
        dy(i)=dis*(y(i-2)-2.0d0*y(i)+y(i+2))-fv-y(i+1)
      end do
      fv=((y(neqn-1)-alpha-1.0d0)*y(neqn-1)+alpha)*y(neqn-1)
      dy(neqn-1)=dis*(y(neqn-3)-y(neqn-1))-fv-y(neqn)
c ----
      do i=2,neqn,2
        dy(i)=eta*y(i-1)+tau*y(i)
      end do
      return
      end
      


