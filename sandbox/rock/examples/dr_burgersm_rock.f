c * * * * * * * * * * * * * * * * * * * * * * * * *
c    Driver for ROCK4 (or ROCK2) at Burger problem
c * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    This driver shows how to use ROCK4 (or ROCK2). It solves a
c    system of ODEs resulting from discretization of the 
c    Burger equation (u=u(x,y)):
c     
c    u_t+(u^2/2)_x=m*u_{xx}  for 0<= t <=2.5, 0<= x <= 1
c                            where m=0.0003
c
c    with initial condition 
c
c    u(x,0)=1.5*x*(1-x)^2
c
c    and  boundary conditions
c
c    u(0,t)=u(1,t)=0 
c
c
c    We discretize the space variables with
c    x_i=i/(N+1) for i=0,1,...,N, with N=500.
c    We obtain a system of 500 equations. 
c    The spectral radius of the Jacobian can be  
c    estimated with the Gerhgorin theorem. Thus  
c    we provide an external function RHO, 
c    giving the spectral radius of the Jacobian 
c    matrix. As output point we choose t_out=2.5
c 
c
c--------------------------------------------------------
c ----- to integrate with rock2.f ----- 
c     include 'rock2.f'  
c
      include 'rock4.f'     
      implicit double precision (a-h,o-z)
c --- parameters for the problem -----
      parameter (neqn=500)
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
      external fburg
c --- common parameters for the problem -----
      common/deltax/deltax,amu,tau
c ----- file for solution -----
       open(8,file='sol.out')
       rewind 8
c -------- dimensions and initialisations --------
        n=neqn
        anp1=n+1
        deltax=1.0d0/anp1
        amu=0.0003d0
        tau=amu/(deltax*deltax)
c--------------------------------------------------------
c     Initialise iwork: 
c      iwork(1)=1  RHO returns an upper bound for the spectral radius
c      iwork(2)=1  The Jacobian is constant (RHO is called once)
c      iwork(3)=0  Return and solution at tend
c      iwork(4)=0  Atol and rtol are scalar
c--------------------------------------------------------
      iwork(1)=1
      iwork(2)=1
      iwork(3)=0
      iwork(4)=0
c ----- initial and end point of integration -----
      t=0.0d0
      tend=2.5d0
c --------- initial values -------
      do i=1,n
        x=i*deltax
        y(i)=1.5d0*x*(1.d0-x)**2
      end do
c ----- required tolerance -----
      rtol=0.01d0
      atol=rtol
c ----- initial step size -----
      h=1.0d-5
c ----- integration -----
      write(6,*) 'Integration of the Burger problem'     
c ----- to integrate with rock2.f ----- 
c      call rock2(neqn,t,tend,h,y,fburg,atol,rtol,work,
c     &           iwork,idid) 
c
c ----- call of the subroutine rock4 -----
       call rock4(neqn,t,tend,h,y,fburg,atol,rtol,work,
     &           iwork,idid)
c ----- print solution -----
      do j=1,498,7
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
      common/deltax/deltax,amu,tau
      rho=0.5d0*(1.d0/deltax)+4.d0*tau+2.d0
      return
      end
c--------------------------------------------------------
c     The subroutine FBURG compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine fburg(neqn,x,y,f)
      implicit double precision (a-h,o-z)
      dimension y(neqn),f(neqn)
      common/deltax/deltax,amu,tau
c ---------- discretisation  -------
c ---------- i=1 --------------
       urig=y(2)
       f(1)=(-urig*urig)/(4.d0*deltax)+
     &       tau*(-2.d0*y(1)+urig)
c ----------- i=2, neqn-1 -----------
       do  i=2,neqn-1
         ulef=y(i-1)
         urig=y(i+1)
         f(i)=(ulef*ulef-urig*urig)/(4.d0*deltax)+
     &         tau*(ulef-2.d0*y(i)+urig)
       end do
c ---------- i=neqn ---------
      ulef=y(neqn-1)
        f(neqn)=(ulef*ulef)/(4.d0*deltax)+
     &         tau*(ulef-2.d0*y(neqn))
      return
      end
c