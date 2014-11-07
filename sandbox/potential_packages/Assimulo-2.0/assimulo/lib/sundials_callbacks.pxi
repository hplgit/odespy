
cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.f to the Sundials
    right-hand-side function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    #cdef ndarray[realtype, ndim=1, mode='c'] rhs #Used for return from the user function
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef realtype* resptr=(<N_VectorContent_Serial>yvdot.content).data
    cdef int i
    
    if pData.dimSens>0: #Sensitivity activated
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                rhs = (<object>pData.RHS)(t,y,<list>pData.sw, p=p)
            else:
                rhs = (<object>pData.RHS)(t,y,p)
                
            #memcpy((<N_VectorContent_Serial>yvdot.content).data,<realtype*>rhs.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = rhs[i]
            
            return CV_SUCCESS
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
        
    else: #No sensitivity
        try:
            if pData.sw != NULL:
                rhs = (<object>pData.RHS)(t,y,<list>pData.sw)
            else:
                rhs = (<object>pData.RHS)(t,y)
            
            #memcpy((<N_VectorContent_Serial>yvdot.content).data,<realtype*>rhs.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = rhs[i]
            
            return CV_SUCCESS
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)

cdef int cv_jac(int Neq, realtype t, N_Vector yv, N_Vector fy, DlsMat Jacobian, 
                void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    #cdef ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
    cdef realtype* col_i=DENSE_COL(Jacobian,0)
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef int i,j
    
    try:
        if pData.sw != NULL:
            #jac=(<object>pData.JAC)(t,(<ndarray>pData.y),<list>pData.sw)
            jac=(<object>pData.JAC)(t,y,<list>pData.sw)
        else:
            jac=(<object>pData.JAC)(t,y)
        
        #This needs further investigations:
        #memcpy(Jacobian.data,<realtype*>jac.data, pData.memSizeJac)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jacobian, i)
            for j in range(Neq):
                col_i[j] = jac[j,i]

        return CVDLS_SUCCESS
    except:
        return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        
cdef int cv_jacv(N_Vector vv, N_Vector Jv, realtype t, N_Vector yv, N_Vector fyv,
				    void *problem_data, N_Vector tmp):
    """
    This method is used to connect the Assimulo.Problem.jacv to the Sundials
    Jacobian times vector function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y  = nv2arr(yv)
    cdef N.ndarray v  = nv2arr(vv)
    cdef N.ndarray fy = nv2arr(fyv)
    cdef int i
    
    cdef realtype* jacvptr=(<N_VectorContent_Serial>Jv.content).data
    
    try:
        jacv = (<object>pData.JACV)(t,y,fy,v)
        
        for i in range(pData.dim):
                jacvptr[i] = jacv[i]
        
        return SPGMR_SUCCESS
    except:
        return SPGMR_ATIMES_FAIL_REC
    

cdef int cv_root(realtype t, N_Vector yv, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    Root-finding function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    #cdef ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef int i
    
    try:
        if pData.sw != NULL:
            root=(<object>pData.ROOT)(t,y,<list>pData.sw) #Call to the Python root function 
        else:
            root=(<object>pData.ROOT)(t,y,None) #Call to the Python root function
            
        #memcpy(gout,<realtype*>root.data,pData.memSizeRoot) #Copy data from the return to the output
        for i in range(pData.dimRoot):
            gout[i]=root[i]
    
        return CV_SUCCESS
    except:
        return CV_RTFUNC_FAIL  # Unrecoverable Error           

cdef int ida_res(realtype t, N_Vector yv, N_Vector yvdot, N_Vector residual, void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.f to the Sundials
    residual function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray[realtype, ndim=1, mode='c'] res #Used for return from the user function
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef N.ndarray yd = nv2arr(yvdot)
    cdef realtype* resptr=(<N_VectorContent_Serial>residual.content).data
    cdef int i
    
    if pData.dimSens!=0: #SENSITIVITY 
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                res=(<object>pData.RHS)(t,y,yd,sw=<list>pData.sw,p=p)  # call to the python residual function
            else:
                res=(<object>pData.RHS)(t,y,yd,p)
            
            #memcpy((<N_VectorContent_Serial>residual.content).data,<realtype*>res.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = res[i]

            return IDA_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            print "Unexpected error, probably due to a programing error in rhs/res function:\n"
            traceback.print_exc()
            return -1
    else: #NO SENSITIVITY
        try:
            if pData.sw != NULL:
                res=(<object>pData.RHS)(t,y,yd,<list>pData.sw)  #Call to the Python residual function
            else:
                res=(<object>pData.RHS)(t,y,yd)
                #res = (<object>pData.RHS)(t,y,yd)
            
            #memcpy((<N_VectorContent_Serial>residual.content).data,<realtype*>res.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = res[i]
            
            return IDA_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            print "Unexpected error, probably due to a programing error in rhs/res function:\n"
            traceback.print_exc()
            return -1
            
cdef int ida_jac(int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, DlsMat Jacobian,
                 void* problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
    cdef realtype* col_i=DENSE_COL(Jacobian,0)
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef N.ndarray yd = nv2arr(yvdot)
    cdef int i,j
    
    try:
        if pData.sw != NULL:
            jac=(<object>pData.JAC)(c,t,y,yd,<list>pData.sw)  # call to the python residual function
        else:
            jac=(<object>pData.JAC)(c,t,y,yd)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jacobian, i)
            for j in range(Neq):
                col_i[j] = jac[j,i]
        return IDADLS_SUCCESS
    except: 
        return IDADLS_JACFUNC_RECVR #Recoverable Error

cdef int ida_root(realtype t, N_Vector yv, N_Vector yvdot, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    root function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef N.ndarray yd = nv2arr(yvdot)
    cdef int i
    
    try:
        if pData.sw != NULL:
            root=(<object>pData.ROOT)(t,y,yd,<list>pData.sw)  #Call to the Python root function
        else:
            root=(<object>pData.ROOT)(t,y,yd,None)  #Call to the Python root function
    
        #memcpy(gout,<realtype*>root.data,pData.memSizeRoot) #Copy data from the return to the output
        for i in range(pData.dimRoot):
            gout[i]=root[i]
        
        return IDA_SUCCESS
    except:
        return IDA_RTFUNC_FAIL  # Unrecoverable Error
    

# Error handling callback functions
# =================================

cdef int cv_err(int error_code, char *module, char *function, char *msg, void *problem_data):
    """
    This method overrides the default handling of error messages.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    
    if error_code > 0 and pData.verbose > 0: #Warning
        print '[CVode Warning]', msg
    
    if pData.verbose > 2: #Verbosity is greater than NORMAL, print warnings and errors
        if error_code < 0: #Error
            print '[CVode Error]', msg
            
cdef int ida_err(int error_code, char *module, char *function, char *msg, void *problem_data):
    """
    This method overrides the default handling of error messages.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    
    if error_code > 0 and pData.verbose > 0: #Warning
        print '[IDA Warning]', msg
    
    if pData.verbose > 2: #Verbosity is greater than NORMAL, print warnings and errors
        if error_code < 0: #Error
            print '[IDA Error]', msg



cdef class ProblemData:
    cdef:
        void *RHS          #Should store the residual or the right-hand-side
        void *ROOT         #Should store the root function
        void *JAC          #Should store the jacobian
        void *JACV         #Should store the jacobian times a vector
        void *SENS         #Should store the sensitivity function
        void *y            #Temporary storage for the states
        void *yd           #Temporary storage for the derivatives
        void *sw           #Storage for the switches
        realtype *p            #Storage for the parameters
        realtype *pbar
        int dim            #Dimension of the problem
        int dimRoot        #Dimension of the roots
        int dimSens        #Dimension of the parameters (For sensitivity)
        int memSize        #dim*sizeof(realtype) used when copying memory
        int memSizeRoot    #dimRoot*sizeof(realtype) used when copying memory
        int memSizeJac     #dim*dim*sizeof(realtype) used when copying memory
        int verbose        #Defines the verbosity
        
#=================
# Module functions
#=================

cdef inline N_Vector arr2nv(x):
    x=N.array(x)
    cdef long int n = len(x)
    cdef N.ndarray[realtype, ndim=1,mode='c'] ndx=x
    cdef void* data_ptr=PyArray_DATA(ndx)
    cdef N_Vector v=N_VNew_Serial(n)
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(realtype))
    return v
    
cdef inline N.ndarray nv2arr(N_Vector v):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    cdef N.ndarray[realtype, ndim=1, mode='c'] x=N.empty(n)
    memcpy(x.data, v_data, n*sizeof(realtype))
    return x

cdef inline realtype2arr(realtype *data, int n):
    """Create new numpy array from realtype*"""
    cdef N.ndarray[realtype, ndim=1, mode='c'] x=N.empty(n)
    memcpy(x.data, data, n*sizeof(realtype))
    return x
