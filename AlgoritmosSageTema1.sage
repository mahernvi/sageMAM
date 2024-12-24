#Metotdos tema 1

import re

def upper_triangular_matrix(size, lista, field=SR):
    m = zero_matrix(field,size,size)
    for i in range(size):
        for j in range(size):
            if i<=j:
                m[i,j] = lista[j - i/2 * (1 + i - 2*size)]
    return m

def lower_triangular_matrix(size, lista, field=SR):
    m = zero_matrix(field,size,size)
    for i in range(size):
        for j in range(size):
            if i>=j:
                m[i,j] = lista[i/2*(i + 1) + j]
    return m

def diagonal_superior(val,size,field):
    m = zero_matrix(field,size,size)
    i = 1
    while i < size:
        for j in range (size-1):
            m[j,i] = val
            i = i+1
    return m

def diagonal_inferior(val,size,field):
    m = zero_matrix(field,size,size)
    i = 1
    while i < size:
        for j in range (size-1):
            m[i,j] = val
            i = i+1
    return m

def num_cond (a,n):
    """Ayuda:
    Parametros: matriz A, numero de norma 1,2, oo; n
    
    Devuelve A.nomr(n)*(~A).norm(n)
    """
    v1 = a.norm(n)
    v2 = (~a).norm(n)
    return v1*v2

def LU_dolittle(a):
    """Ayuda:
    Parametros:
        Matriz
    
    Devuelve:
        | Matrices P, L, U
        | P*A == L*U
    """
    p,l,u = a.LU()
    return ~p,l,u
    
def LU_crout(a):
    """Ayuda:
    Parametros:
        Matriz
    
    Devuelve:
        | Matrices P, L, U
        | P*A == L*U
    """
    size = a.nrows()
    p,l,u = LU_dolittle(a)
    aux = []
    for i in range (size):
        aux.append(u[i,i])
    d = diagonal_matrix(aux)    
    l = l*d
    u = ~d*u
    return p,l,u

def separar_matriz(a):
    """Ayuda:
    Parametros:
        Matriz
    
    Devuelve:
        | D = diagonal de A
        | E = triangular inferio A
        | F = triangular superior A
    """
    if not a.is_square():
        raise ValueError("La matriz no es cuadrada")
    size = a.nrows()
    daux = []
    for i in range(size):
        daux.append(a[i,i])
    d = diagonal_matrix(daux)
    eaux = []
    faux = []
    aux = a - d
    for i in range(size):
        for j in range(size):
            if i<=j:
                faux.append(-aux[i,j])
    for i in range(size):
        for j in range(size):
            if i>=j:
                eaux.append(-aux[i,j])
    e = lower_triangular_matrix(size,eaux,SR)
    f = upper_triangular_matrix(size,faux,SR)
    return d, e, f

def jacobi(a):
    """Ayuda:
    Parametros: 
        Matriz
    
    Devuelve:
        Matriz de jacobi
    """
    d,e,f = separar_matriz(a)
    m = d
    n = e+f
    j = ~m*(n)
    
    return  j

def gauss_seidel(a):
    """Ayuda:
    Parametros: 
        Matriz
    
    Devuelve:
        Matriz de Gauss-Seidel
    """
    d,e,f = separar_matriz(a)
    m = (d-e)
    n = f
    l = ~m*n
    
    return l

def SOR(a,w):
    """Ayuda:
    Parametros: 
        | Matriz, parametro w
        | w = 1 es Gauss-Seidel
    
    Devuelve:
        Matriz SOR en funcion de w
    """
    d,e,f = separar_matriz(a)
    m = d -w*e
    n = (1-w)*d + w*f
    l = ~m*n
    
    return l

def valor_optimo_SOR(a):
    """Ayuda:
    Parametros:
        Matriz
    Devuelve:
        Valor optimo de w
    """
    l = gauss_seidel(a)
    k = radio_espectral(l)
    return 2/(1+sqrt(1-k))

def radio_espectral(a):
    """Ayuda:
    Parametros:
        | Matriz
    
    Devuelve:
        Radio espectral de la matriz
    """
    size = a.nrows()
    vaps = a.eigenvalues()
    absVaps = []
    for i in range(size):
        absVaps.append(abs(vaps[i]))
        
    re = N(max(absVaps))
    
    return re

def vector_nulo (size):
    v = []
    for i in range(size):
        v.append(0)
    
    return vector(v)

def sol_jacobi(a,b,iteraciones):
    """Ayuda:
    Parametros:
        Matriz A, vector b, num iteraciones
        
    Devuelve:
        | Solución aproximada mediante jacobi
        | Si el radio espectral es mayor que uno
        | no lo hace
    """
    size = a.nrows()
    j = jacobi(a)
    if radio_espectral(j) > 1:
        print("Radio espectral mayor que uno")
    else:
        d,e,f = separar_matriz(a)
        v = vector_nulo(size)
        for i in range(iteraciones):
            v = j*v + ~d*b
        return N(v)

def sol_gauss_seidel(a,b,iteraciones):
    """Ayuda:
    Parametros:
        Matriz A, vector b, num iteraciones
        
    Devuelve:
        | Solución aproximada mediante gauss_seidel
        | Si el radio espectral es mayor que uno
        | no lo hace
    """
    size = a.nrows()
    l = gauss_seidel(a)
    if radio_espectral(l) > 1:
        print("Radio espectral mayor que uno")   
    else:
        d,e,f = separar_matriz(a)
        v = vector_nulo(size)
        for i in range(iteraciones):
            v = l*v + ~(d-e)*b
        return N(v)
             
def sol_SOR(a,b,w,iteraciones):
    """Ayuda:
    Parametros:
        Matriz A, vector b, parametro w, num iteraciones
        
    Devuelve:
        | Solución aproximada mediante SOR
        | Si el radio espectral es mayor que uno
        | no lo hace
    """
    size = a.nrows()
    l = SOR(a,w)
    if radio_espectral(l) > 1:
        print("Radio espectral mayor que uno") 
    else:
        d,e,f = separar_matriz(a)
        v = vector_nulo(size)
        for i in range (iteraciones):
            v = l*v+~(d-w*e)*b
        return N(v)




