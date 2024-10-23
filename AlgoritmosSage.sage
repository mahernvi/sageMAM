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


def biseccion(f,a,b):
    """Ayuda:
       
    Parametros:
        Necesitas declararte una función f(x)
        
        Ejemplo:
        
        | f(x) = x*sin(x) -1 
        
    Devuelve:
        | Una lista con el intervalo nuevo (a,b)
        | el valor de la solución y el nº de iteraciones
        
    """
    
    if f(a)>0 and f(b)<0 or f(a)<0 and f(b)>0:
        print("El intervalo sirve: ")
        EPSILON = 0.0001
        i = 1
        m = (a+b)/2
        while abs(N(f(m))) > EPSILON:
            if f(a) > 0:
                if f(m) > 0:
                    a = N(m)
                else:
                    b = N(m)
            else:
                if f(m) < 0:
                    a = N(m)
                else:
                    b = N(m)
                m = N((a+b)/2)
                i = i+1
        return [(a,b),m,i]
    else:
        print("El intervalo no sirve")
        
def regula_falsi(f,a,b):
    """Ayuda:
       
    Parametros:
        Necesitas declararte una función f(x)
        
        Ejemplo:
        
        | f(x) = x*sin(x) -1  
        
    Devuelve:
        | Una lista con el intervalo nuevo (a,b)
        | el valor de la solución y el nº de iteraciones
        
    """
    
    if f(a)>0 and f(b)<0 or f(a)<0 and f(b)>0:
        print("El intervalo sirve: ")
        EPSILON = 0.0001
        i = 1
        m = (a*f(b)-b*f(a))/(f(b)-f(a))
        while abs(N(f(m))) > EPSILON:
            if f(a) > 0:
                if f(m) > 0:
                    a = N(m)
                else:
                    b = N(m)
            else:
                if f(m) < 0:
                    a = N(m)
                else:
                    b = N(m)
                m = N((a*f(b)-b*f(a))/(f(b)-f(a)))
                i = i+1
        return [(a,b),m,i]
    else:
        print("El intervalo no sirve")
        
def es_contractiva (g,a,b):
    """Ayuda:
    Parametros: funcion g, extremo a, extremo b
    Ejemplo:
    | g(x) = cos(x), intervalo [a,b]
    | es_contractiva(g,0,1)
    Devuelve:
    | Si los extremos estan en el intervalo devuelve un dibujo
    | Si los extremos no estan en el intervalo no devuelve nada
    """
    c = RealSet([a,b])
    if g(a) in c and g(b) in c:
        print("Los extremos estan en el intervalo")
        fa(x) = g(a)
        fb(x) = g(b)
        Gf=plot(g(x),x,a,b)
        Ga = plot(fa(x),x,a,b,color = 'red')
        Gb = plot(fb(x),x,a,b,color = 'red')
        show(Gf+Ga+Gb)
        
    else:
        print("Los extremos no estan en el intervalo")
        
def posible_puntofijo (g,p):
    """Ayuda:
    Paramteros: funcion g, punto p
    Ejemplo:
    | g(x) = g(x) = 1+x-x^2/4
    | posible_puntofijo(g,2)
    Devuelve:
    | p atractor o repulsor
    """
    dg(x) = diff(g(x),x)
    if abs(dg(p)) < 1:
        print(p); print( "punto fijo atractor")
    else:
        print(p); print( "punto fijo repulsor")