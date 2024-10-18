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

def LU_dolittle(a):
    p,l,u = A.LU()
    return [p,l,u]
    
def LU_crout(a):
    size = a.nrows()
    p,l,u = A.LU()
    aux = []
    for i in range (size):
        aux.append(u[i,i])
    d = diagonal_matrix(aux)    
    l = l*d
    u = ~d*u
    return [p,l,u]

def jacobi(a):
    size = a.nrows()
    e = []
    f = []
    d = []
    
    for i in range(size):
        d.append(a[i,i])
        
    m = diagonal_matrix(d)
    aux = a-m
    for i in range(size):
        for j in range(size):
            if i<=j:
                e.append(-aux[i,j])
                
    for i in range(size):
        for j in range(size):
            if i>=j:
                f.append(-aux[i,j])
                
    n = upper_triangular_matrix(size,e,SR) + lower_triangular_matrix(size,f,SR)
    
    j = ~m*(n)
    
    return  j

def gauss_seidel(a):
    size = a.nrows()
    e = []
    f = []
    d = []
    
    for i in range(size):
        d.append(a[i,i])
        
    daux = diagonal_matrix(d)
    aux = a-daux
    for i in range(size):
        for j in range(size):
            if i<=j:
                e.append(-aux[i,j])
                
    for i in range(size):
        for j in range(size):
            if i>=j:
                f.append(aux[i,j])
    
    n = upper_triangular_matrix(size,e,SR)
    m = daux + lower_triangular_matrix(size,f,SR)
    
    l = ~m*n
    
    return l

def radio_espectral(a):
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

def sol_jacobi(a,b):
    size = a.nrows()
    j = jacobi(a)
    if radio_espectral(j) < 1:
        daux = []
        for i in range(size):
            daux.append(a[i,i])
        d = ~diagonal_matrix(daux)
        v = vector_nulo(size)
        for i in range(100):
            v = j*v+d*b
        return N(v)
    else:
        print("Radio espectral mayor que uno")

def sol_gauss_seidel(a,b):
    size = a.nrows()
    l = gauss_seidel(a)
    if radio_espectral(l) < 1:
        e = []
        f = []
        d = []
    
        for i in range(size):
            d.append(a[i,i])
        
        daux = diagonal_matrix(d)
        aux = a-daux
        for i in range(size):
            for j in range(size):
                if i<=j:
                    e.append(-aux[i,j])
                
        for i in range(size):
            for j in range(size):
                if i>=j:
                    f.append(aux[i,j])
    
        n = upper_triangular_matrix(size,e,SR)
        m = daux + lower_triangular_matrix(size,f,SR)
        
        v = vector_nulo(size)
        for i in range(100):
            v = l*v+~m*b
        return N(v)
    else:
        print("Radio espectral mayor que uno")        


def f(x): #la funcion la cual se quiere evaluar
    return N(x*sin(x) -1)

def biseccion(a,b):
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
        return [a,b,m,i]
    else:
        print("El intervalo no sirve")
        
def regula_falsi(a,b):
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
        return [a,b,m,i]
    else:
        print("El intervalo no sirve")