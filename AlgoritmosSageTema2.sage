# Metotdos tema 2

def biseccion(f,a,b):
    """
    Ayuda: f funcion, a intervalo inferior, b intervalo superior
    Devuelve nuevos intervalos, m solucion r i num iter
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
        return a,b,m,i
    else:
        print("El intervalo no sirve")

def regula_falsi(f,a,b):
    """
    Ayuda: f funcion, a intervalo inferior, b intervalo superior
    Devuelve nuevos intervalos, m solucion r i num iter
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
        return a,b,m,i
    else:
        print("El intervalo no sirve")

def existePtoFijo (g,a,b):
    """
    Ayuda: g funcion, a intervalo inferior, b intervalo superior
    Devuelve booleano
    """
    c = RealSet.closed(a,b)
    cond1 = g(c.sup()) in c and g(c.inf()) in c
    cond2 = (abs(g(c.sup())) <= 1) and (abs(g(c.inf())) <= 1)
    return cond1 and cond2

def schroeder (g,p):
    """
    Ayuda: obtiene orden de convergencia de g tq g(p)=p
    g funcion, p punto fijo
    Devuelve orden de convergencia y cota de error
    """
    gaux = diff(g(x),x)
    k = 0
    while gaux(p) == 0:
        k = k+1
        gaux = diff(g(x),x)
    
    cota = gaux(p)/factorial(p)
    return k, cota

def condicionesFourier(a,b,f):
    """
    Parámetros ,a,b puntos del intervalo, f funcion
    
    Devuelve: booleano
    """
    df(x) = diff(f(x),x)
    ddf(x) = diff(df(x),x)
    cond = f(a)*f(b) < 0 and df(x) != 0 and sign(ddf(a)) == sign(ddf(b)) 
    return cond

def convergenciaGlobal (a,b,f):
    """
    Parámetros a,b puntos del intervalo, f funcion
    
    Devuelve: booleano
    """
    df(x) = diff(f(x),x)
    auxList = []
    auxList.append(abs(f(b)/df(b)))
    auxList.append(abs(f(a)/df(a)))
    return condicionesFourier(a,b,f) and max(auxList) <= abs(b-a) 

def convergenciaSemilocal (xi,a,b,f):
    """
    Parámetros: x punto a iterar, x pertenece [a,b]
    
    Devuelve booleano
    """
    df(x) = diff(f(x),x)
    ddf(x) = diff(f(x),x)
    return condicionesFourier(a,b,f) and f(x)*ddf(xi) >= 0

def newton(xi,f,it):
    """
    Parámetros: xi punto de salida, f funcion, it iteraciones
    
    Devuelve sol aprox f;
    """
    dif(x) = diff(f(x),x)
    for i in range (it):
        xi = xi - f(xi)/dif(xi)
    
    return N(xi)

def secante (x1,x2,f,it):
    """
    Parametros f funcion x1, x2 puntos iniciales, it iteraciones
    Devuelve sol aproximada
    """
    for i in range (it):
        aux = x2 - f(x1)*((x2-x1)/(f(x2)-f(x1)))
        x1 = x2
        x2 = aux
    
    return N(x2)

def sist2ecuNolineales (f,v,it):
    """
    La funcion tiene que estar definida como f(x,y) = ()
    v es un vector de dos componentes
    """
    x,y = var("x,y")
    jac(x,y) = jacobian(f(x,y),(x,y))
    for i in range (it):
        v = v - f(v[0],v[1])*~(jac(v[0],v[1]))
    
    return N(v)