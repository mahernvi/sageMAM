def vector_inicial (size):
    """
    Ayuda: size dim vector
    """
    vaux = []
    for i in range (size):
        vaux.append(1)
    v = vector(vaux)
    return v

def phi (v,size):
    """
    Función auxiliar para el metodod de la potencia
    Ayuda: v vector,size tamaño del vector
    """
    aux = 0;
    for i in range (size):
        aux = aux + v[i]
    
    return aux;

def metodoPotencia (a,vx,size,it):
    """
    Ayuda: a matriz,vx vector inicial ,size dim, it iteraciones
    Devuelve: valor propio dominate aproximado y vector propio correspondiente
    """
    for i in range (it):
        vy = a*vx
        r = phi(vy,size)/phi(vx,size)
        vx = vy/vy.norm(oo)
    
    return N(r),N(vx)

def metodoPotenciaInversa (a,vx,size,it):
    """
    Ayuda: a matriz,vx vector inicial ,size dim, it iteraciones
    Devuelve: valor propio más pequeño aproximado y vector propio correspondiente
    """
    
    for i in range (it):
        vy = ~a*vx
        r = phi(vy,size)/phi(vx,size)
        vx = vy/vy.norm(oo)
        
    return N(r),N(vx)

def metodoPotenciaDesplazamiento (a,vx,desplazamiento,size,it):
    """
    Ayuda: a matriz,vx vector inicial ,size dim, it iteraciones
    Devuelve: valor propio dominante más cercano al desplazamiento y vector propio correspondiente
    """
    for i in range (it):
        vy = (a-desplazamiento*identity_matrix(size))*vx
        r = phi(vy,size)/phi(vx,size)
        vx = vy/vy.norm(oo)
        
    return N(r),N(vx)

def metodoPotenciaInversaDesplazamiento (a,vx,desplazamiento,size,it):
    """
    Ayuda: a matriz,vx vector inicial ,size dim, it iteraciones
    Devuelve: valor propio más pequeño y más cercano al desplazamiento y vector propio correspondiente
    """
    for i in range (it):
        vy = ~(a-desplazamiento*identity_matrix(size))*vx
        r = phi(vy,size)/phi(vx,size)
        vx = vy/vy.norm(oo)
        
    return N(r),N(vx)

