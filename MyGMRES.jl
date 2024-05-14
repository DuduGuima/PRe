using LinearAlgebra

function Arnold(A,n,b=A[:,1])
    Q = zeros(size(A)[2],n)
    H = zeros(n,n)
    v=b/norm(b)
    Q[:,1] = v
    for k=1:n
        if k>1
            Q[:,k] = v / H[k,k-1]
        end
        v = A*Q[:,k]
        for j=1:k
            H[j,k] = (Q[:,j]')*v
            v-= H[j,k]*Q[:,j]
        end
        if k<n
            H[k+1,k] = norm(v)
        end
        
    end
    return Q,H
end

function MyGMRES(A,b,maxiter=size(A,2))
    Q=zeros(size(A,1),maxiter)
    H=zeros(maxiter,maxiter)
    x=zeros(size(b))
    bheta = norm(b)
    v = b/bheta
    Q[:,1] = v
    for k=1:maxiter
        #Q,H=Arnold(A,k,b)#melhorar essa insercao de matrizes
        ###Arnold's iteration inside GMRES to use Q,H from past iterations
        if k>1
            Q[:,k] = v/H[k,k-1]
        end
        v = A*Q[:,k]
        for j=1:k
            H[j,k] = (Q[:,j]')*v
            v-= H[j,k]*Q[:,j]
        end
        if k<maxiter
            H[k+1,k]=norm(v)
        end

        e1=Matrix(I,k,1)
        y = H[1:k,1:k]\(bheta*e1)        
        x=Q[:,1:k]*y
        
        if norm(A*x - b) < 0.001#ta zoado esse calculo aqui tb
            return x
        end
    end
    return x
end