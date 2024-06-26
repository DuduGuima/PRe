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

function MyGMRES(A,b,maxiter=size(A,2),restart = length(b),see_r = false)
    it=0
    bheta = norm(b)
    m=restart
    x=zeros(size(b))
    res = bheta
    while it<maxiter
        Q=zeros(size(A,1),m+1)
        H=zeros(m+1,m)
        J = Vector{Any}(undef,m)
        e1=zeros(m+1,1)
        e1[1]=bheta
        v = b/bheta
        Q[:,1] = v
        for k=1:m
            it+=1
            if see_r
                println("Iteration: ",it," Current residual: ",res)
            end
            ###Arnold's iteration inside GMRES to use Q,H from past iterations
            #----------------------------------------------
            v = A*Q[:,k]
            for j=1:k
                H[j,k] = (Q[:,j]')*v
                v-= H[j,k]*Q[:,j]
            end
            H[k+1,k]=norm(v)
            Q[:,k+1] = v/H[k+1,k]

            #---------------------------#

            ###Givens rotation
            #-----------------------------------
            #We use this so that H is upper triangular, which speeds up the
            #least square problem we'll have to solve afterwards, to find y
            for j=1:k-1
                H[1:k+1,k] = J[j] * H[1:k+1,k]
            end
            #this givens functions isn't very intuitive
            #given two integers k and k+1, it will return an object with 4 numbers:(k,k+1,s,c)
            #where s and c are the numbers needed to rotate a vector X in a way such that will
            #anulate H[k+1,k] and change H[k,k].
            #it seems the object is pre-programmed to do calculations on vectors where k and k+1
            #make sense, so thats why we have to create e1 as a maxiter + 1 column vector 
            J[k], = givens(H[k,k],H[k+1,k],k,k+1)
            #the loop above only fixes the first k elements, to change the last one:
            H[1:k+1,k] = J[k] * H[1:k+1,k]
            #since now H is rotated in the final linear system, we need to change the RHS too:
            e1=J[k]*e1
            #-------------------------------------
            y = H[1:k,1:k]\e1[1:k] #-> ta feio isso aqui paizao
            x=Q[:,1:k]*y#-> isso aqui entao nem se fala
            res=norm(A*x - b)
            if res < 0.001#ta zoado esse calculo aqui tb
                println("Finished at iteration: ",it+1," Final residual: ",res)
                return x
            end
        end
    end
    println("Maximum iteration reached")
    return x
end