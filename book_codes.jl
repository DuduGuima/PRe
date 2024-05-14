
include("./MyGMRES.jl")

function Arnoldi_book(A,q,kmax,b)
    H=zeros(kmax,kmax)
    r=copy(b)
    q[:,1] = r/norm(r)
    for k=1:kmax
        if k>1
            q[:,k] = r/H[k,k-1]
        end
        r=A*q[:,k]
        for i=1:k
            H[i,k]=dot(q[:,i],r)
            r-=H[i,k] *q[:,i]
        end

        if k<kmax
            H[k+1,k] = norm(r)
        end
    end
    return H
end

# function gmres_book(A,b,maxiter)
#     for i=1:maxiter
#         z=M\(A*q[i])
#         for k=1:i 
#             H[k,i]=z' *q[k]
#             z -= H[k,i] *q[k]
#         end
#         H[i+1,i]=norm(z)
#         q[i+1]=z/H[i+1,i]

#         for k=1:i-1
#             H[1:i+1,i]=J[k] * H[1:i+1,i]
#         end
#         J[i], = givens(H[i,i],H[i+1,i],i,i+1)
#         H[1:i+1,i]=J[i]*H[1:i+1,i]
#         s=J[i]*s

#         it+=1

#         res=abs(s[i+1])
#         if res<0.001
#             x=gmres_update(x,s,v,i,H)
#             return x
#         end
#     end
# end