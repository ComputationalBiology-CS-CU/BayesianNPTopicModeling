import numpy as np
import numpy.random as rd

#size of plates
M = 100000
N = 2000
T = 17
S = 200
K = 500
J = 50

print 'M=%d, N=%d, T=%d, S=%d, K=%d, J=%d' %(M,N,T,S,K,J)

#hyperparameters
alpha = 6
omega = 15
gamma = 2
eta = [1,1,1]
print 'alpha=%f, omega=%f, gamma=%f, eta=[%d,%d,%d]' %(alpha,omega,gamma,eta[0],eta[1],eta[2])

#initialize latent variables (empty or with zeros)
#stick-breaking
v = []
pi = np.zeros((S,J))
stick_v = []
stick_pi = np.zeros((S,J))
#topic
beta = np.zeros((K,M,3))
lambd = np.zeros((K,N,T))
#topic assignment
c = np.zeros((S,J))
zE = np.zeros((S,N))
zO = np.zeros((S,M))
#observation
E = np.zeros((S,N,T))
O = np.zeros((S,M))

#stick-breaking process
print 'stick-breaking for v......'
for k in range(K):
    v_k = rd.beta(1,omega)
    sb = v_k
    for i in range(len(v)):
        sb *= 1-v[i]
    v.append(v_k)
    stick_v.append(sb)
print 'stick-breaking for pi......'
for s in range(S):
    pi_s = []
    for j in range(J):
        pi_sj = rd.beta(1,alpha)
        sb = pi_sj
        for i in range(len(pi_s)):
            sb *= 1-pi_s[i]
        pi_s.append(pi_sj)
        stick_pi[s,j] = sb
        pi[s,j] = pi_sj

#topic assignment: categorical drawing
print 'categorical drawing......'
for s in range(S):
    if (s+1)%10 == 0:
        print '%s:%d out of %d' %('s',s+1,S)
    for j in range(J):
        c[s,j] = np.asscalar(np.nonzero(rd.multinomial(1,stick_v))[0])

    for n in range(N):
        zE[s,n] = np.asscalar(np.nonzero(rd.multinomial(1,stick_pi[s,:]))[0])

    for m in range(M):
        zO[s,m] = np.asscalar(np.nonzero(rd.multinomial(1,stick_pi[s,:]))[0])

#topic: dirichlet and gamma drawing
print 'topic drawing......'
for k in range(K):
    if (k+1)%10 == 0:
        print '%s:%d out of %d' %('k',k+1,K)
    for m in range(M):
        beta[k,m,:] = rd.dirichlet(eta)
    for n in range(N):
        for t in range(T):
            lambd[k,n,t] = rd.gamma(gamma)
        
#observation drawing
print 'observation drawing......'
for s in range(S):
    if (s+1)%10 == 0:
        print '%s:%d out of %d' %('s',s+1,S)
    for n in range(N):
        for t in range(T):
            E[s,n,t] = rd.poisson(lam=lambd[c[s,zE[s,n]],n,t])
    for m in range(M):
        O[s,m] = np.asscalar(np.nonzero(rd.multinomial(1,beta[c[s,zO[s,m]],m]))[0])
        
#save the sampled data
print 'v...'
np.savetxt("v.csv", v, fmt='%.8e', delimiter="\t")
print 'pi...'
np.savetxt("pi.csv", pi, fmt='%.8e', delimiter="\t")
print 'beta0...'
np.savetxt("beta0.csv", beta[:,:,0], fmt='%.8e', delimiter="\t")
print 'beta1...'
np.savetxt("beta1.csv", beta[:,:,1], fmt='%.8e', delimiter="\t")
print 'beta2...'
np.savetxt("beta2.csv", beta[:,:,2], fmt='%.8e', delimiter="\t")

print 'c...'
np.savetxt("c.csv", c, fmt='%i', delimiter="\t")
print 'zE...'
np.savetxt("zE.csv", zE, fmt='%i', delimiter="\t")
print 'zO...'
np.savetxt("zO.csv", zO, fmt='%i', delimiter="\t")
print 'O...'
np.savetxt("O.csv", O, fmt='%i', delimiter="\t")

print 'lambda, E...'
for t in range(T):
    print t
    np.savetxt("lambda"+str(t)+".csv", lambd[:,:,t], fmt='%.8e', delimiter="\t")
    np.savetxt("E"+str(t)+".csv", E[:,:,t], fmt='%i', delimiter="\t")