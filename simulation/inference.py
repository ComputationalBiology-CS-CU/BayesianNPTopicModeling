from scipy.special import digamma as digamma
from scipy.special import beta as beta
from scipy.misc import factorial as fact
import numpy as np
import numpy.random as rd

def stirling_approx(x):
    y = np.copy(x)
    small_index = y < 2
    y[small_index] = 0
    y[np.invert(small_index)] = y[np.invert(small_index)]*np.log(y[np.invert(small_index)]) - y[np.invert(small_index)] + 0.5*np.log(2*math.pi*y[np.invert(small_index)])
    return y

def nonzero_log(x):
    y = np.copy(x)
    zeros = y < 0.0001
    y[zeros] = 0
    y[np.invert(zeros)] = np.log(y[np.invert(zeros)])
    return y

def nonoverflow_softmax(x):
    max_x = np.amax(x,axis=2)[:,:,None]
    return np.exp(x-max_x)/np.sum(np.exp(x-max_x), axis=2)[:,:,None]

#size of plates
M = 100000
N = 2000
T = 17
S = 200
K = 500
J = 50

#hyperparameters
alpha = 1
omega = 1
gamma = 1
eta = [1,1,1]


#initialization
#dimension preassignment
v_alpha = np.zeros((S,J,2)) #beta, sj
v_omega = np.zeros((K,2)) #beta, k
v_gamma = np.zeros((K,N,T,2)) #gamma, knt
v_eta = np.ones((K,M,3)) #3-dirichlet, km
phiE = np.zeros((S,N,J)) #J-categorical, sn
phiO = np.zeros((S,M,J)) #J-categorical, sm
v_nu = np.zeros((S,J,K)) #K-categorical, sj

init_random = True
if init_random:
    #initialize variational parameters randomly
    
    for k in range(K):
        if (k+1)%10 == 0:
            print '%s:%d out of %d' %('k',k+1,K)
        v_omega[k,:] = rd.uniform(0.01,10,2)
        for n in range(N):
            for t in range(T):
                v_gamma[k,n,t,:] = rd.uniform(0.01,10,2)
        
    for s in range(S):
        if (s+1)%10 == 0:
            print '%s:%d out of %d' %('s',s+1,S)
        for j in range(J):
            v_alpha[s,j,:] = rd.uniform(0.01,10,2)
            v_nu[s,j,:] = rd.dirichlet(np.ones(K))
        for n in range(N):
            phiE[s,n,:] = rd.dirichlet(np.ones(J))
        for m in range(M):
            phiO[s,m,:] = rd.dirichlet(np.ones(J))  

#reading data
folder_path = './data/'

#beta = np.zeros((K,M,3))
E = np.zeros((S,N,T))
#lambd = np.zeros((K,N,T))
#for i in range(3):
#    beta[:,:,i] = np.loadtxt(folder_path+'beta'+str(i)+'.csv', delimiter='\t')
    
for t in range(T):
    E[:,:,t] = np.loadtxt(folder_path+'E'+str(t)+'.csv', delimiter='\t')
#    lambd[:,:,t] = np.loadtxt(folder_path+'lambda'+str(t)+'.csv', delimiter='\t')
    
#c = np.loadtxt(folder_path+'c.csv', delimiter='\t')
#v = np.loadtxt(folder_path+'v.csv', delimiter='\t')
O = np.loadtxt(folder_path+'O.csv', delimiter='\t')
#pi = np.loadtxt(folder_path+'pi.csv', delimiter='\t')
#zE = np.loadtxt(folder_path+'zE.csv', delimiter='\t')
#zO = np.loadtxt(folder_path+'zO.csv', delimiter='\t')

#O with one-hot encoding
O_onehot = np.zeros((S,M,3))
for l in range(3):
    O_onehot[:,:,l] = (O == l).astype(int)
    
alpha_normalizer = S*J*np.log(beta(1,alpha))
log_fact_E = np.sum(stirling_approx(E), axis=2)

#updates
alpha_normalizer = S*J*np.log(beta(1,alpha))
log_fact_E = np.sum(stirling_approx(E), axis=2)

conv_criterion = 0.001
max_step = 600
step = 0;
prev_ELBO = -1000000
while True:
    print '%s: %d' %('step', step)
    
    #v_alpha
    sum_phiO_m = np.sum(phiO, axis=1)
    sum_phiE_n = np.sum(phiE, axis=1)
    v_alpha[:,:,0] = 1 + sum_phiO_m + sum_phiE_n

    sum_phiO_m_j = np.sum(sum_phiO_m, axis=1)
    sum_phiE_n_j = np.sum(sum_phiE_n, axis=1)
    for j in range(J-1):
        sum_phiO_m_j = sum_phiO_m_j - sum_phiO_m[:,j]
        sum_phiE_n_j = sum_phiE_n_j - sum_phiE_n[:,j]
        v_alpha[:,j,1] = alpha + sum_phiO_m_j + sum_phiE_n_j
    v_alpha[:,J-1,1] = np.zeros(S) + alpha
    
    #v_omega
    sum_nu_sj = np.sum(np.sum(v_nu, axis=0), axis=0)
    v_omega[:,0] = 1 + sum_nu_sj
    
    sum_nu_sj_k = np.sum(sum_nu_sj, axis=0)
    for k in range(K-1):
        sum_nu_sj_k = sum_nu_sj_k - sum_nu_sj[k]
        v_omega[k,1] = omega + sum_nu_sj_k
    v_omega[K-1,1] = omega
    
    #v_gamma
    prod_nu_phiE_j = np.zeros((K,N,S))
    for s in range(S):
        prod_nu_phiE_j[:,:,s] = np.dot(phiE[s,:,:], v_nu[s,:,:]).T
    
    prod_E_nuphiE_s = np.zeros((K,N,T))
    for n in range(N):
        prod_E_nuphiE_s[:,n,:] = np.dot(prod_nu_phiE_j[:,n,:], E[:,n,:])
    
    v_gamma[:,:,:,0] = gamma + prod_E_nuphiE_s
    v_gamma[:,:,:,1] = np.ones((K,N,T)) + np.sum(prod_nu_phiE_j, axis=2)[:,:,None]
    
    #v_eta
    prod_nu_phiO_j = np.zeros((K,M,S))
    for s in range(S):
        prod_nu_phiO_j[:,:,s] = np.dot(phiO[s,:,:], v_nu[s,:,:]).T

    prod_O_nuphiO_s = np.zeros((K,M,3))
    for m in range(M):
        prod_O_nuphiO_s[:,m,:] = np.dot(prod_nu_phiO_j[:,m,:], O_onehot[:,m,:])

    v_eta = eta + prod_O_nuphiO_s
    
    #phiE
    #digamma of v_alpha part (reusable in phiO) (S,J)
    #di_alpha reusable
    di_alpha1 = digamma(v_alpha[:,:,0])
    di_alpha2 = np.hstack((np.zeros((S,1)), digamma(v_alpha[:,:,1])[:,0:J-1])) #discard last column since not used
    di_alpha12 = digamma(np.sum(v_alpha, axis=2))

    di_alpha = di_alpha1
    di_alpha2_sum = np.zeros(S)
    di_alpha12_sum = np.zeros(S)
    for j in range(J):
        di_alpha2_sum = di_alpha2_sum + di_alpha2[:,j]
        di_alpha12_sum = di_alpha12_sum + di_alpha12[:,j]
        di_alpha[:,j] = di_alpha[:,j] + (di_alpha2_sum - di_alpha12_sum)
    
    #sum over t part (reusable in v_nu) (S,K,N)
    #t_sum reusable
    t_sum12 = np.zeros((S,K,N))
    zeros_index_gamma1 = v_gamma[:,:,:,1] < 0.0001
    nonzero_gamma1 = v_gamma[:,:,:,1]
    nonzero_gamma1[zeros_index_gamma1] = 1
    for n in range(N):
        di_log_gamma = (digamma(v_gamma[:,n,:,0]) - nonzero_log(v_gamma[:,n,:,1])).T # (T,K)
        t_sum12[:,:,n] = np.dot(E[:,n,:], di_log_gamma) # (S,K)
    t_sum3 = np.sum(v_gamma[:,:,:,0]/nonzero_gamma1, axis=2)
    t_sum = t_sum12 - t_sum3
    t_sum[zeros_index_gamma1] = 0
    
    #sum over k part (S,N,J)
    k_t_sum = np.zeros((S,N,J))
    for s in range(S):
        k_t_sum[s,:,:] = np.dot(v_nu[s,:,:], t_sum[s,:,:]).T
    
    phiE = nonoverflow_softmax(k_t_sum + di_alpha[:,None,:]) #normalize
    
    #phiO
    #sum over l part (reusable in v_nu) (S,K,M)
    #l_sum reusable
    di_eta = digamma(v_eta) - digamma(np.sum(v_eta, axis=2))[:,:,None]
    l_sum = np.zeros((S,K,M))
    for m in range(M):
        l_sum[:,:,m] = np.dot(O_onehot[:,m,:], di_eta[:,m,:].T)
    
    #sum over k part (S,M,J)
    k_l_sum = np.zeros((S,M,J))
    for s in range(S):
        k_l_sum[s,:,:] = np.dot(v_nu[s,:,:], l_sum[s,:,:]).T
    
    phiO = nonoverflow_softmax(k_l_sum + di_alpha[:,None,:]) #normalize
    
    #v_nu
    #digamma omega part
    di_omega1 = digamma(v_omega[:,0])
    di_omega2 = np.hstack(([0],digamma(v_omega[:,1])[0:K-1]))
    di_omega12 = digamma(np.sum(v_omega, axis=1))
    
    di_omega = di_omega1 #(1,K)vector
    di_omega2_sum = 0
    di_omega12_sum = 0
    for k in range(K):
        di_omega2_sum = di_omega2_sum + di_omega2[k]
        di_omega12_sum = di_omega12_sum + di_omega12[k]
        di_omega[k] = di_omega[k] + (di_omega2_sum - di_omega12_sum)
        
    #sum over n,t and m,l part
    n_t_sum = np.zeros((S,J,K))
    m_l_sum = np.zeros((S,J,K))
    for s in range(S):
        n_t_sum[s,:,:] = np.dot(t_sum[s,:,:], phiE[s,:,:]).T
        m_l_sum[s,:,:] = np.dot(l_sum[s,:,:], phiO[s,:,:]).T
        
    v_nu = nonoverflow_softmax(di_omega[None,None,:] + n_t_sum + m_l_sum)
    #max_v_nu_k = np.amax(v_nu, axis=2)
    #log_v_nu = v_nu - max_v_nu_k[:,:,None] - np.log(np.sum(np.exp(v_nu - max_v_nu_k[:,:,None]), axis=2))[:,:,None]
    #v_nu = np.exp(log_v_nu) #log-sum-exp trick
    
    #calculate ELBO and check convergence
    ELBO_sj_nu = np.sum((di_omega[None,None,:] - nonzero_log(v_nu)) * v_nu, axis=2) #(S,J)
    ELBO_sj_alpha = digamma(v_alpha[:,:,1]) * (alpha - v_alpha[:,:,1]) - di_alpha1 * (v_alpha[:,:,0] - 1) - di_alpha12 * (np.sum(v_alpha, axis=2) - alpha - 1) + beta(v_alpha[:,:,0], v_alpha[:,:,1])
    ELBO_sj = np.sum(ELBO_sj_nu + ELBO_sj_alpha) - alpha_normalizer
    
    ELBO_sm = np.sum(phiO * (k_l_sum + di_alpha[:,None,:] - nonzero_log(phiO)))
    ELBO_sn = np.sum(phiE * (k_t_sum + di_alpha[:,None,:] - log_fact_E[:,:,None] - nonzero_log(phiE)))
    
    ELBO = ELBO_sj + ELBO_sm + ELBO_sn   
    converg = abs(prev_ELBO - ELBO)/abs(prev_ELBO)
    print 'ELBO: %e, convergence: %e' %(ELBO, converg)
    
    if converg < conv_criterion:
        print 'converged.'
        break
    if step > max_step:
        print 'max step reached.'
        break
    prev_ELBO = ELBO
    step = step + 1

#save variational parameters
np.save(folder_path+"output/"+"v_alpha.dat", v_alpha)
np.save(folder_path+"output/"+"v_omega.dat", v_omega)
np.save(folder_path+"output/"+"v_gamma.dat", v_gamma)
np.save(folder_path+"output/"+"v_eta.dat", v_eta)
np.save(folder_path+"output/"+"phiO.dat", phiO)
np.save(folder_path+"output/"+"phiE.dat", phiE)
np.save(folder_path+"output/"+"v_nu.dat", v_nu)
np.savetxt(folder_path+"output/"+"ELBO.txt", ELBO)