import base_modules
from base_modules import *




def get_payoff_diff(i, N, m):
    
    u = ((m[0,0]+m[1,1]-(m[0,1]+m[1,0]))/(N-1))
    v = ((N*(m[0,1]-m[1,1]) - (m[0,0]-m[1,1]))/(N-1))
    
    delta_pi = (u*i) + v
    
    return delta_pi
    


def get_T_plus(i, N, m, b) :
    
    delta_pi = get_payoff_diff(i, N, m)

    p_plus = 1.0 / (1 + np.exp(-1 * b * delta_pi))
    
    Tplus = ((i * (N-i))/(N*N))*p_plus
    
    return Tplus



def get_T_minus(i, N, m, b) :
    
    u = ((m[0,0]+m[1,1]-(m[0,1]+m[1,0]))/(N-1))
    v = ((N*(m[0,1]-m[1,1]) - (m[0,0]-m[1,1]))/(N-1))
    
    delta_pi = (u*i) + v
    
    p_minus = 1.0 / (1 + np.exp(1 * b * delta_pi))
    
    Tminus = ((i * (N-i))/(N*N))*p_minus

    
    return Tminus


def H_n(N) :
    
    harmonic = 1.00
    
    for i in range(2, N + 1) :
        harmonic += 1 / i
    
    return harmonic



def weak_fixation_time(N, m, b) : 
    
    return 1 + (0.5*(((((N*(m[0,1]-m[1,1]))-m[0,0]+m[1,1])/(N-1))*(((N-1)/(H_n(N-1)))-1))*b))     

def weak_fixation_time_A(N, b, m) : 

    return 1 - (((m[0,0]-m[0,1]-m[1,0]+m[1,1])/(N-1))*(((N*N)+N-6)/36)*b)


def get_gamma(i, N, b, delta_pi):
    
    gamma = np.exp(-b*delta_pi)
    
    return gamma
    


def get_denom(i, N, b, m):
    
    sum_ = 0
    for k in range(1, N):
        prod = 1
        for l in range(1,k+1):
            prod *= get_gamma(l, N, b, get_payoff_diff(l, N, m))
            
        sum_ += prod
        
    denom = 1 + sum_
    
    return denom


def get_num(i, N, b, m):
    
    sum_ = 0
    for k in range(1, i):
        prod = 1
        for l in range(1,k+1):
            prod *= get_gamma(l, N, b, get_payoff_diff(l, N, m))
            
        sum_ += prod
        
    num = 1 + sum_
    
    return num
    
    
            
            
def get_phi_i(i, N, b, m):
    
    phi_i = get_num(i, N, b, m) / get_denom(i, N, b, m)
    
    return phi_i



def get_theory_fixation_time(N, b, m):
    
    u = ((m[0,0]+m[1,1]-(m[0,1]+m[1,0]))/(N-1))
    v = ((N*(m[0,1]-m[1,1]) - (m[0,0]-m[1,1]))/(N-1))
    
    sum2 = 0
    for k in range(1, N):
        sum1 = 0
        for l in range(1, k+1):
            g = 0 
            for w in range(l+1, k+1):
                g += get_payoff_diff(w, N, m)
                
            f = ((N*N)/(l*(N-l))) * (1 + np.exp(-b*((u*l)+v))) * np.exp(-b*g)
            sum1 += f
        
        sum2 += sum1
        
    theory_fixation_time = get_phi_i(1, N, b, m) * sum2
    
    return theory_fixation_time


def get_theory_fixation_time_A(N, b, m):
    
    u = ((m[0,0]+m[1,1]-(m[0,1]+m[1,0]))/(N-1))
    v = ((N*(m[0,1]-m[1,1]) - (m[0,0]-m[1,1]))/(N-1))
    
    sum2 = 0
    for k in range(1, N):
        sum1 = 0
        for l in range(1, k+1):
            g = 0 
            for w in range(l+1, k+1):
                g += get_payoff_diff(w, N, m)
                
            f = (get_phi_i(l,N,b,m)) * ((N*N)/(l*(N-l))) * (1 + np.exp(-b*((u*l)+v))) * np.exp(-b*g)
            sum1 += f
        
        sum2 += sum1

        
        theory_fixation_time_A = sum2
    
    return theory_fixation_time_A



def get_theory_tau(N, b, m):
    
    theory_tau = get_theory_fixation_time(N, b, m)/get_theory_fixation_time(N, 0, m)
    
    return theory_tau


def get_theory_tau_A(N, b, m):

    theory_tau_A = get_theory_fixation_time_A(N, b, m)/get_theory_fixation_time_A(N, 0, m)

    return theory_tau_A






        
        