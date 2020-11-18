import math
from scipy.special import comb
from scipy.stats import norm

L = 256
T = 256
W = 128
ALPHA = 150.99775853018372
SIGMA_BASE = 9.553147931219351

def p_t_w(t, w, l):
  return comb(t-1, w-1, exact=True) * pow(2, -l)

def p_z_t_w(z, t, w):
  mu = (t-1)*(3*ALPHA+11) + w*(3*ALPHA + 15/2)
  sigma = 3*SIGMA_BASE*math.sqrt(t+w-1)
  return norm.pdf(z, loc=mu, scale=sigma)

def p_z(z):
  ret = 0
  for t in range(L, L-30, -1):
    for w in range(t, 0, -1):
      ret += (p_z_t_w(z, t, w) * p_t_w(t, w, L))
  return ret


# 11.176580436192054
def h_z_t_w():
  ret = 0
  for t in range(L, 0, -1):
    for w in range(t, 0, -1):
      sigma = 3*SIGMA_BASE*math.sqrt(t+w-1)
      ret += (p_t_w(t, w, L) * math.log2(math.sqrt(2*math.pi*math.e*(sigma**2))))
  return ret

def h_z():
  ret = 0
  for z in range(150000, 200000):
    print(z)
    tmp = p_z(z)
    ret += (-tmp * math.log2(tmp))
    print("ret = {0}".format(ret))
  return ret

# 176797.43127558831
def calc_mu():
  ret = 0
  for t in range(L, 0, -1):
    for w in range(t, 0, -1):
      mu = (t-1)*(3*ALPHA+11) + w*(3*ALPHA + 15/2)
      ret += (mu * p_t_w(t, w, L))
      print(ret)
  return ret

# 560.1436784349372
def calc_sigma():
  ret = 0
  for t in range(L, 0, -1):
    for w in range(t, 0, -1):
      sigma = 3*SIGMA_BASE*math.sqrt(t+w-1)
      ret += ((sigma**2) * p_t_w(t, w, L))
      print(math.sqrt(ret))
  return math.sqrt(ret)

if __name__ == '__main__':
  print(h_z_t_w())

