import math
from scipy.special import comb
from scipy.stats import norm

L = 256
T = 256
W = 128
ALPHA = 150.99775853018372
SIGMA_BASE = 9.553147931219351
CM = (L//32)**2
LB = 160000
UB = 200000

def p_t_w(t, w):
  return comb(t-1, w-1, exact=True) * pow(2, -L)

def p_y_t_w(y, t, w):
  mu = (t-1)*(3*ALPHA+11) + w*(3*ALPHA + 15/2)
  sigma = 3*SIGMA_BASE*math.sqrt(t+w-1)
  return norm.pdf(y, loc=mu, scale=sigma)

def p_z_t_w(z, t, w):
  return p_y_t_w(z/CM, t, w)

def p_z(z):
  if z % CM != 0:
    return 0
  ret = 0
  for t in range(L, L//2, -1):
    for w in range(t, 0, -1):
      ret += (p_z_t_w(z, t, w) * p_t_w(t, w))
  return ret

# 11.176580436192054
# def h_z_t_w():
#   ret = 0
#   for t in range(L, L//2, -1):
#     for w in range(t, 0, -1):
#       sigma = 3*SIGMA_BASE*math.sqrt(t+w-1)
#       ret += (p_t_w(t, w) * math.log2(2*math.pi*math.e*(sigma**2)) / 2 )
#     print(ret)
#   return ret


# 11.176424250835314
# t = 232, w = 73
# 11.176424250835314
# t = 231, w = 231
def h_y_t_w():
  ret = 0
  for t in range(L, L//2, -1):
    switch = False
    for w in range(t, 0, -1):
      print("t = {0}, w = {1}".format(t, w))
      entropy = 0
      for y in range(LB, UB, 1):
        tmp = p_y_t_w(y, t, w)
        if tmp != 0:
          entropy += (-tmp * math.log2(tmp))
      tmp2 = (p_t_w(t, w) * entropy)
      if tmp2 != 0:
        switch = True
        ret += tmp2
      elif switch and tmp2 == 0:
        break
      print(ret)
    print(ret)
  return ret

# 13.953367137522822
def h_z():
  # 9673984
  # ret = 4.3574652903916986e-08
  # start: 160000*CM, end: 200000*CM, CM=64、初期ret=0の場合、
  # ---------------------------
  # ret = 13.953367070591918
  # 12727296
  # ---------------------------

  ret = 0
  for z in range(LB*CM, UB*CM, CM):
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
      ret += (mu * p_t_w(t, w))
      print(ret)
  return ret

# 560.1436784349372
def calc_sigma():
  ret = 0
  for t in range(L, 0, -1):
    for w in range(t, 0, -1):
      sigma = 3*SIGMA_BASE*math.sqrt(t+w-1)
      ret += ((sigma**2) * p_t_w(t, w))
      print(math.sqrt(ret))
  return math.sqrt(ret)

if __name__ == '__main__':
  print(h_y_t_w())


