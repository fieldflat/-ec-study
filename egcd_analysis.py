import math
import random
from tqdm import tqdm
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# =============================
# global constants
# =============================
LOOP = 10000
MODULO_P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f  # 256bit
ITERATION = 255+127-1

# =============================
# variables for analysis
# =============================
count = 0
count_list = []

# Extended Euclidean algorithm
def egcd(a: int, b: int):
  global count
  (x, lastx) = (0, 1)
  (y, lasty) = (1, 0)
  while b != 0:
    count += 3
    q = a//b
    (a, b) = (b, a - q*b)
    (x, lastx) = (lastx - q*x, x)
    (y, lasty) = (lasty - q*y, y)
  return (lastx, lasty, a)

if __name__ == '__main__':
  for _ in tqdm(range(LOOP)):
    for _ in range(ITERATION):
      a = random.randint(1, MODULO_P)
      b = random.randint(1, MODULO_P)
      egcd(a, b)
    count_list.append(count)
    count = 0
  
  # display histgram
  plt.hist(count_list, bins=len(set(count_list)), density=True)

  # normal distribution
  # average is ok
  alpha = (12*math.log(2)*math.log(2**256))/(math.pi**2) + 0.06535

  mu = 3*alpha*ITERATION
  sigma = 3*9.5*math.sqrt(ITERATION)

  # mu = alpha*ITERATION
  # sigma = 9.5*math.sqrt(ITERATION)

  print("LOOP = {0}".format(LOOP))
  print("MODULO_P = {0}".format(MODULO_P))
  print("ITERATION = {0}".format(ITERATION))
  print("alpha = {0}".format(alpha))
  print("mu = {0}".format(mu))
  print("sigma = {0}".format(sigma))
  print("Average of Multipliucation (experimentally) = {0}".format(sum(count_list)/len(count_list)))
  print("Average of Multipliucation (theoretical) = {0}".format(mu))

  X = np.arange(mu-sigma*5, mu+sigma*5, 1)
  Y = norm.pdf(X, loc=mu, scale=sigma)
  plt.plot(X, Y, 'r-')
  plt.savefig("egcd_analysis.png")

  plt.show()

