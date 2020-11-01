import math
import random
from tqdm import tqdm
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# =============================
# global constants
# =============================
LOOP = 100000
MODULO_P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f  # 256bit

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
    count += 1
    q = a//b
    (a, b) = (b, a - q*b)
    (x, lastx) = (lastx - q*x, x)
    (y, lasty) = (lasty - q*y, y)
  return (lastx, lasty, a)

if __name__ == '__main__':
  for _ in tqdm(range(LOOP)):
    a = random.randint(1, MODULO_P//100)
    b = random.randint(1, MODULO_P//100)
    egcd(a, b)
    count_list.append(count)
    count = 0
  
  # display histgram
  plt.hist(count_list, bins=len(set(count_list)))

  # normal distribution
  mu = (12*math.log(2)*math.log(2**255))/(math.pi**2) + 1.467
  sigma = 9.5

  X = np.arange(100, 200, 0.1)
  Y = norm.pdf(X, loc=mu, scale=sigma)
  # プロットする
  # fig = plt.figure()
  # ax = fig.add_subplot(1, 1, 1)
  # ax.grid(color='gray')
  # ax.plot(X, Y, color='blue')
  plt.plot(X, Y*LOOP, 'r-')

  plt.show()

