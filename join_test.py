import numpy as np
import time

arr = np.random.randint(0, 1000, size=(1000000,))

# 方法 1: astype
start = time.time()
"\t".join(arr.astype(str))
print("astype:", time.time() - start)

# 方法 2: generator
start = time.time()
"\t".join(str(x) for x in arr)
print("generator:", time.time() - start)

# 方法 3: np.char
start = time.time()
np.char.join("\t", arr.astype(str))
print("np.char:", time.time() - start)
