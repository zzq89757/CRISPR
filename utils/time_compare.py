import pandas as pd
import numpy as np
import time

# 创建示例 DataFrame
df = pd.DataFrame({
    "col1": np.random.randint(0, 1000, 10**6),
    "col2": np.random.randint(0, 1000, 10**6),
    "col3": np.random.randint(0, 1000, 10**6),
})

# 测试 iloc 方法
start = time.time()
for idx in range(1000):
    _ = df.iloc[idx]
print("iloc time:", time.time() - start)

start = time.time()
# 转换为 NumPy 数组
col1 = df["col1"].to_numpy()
col2 = df["col2"].to_numpy()
col3 = df["col3"].to_numpy()

# 测试 NumPy 方法
for idx in range(1000):
    _ = (col1[idx], col2[idx], col3[idx])
print("NumPy time:", time.time() - start)
