def array1_coverage(array1, array2, range_val=3):
    i, j = 0, 0
    coverage = [False] * len(array1)  # 初始化结果为 False

    while i < len(array1):
        # 检查 array2[j] 是否在 array1[i] 的上下游范围内
        while j < len(array2) and array2[j] < array1[i] - range_val:
            # 如果 array2[j] 小于 array1[i] 的左侧范围，则右移 j
            j += 1

        if j < len(array2) and array2[j] <= array1[i] + range_val:
            # 如果 array2[j] 在 array1[i] 的范围内，则标记覆盖为 True
            coverage[i] = True

        # 无论是否覆盖，继续处理下一个 array1[i]
        i += 1

    return coverage

# 示例数据
array1 = [10, 11, 12, 30, 50, 70, 100]
array2 = [15, 20, 35, 55, 90]

# 调用函数
result = array1_coverage(array1, array2)
print("结果:", result)
