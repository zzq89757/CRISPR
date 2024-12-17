class Solution:
    def findRepeatedDnaSequences(self, s: str):
        if len(s) < 10:
            return []

        m = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N':4}
        res = []
        seen_once = set()
        seen_twice = set()
        val = 0
        mask = (1 << 3) - 1  # mask等于二进制的3个1

        # 类似于滑动窗口先把前2个字母组合
        for i in range(2):
            val = (val << 3) | m[s[i]]
        seen_once.add(val)  # 置位

        # 开始滑动窗口
        for i in range(2, len(s)):
            val = ((val << 3) & mask) | m[s[i]]  # 去掉左移的一个字符再加上一个新字符
            if val in seen_twice:
                continue  # 出现过两次跳过
            if val in seen_once:
                res.append(s[i - 1:i + 1])
                seen_twice.add(val)
            else:
                seen_once.add(val)
        print(seen_once)
        return res

if __name__ == "__main__":
    s = Solution()
    s.findRepeatedDnaSequences("NNCGATCGATGTCGATGCTGATGCTG")  