#include <vector>
#include <string>
#include <unordered_map>
#include <bitset>
#include <iostream>
using namespace std;
class Solution
{
public:
    std::vector<string> findRepeatedDnaSequences(string s)
    {   // 10bp 不同碱基组合为0到(1 << 20 - 1) 因此建立该长度的binset 出现了就将对应位置设置为1 之后再test test发现重复后放入s2 后续遇到时跳过
        // 对应二进制00, 01, 10, 11.那么10个组合只要20位就够了。
        unordered_map<char, int> m{{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
        vector<string> res;
        bitset<1 << 20> s1, s2;            // 那么所有组合的值将在0到(1 << 20 - 1)之间
        int val = 0, mask = (1 << 20) - 1; // mask等于二进制的20个1
        // 类似与滑动窗口先把前10个字母组合
        for (int i = 0; i < 10; ++i)
            val = (val << 2) | m[s[i]];
        s1.set(val); // 置位 不同的10bp的十进制数会把s1对应的位置设为1
        for (int i = 10; i < s.size(); ++i)
        {
            val = ((val << 2) & mask) | m[s[i]]; // & mask 操作保留 val 的后 20 位 去掉左移的一个字符再加上一个新字符
            if (s2.test(val))
                continue; // 出现过两次跳过
            if (s1.test(val))
            {
                res.push_back(s.substr(i - 9, 10));
                s2.set(val);
            }
            else
                s1.set(val);
        }
        return res;
    }
};

int main()
{
    char const *s = "AAAAACCCCCAAAAACCCCCCAAAAAGGGTTT";
    Solution ss;
    vector<string> res = ss.findRepeatedDnaSequences(s);
    for (int i = 0; i < res.size(); i++)
    {
        cout << res.at(i) << "\n";
    }
}