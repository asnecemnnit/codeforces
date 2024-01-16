from math import gcd as GCD

output_list = []


def LCM(a, b, gcd):
    return int(a * b / gcd)


def solve(output_list):
    a, b = map(int, input().split())

    gcd = GCD(a, b)
    lcm = LCM(a, b, gcd)

    if lcm > b:
        ans = lcm
    else:
        a = a / gcd
        b = b / gcd

        ans = gcd * b * b

    output_list.append(int(abs(ans)))


t = int(input())

while t > 0:
    solve(output_list)
    t -= 1

for output in output_list:
    print(output)
