output_list = []


def solve(output_list):
    lmax = 0
    hmin = 1e9
    ignore_list = []
    n = int(input())
    while n > 0:
        a, x = map(int, input().split())

        if a == 3:
            ignore_list.append(x)
        elif a == 1:
            lmax = max(lmax, x)
        elif a == 2:
            hmin = min(hmin, x)

        if lmax <= hmin:
            ans = hmin - lmax + 1
        else:
            ans = 0

        n -= 1

    if ans != 0:
        for val in ignore_list:
            if val >= lmax and val <= hmin:
                ans -= 1
    output_list.append(int(ans))


t = int(input())

while t > 0:
    solve(output_list)
    t -= 1

for output in output_list:
    print(output)
