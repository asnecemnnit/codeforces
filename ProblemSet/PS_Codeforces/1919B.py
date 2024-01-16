output_list = []


def solve(output_list):
    ans = 0
    n = int(input())
    s = str(input())

    for c in s:
        if c == "+":
            ans += 1
        else:
            ans -= 1

    output_list.append(int(abs(ans)))


t = int(input())

while t > 0:
    solve(output_list)
    t -= 1

for output in output_list:
    print(output)
