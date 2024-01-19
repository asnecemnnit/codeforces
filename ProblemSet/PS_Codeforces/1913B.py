output_list = []


def solve(output_list):
    cnt0 = cnt1 = 0
    s = str(input())

    for c in s:
        if c == "0":
            cnt0 += 1
        else:
            cnt1 += 1

    idx = 0
    for c in s:
        if c == "0":
            if cnt1 > 0:
                cnt1 -= 1
                idx += 1
            else:
                break
        else:
            if cnt0 > 0:
                cnt0 -= 1
                idx += 1
            else:
                break
    output_list.append(len(s) - idx)
    return


t = int(input())

while t > 0:
    solve(output_list)
    t -= 1

for output in output_list:
    print(output)
