import math

L1 = 0.3
L2 = 0.3
r1 = 0.1
r2 = 0.1

coupling_factor = 0.99
M = coupling_factor*math.sqrt(L1 * L2)

L = [[L1, M], [M, L2]]
R = [[r1, 0.0], [0.0, r2]]
V = [24.0, 0.0]

for count1 in range(len(L)):
    if not L[count1][count1]:
        for count2 in range(count1+1, len(L)):
            if L[count2][count1]:
                L[count1], L[count2] = L[count2], L[count1]
                R[count1], R[count2] = R[count2], R[count1]
                V[count1], V[count2] = V[count2], V[count1]
                break

    if L[count1][count1]:
        for count2 in range(count1+1, len(L)):
            comm_factor = L[count2][count1]/L[count1][count1]
            for count3 in range(len(L[count1])):
                L[count2][count3] -= L[count1][count3]*comm_factor
                R[count2][count3] -= R[count1][count3]*comm_factor
            V[count2] -= V[count1]*comm_factor


for count1 in range(len(L)):
    print(L[count1])
print()

for count1 in range(len(R)):
    print(R[count1])
print()

print(V)
