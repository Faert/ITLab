import numpy as np
import matplotlib.pyplot as plt

with open("Shor_result.txt") as Shor_result:
    all_result = [row.strip() for row in Shor_result]

result = all_result[1].split()
result = [int(i) for i in result]

#size, start, end, count
info_condition_exp = all_result[0].split()
info_condition_exp = [int(i) for i in info_condition_exp]

#a, N, error, res
info_Shor = all_result[2].split()
info_Shor = [int(i) for i in info_Shor]

#plt.plot(result)
a = [i for i in range(1 << (info_condition_exp[0]-2) // 2)]
plt.bar(a, result)
title = 'Shor(' + str(info_Shor[0]) + ', ' + str(info_Shor[1]) + ') with error = ' + str(info_Shor[2]) + '%'
plt.title(title)
plt.show()
