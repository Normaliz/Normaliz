import random
import os
import time

def random_ineq(n,d,ineq,grad):
    f = open('random_poly_'+str(n)+'_'+str(d)+'.in','w')
    print(f.name)
    f.write(str(n)+'\n')
    f.write(str(d)+'\n')
    x=[[0 for x in range(d)] for y in range(n)]
    for i in range(n):
        for j in range(d-1):
            x[i][j] = random.randint(-ineq,ineq)
            f.write(str(x[i][j])+' ')
        x[i][d-1] = random.randint(1,grad)
        f.write(str(x[i][d-1]))
        f.write('\n')
        print(x[i])
    f.write('inequalities\n')
    f.write('1\n')
    f.write(str(d)+'\n')
    for i in range(d-1):
        f.write('0 ')
    f.write('1\n')
    f.write('grading')
    f.close()
    start = time.time()
    os.system('/home/math/rsieg/normaliz/normaliz/source/normaliz -x=20 -cr1 '+f.name)
    end = time.time()
    print('Runtime: '+str(end-start)+' sec')
    # return x

