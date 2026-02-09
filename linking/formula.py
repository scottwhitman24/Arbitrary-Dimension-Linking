import math


def vol(n,m):
    term1=1
    term2=1
    if m%2==0:
        term1=math.pi/2
    if m%2==1:
        term1=1
    if n%2==0:
        term2=math.pi/2
    if n%2==1:
        term2=1
    print("m")
    print(m)
    for i in range(1,m//2):
        print(m-2*i)
        term1=((m-2*i-1)/(m-2*i))*term1
    for i in range(1,n//2):
        print(n-2*i)
        term2=((n-2*i-1)/(n-2*i))*term2
    print("terms")
    print(term1)
    print(term2)
    return 4*spherearea(m-1)*spherevol(n-m)*(term1/m+term2/n)-4/n*spherearea(m-1)*spherevol(n-m+1)

def voleven(n,m):
    k=n-m+1
    term1=sinint(m-2)/m+sinint(n-2)/n
    term2=0
    for i in range(1,k-1):
        term2+=sinint(2*(k-i)+m-1)/((2*(k-i)+1)*(2*(k-i)+m-1))*P(i,2*k)
    term3=0
    for i in range(1,k-3):
        term3+=sinint(2*(k-2-i)+m-1)/((2*(k-2-i)+1)*(2*(k-2-i)+m-1))*P(i,2*(k-2))
    term4=(math.pi/4+7/6)*P(k-3,2*(k-2))/(2*k)
    '''print("|")
    print(term1)
    print(term2)
    print(term3)
    print(term4)
    print("|")'''
    return 4*spherearea(m-1)*spherevol(k-1)*(term1+2*(term2-term3)+term4)
def spherevol(n):
    result=1
    if n%2==0:
        result=1
    else:
        result=2
    for i in range(n//2):
        result*=2*math.pi/(n-2*i)
    return result
def spherearea(n):
    return n*spherevol(n)

def P(i,n):
    prod=1
    for j in range(i):
        prod*=(n-2*j-1)/(n-2*j)
    return prod
def sinint(n):
    if n%2==0:
        return P(n//2-1,n)*math.pi/4
    else:
        return P(n//2-1,n)*math.pi/2



'''for m in range(3,19):
    print(str(m)+"|",end='')
print("\n")
for n in range(3,20):
    for m in range(2,n-1):
        print("n="+str(n))
        print("m="+str(m))
        print(str(vol(n,m)))
    print("\n")

print(1+3+1+2+2+2+1+3+3+2+2+3+2+5+1+3+4+2+1+1+1+2+4+5+2+2+5+4+2+4+1+1+0+4+5+4+2+6+0+2+3+1+3+3+2+1+1+5+1+1+4+2+1+3+2+6+1+1+0+0+1+3+2+2+3+2+4+2+2+2+1+1+5+2+3+1+1+1+1+1+1+1+2+6+1+3+1+3+0+2+0+6+5)
'''
for n in range(3,20):
    print(str(n)+":",end="\t")
    for m in range(1,n-1):
        if (n-m+1)%2==0:
            print(str(m)+"->",end="")
            print(voleven(n,m),end="\t")
    print()
for n in range(3,50):
    print(str(n)+":",end="\t")
    for m in range(1,n-1):
        if (n-m+1)%2==0:
            print(str(m)+"->",end="")
            print(voleven(n,m)/(2**n*spherevol(n)),end="\t")
    print()