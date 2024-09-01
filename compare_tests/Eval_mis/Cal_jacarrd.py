import sys

'''
J(A,B) = (A∩B) / (AUB)
'''
comm = sys.argv[1]
a = sys.argv[2]
b = sys.argv[3]
comm, a, b = int(comm), int(a), int(b)
def cal_jaccarrd(comm, a, b):
    return comm / (comm + a + b)

print(cal_jaccarrd(comm, a, b))

def cal_jaccarrd2(comm1, comm2, a, b):
    '''
    将comm1和a放缩为comm2和a_
    '''
    if comm1 == 0 or comm2 == 0: return 0
    a_ = (comm2 / comm1) * a
    return comm2(comm2 + a_ + b)
