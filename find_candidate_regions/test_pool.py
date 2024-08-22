from multiprocessing import Pool, Queue, Manager
import os
import time
 
 
def compiler(q1, q2):
    """
    从q1队列中获取输入 -->  执行编译，得到结果 -->  将推理命令行加入q2队列
    :param q1: compile queue
    :param q2: infer queue
    :return: call back info
    """
    # 当一个队列为空的时候如果再用get取则会堵塞
    # 也可以使用q.get_nowait() ==相当于== q.get(False)
    compile_codeline = q1.get_nowait()  # 删除队列中的item，并返回
 
    print("subprocess ID {} is compiling".format(os.getpid()))
    time.sleep(5)
    infer_codeline = "infer" + str(compile_codeline)
 
    q2.put(infer_codeline)
    res_info = str(compile_codeline) + " complete compile."
    return res_info
 
 
def infer(q2):
    infer_codeline = q2.get(True)
    print("subprocess ID {} is infering".format(os.getpid()))
    time.sleep(1)
    res_info = str(infer_codeline) + " complete infer."
    return res_info
 
 
def call_back(res):
    print(res)
 
 
def error_call_back(error_code):
    print(error_code)
 
 
def multiProcess():
    """
    1. 异步非阻塞式， pool.apply_async()
    不用等待当前进程执行完毕，随时跟进操作系统调度来进行进程切换，即多个进程并行执行，提高程序的执行效率。
    2. BUG
    Queue objects should only be shared between processes through inheritance 队列只能通过集成在进程之间共享
    解决方法：如果想要在进程池中使用队列，需要使用mutiprocessing的Manager类
    """
    manager = Manager()
    q1 = manager.Queue()
    q2 = manager.Queue()
 
    cases = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for i in cases:
        q1.put(i)
    case_number = q1.qsize()
    print("the number of compile model is {}".format(case_number))
 
    compile_pool = Pool(processes=4)
    for i in range(case_number):
        compile_pool.apply_async(compiler, args=(q1, q2,), callback=call_back, error_callback=error_call_back)
    compile_pool.close()  # 不能再继续添加新的子进程
    compile_pool.join()  # 同步函数
 
    infer_pool = Pool(1)
    infer_number = q2.qsize()
    print("the number of infer model is {}".format(infer_number))
    for i in range(infer_number):
        infer_pool.apply_async(infer, args=(q2,), callback=call_back, error_callback=error_call_back)
    infer_pool.close()
    infer_pool.join()
 
 
multiProcess()
 