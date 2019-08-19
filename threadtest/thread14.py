from multiprocessing import Process
import time


# 计算密集型
def task1():
    time.sleep(1)


def task2():
    time.sleep(2)


if __name__ == '__main__':
    start = time.time()
    task1()
    task2()
    stop = time.time()
    print('串行执行的时间：%s' % (stop - start))
    p1 = Process(target=task1,)
    p2 = Process(target=task2,)
    start = time.time()
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    stop = time.time()
    print('多进程执行的时间：%s' % (stop-start))
