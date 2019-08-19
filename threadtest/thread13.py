from threading import Thread
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
    t1 = Thread(target=task1,)
    t2 = Thread(target=task2,)
    start = time.time()
    t1.start()
    t2.start()
    t1.join()
    t2.join()
    stop = time.time()
    print('多线程执行的时间：%s' % (stop-start))
