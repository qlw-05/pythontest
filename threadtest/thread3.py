from threading import Event, current_thread, Thread
import time
event = Event()


def check():
    print('%s 正在检测服务是否正常....' % current_thread().name)
    time.sleep(5)
    event.set()  # success由false改成了True


def connect():
    print('%s 等待连接...' % current_thread().name)
    count = 0
    flag = True
    while flag:
        event.wait(1)  # 取代了上面注释掉的代码
        count += 1
        if event.is_set():
            flag = False
            print('%s 开始连接...' % current_thread().name)
        else:
            if count >= 3:
                flag = False
                print('%s 三次连接失败' % current_thread().name)
            else:
                print('%s 一次连接失败' % current_thread().name)


if __name__ == '__main__':
    t1 = Thread(target=connect)
    t2 = Thread(target=connect)
    t3 = Thread(target=connect)
    c1 = Thread(target=check)
    t1.start()
    t2.start()
    t3.start()
    c1.start()
