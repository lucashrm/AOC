import os
from time import sleep

def myfunc():
    fifo_fn = r'\\.\pipe\ixon_fifo_in'
    try:
        fd = os.open(fifo_fn, os.O_WRONLY)
        res = os.write(fd, b'quit')
    except:
        print("could not open the fifo in write mode (%s)" %(fifo_fn))

myfunc()