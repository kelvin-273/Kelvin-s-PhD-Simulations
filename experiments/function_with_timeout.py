from time import sleep
import asyncio
import eugene.utils as eu

import multiprocessing as mp


def target_func(queue, *args, **kwargs):
    try:
        result = your_function(*args, **kwargs)
        queue.put(("ok", result))
    except Exception as e:
        queue.put(("err", e))


def run_with_timeout(func, timeout=300, *args, **kwargs):
    queue = mp.Queue()
    p = mp.Process(target=target_func, args=(queue, *args), kwargs=kwargs)
    p.start()
    p.join(timeout)

    if p.is_alive():
        p.terminate()
        p.join()
        return "timeout"

    if not queue.empty():
        status, value = queue.get()
        if status == "ok":
            return value
        else:
            raise value

    return None
