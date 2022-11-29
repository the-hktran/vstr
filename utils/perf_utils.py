import time
import numpy as np


class TIMER:
    """ Class for a simple timer.

    Examples:
    >>> timer = TIMER(3)
    >>> # set overhead
    >>> ovhd = timer.estimate_overhead()
    >>> timer.set_overhead(ovhd)
    >>> # use timers
    >>> for i in range(10):
    >>>     timer.start(0)  # start 1st timer
    >>>     func1(...)
    >>>     timer.stop(0)   # stop 1st timer
    >>>     timer.start(1)
    >>>     ...
    >>> # read timing
    >>> print("Timer 0: {:.3f}" % (timer.read(0)))
    >>> # formatted timing report
    >>> timer_names = ["func1", "fancy_func", "func3"]
    >>> timer.report(timer_names)
    """
    def __init__(self, n=None, name=None, stdout=None,
        logger_=None, logger_last_=None, verbose=3):

        if n is None: n = 1
        if n <= 0 or not isinstance(n,int):
            raise ValueError("n (# of timers) must be a positive integer.")
        self.n = n
        self.name = "" if name is None else name

        self.__overhead = 0. # set it using method "set_overhead"

        self.reset()

        if stdout is None:
            import sys
            self.stdout = sys.stdout
        else:
            self.stdout = stdout
        self.verbose = verbose
        #self.logger = logger.note if logger_ is None else logger_
        #self.logger_last = self.logger if logger_last_ is None else logger_last_

    def reset(self):
        self.__tstart = [0.] * self.n
        self.__telapsed = [0.] * self.n
        self.__telapsed_last = [0.] * self.n

    def start(self, i):
        self.__tstart[i] = time.time()

    def stop(self, i):
        self.__telapsed_last[i] = time.time() - self.__tstart[i] - \
            self.__overhead
        self.__telapsed[i] += self.__telapsed_last[i]

    def add(self, i, dt):
        self.__telapsed[i] += dt

    def read(self, i, last=False):
        if last:
            return self.__telapsed_last[i]
        else:
            return self.__telapsed[i]

    def read_tot(self, last=False):
        if last:
            return np.sum(self.__telapsed_last)
        else:
            return np.sum(self.__telapsed)

    def estimate_overhead(self, nrepeat=None):
        if nrepeat is None: nrepeat = int(1E6)
        timer = TIMER(n=1)
        a = 0
        for i in range(nrepeat):
            a += 1
            timer.start(0)
            timer.stop(0)

        return timer.read(0) / nrepeat

    def set_overhead(self, overhead):
        self.__overhead = overhead

    def get_overhead(self):
        return self.__overhead

    def report(self, tnames=None, trange=None, last=False, comments=None):
        if tnames is None:
            tnames = ["Timer %-2d" % i for i in range(self.n)]
        if trange is None:
            trange = list(range(self.n))
        assert(len(tnames) == len(trange))
        #logger_ = self.logger_last if last else self.logger
        if not comments is None:
        #    logger_(self, "%s" % comments)
            print("%s" % comments)
        ts = [self.read(i,last=last) for i in trange]
        t_tot = np.sum(ts)
        for ti,tname in zip(ts,tnames):
            #logger_(self, "  %20s  :::  %11.3f (%6.2f%%)",
            #    tname.ljust(20), ti, ti/t_tot*100.)
            print("  %20s  :::  %11.3f (%6.2f%%)" % (tname.ljust(20), ti, ti/t_tot*100.))

        total_name = "Total " + self.name
        #logger_(self, "  %20s  :::  %11.3f (100.00%%)\n",
        #    total_name.ljust(20), t_tot)
        print("  %20s  :::  %11.3f (100.00%%)\n" % (total_name.ljust(20), t_tot))


class TIMER_:
    """Void timer.
    """
    def __init__(self):
        pass
    def start(self, x):
        pass
    def stop(self, x):
        pass


def check_mem(obj, thr=1E-3):
    """Only print items > thr [MB]
    """
    from pympler.asizeof import asizeof as asz
    from frankenstein.tools.io_utils import prtvar
    prtvar("Object class", str(obj.__class__), "{:s}")
    hasdict = hasattr(obj, "__dict__")
    if hasdict:
        sBdict_sig = dict()
        thrB = thr*1024**2.
        sBtot = 0
        for key in obj.__dict__:
            s = asz(obj.__dict__[key])
            sBtot += s
            if s > thrB:
                sBdict_sig[key] = s
        sorted_sB_sig = sorted(sBdict_sig.items(), key=lambda kv: kv[1])

        for kv in sorted_sB_sig:
            key, sB = kv
            sMB = sB / 1024**2.
            prtvar(key, ("%15d B"%sB).ljust(17) + "  " +
                ("%15.2f MB"%(sMB)).ljust(18), "{:s}")
        prtvar("Total", ("%15d B"%sBtot).ljust(17) + "  " +
            ("%15.2f MB"%(sBtot/1024**2.)).ljust(18), "{:s}")
        print(flush=True)
    else:
        sB = asz(obj)
        sMB = sB / 1024**2.
        prtvar("Total", ("%15d B"%sB).ljust(17) + "  " +
            ("%15.2f MB"%(sMB)).ljust(18), "{:s}")
        print(flush=True)


def current_memory(mem0=None):
    """ Return current memory usage in MB
    """
    if mem0 is None: mem0 = 0
    import os
    import psutil
    p = psutil.Process(os.getpid())
    return p.memory_info()[0] / 1024**2. - mem0
