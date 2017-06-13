from cStringIO import StringIO

class LineReader(object):
    """Adapt a line iterator into something that also supports read.

    Note that seek is not supported! The point of this is to use the least
    amount of memory, and seek would require keeping the entire content in
    memory so there's no point in being clever. Also note that interleaved calls
    to read and readline/iteration will likely drop some data, so don't do that.

    """

    def __init__(self, line_iter):
        self.line_iter = line_iter
        self.buf = None

    def read(self, size=-1):
        # buf is a persistent sliding window over the current line that is
        # advanced as read() is called successively. When exhausted, it's
        # reloaded from the next line. data holds the results of the current
        # read as we assemble it.
        data = StringIO()
        if size < 0:
            data.writelines(self.line_iter)
        while data.tell() < size:
            if self.buf is None or len(self.buf) == 0:
                try:
                    line = next(self.line_iter)
                except StopIteration:
                    break
                self.buf = buffer(line)
            remainder = size - data.tell()
            if remainder <= len(self.buf):
                data.write(self.buf[:remainder])
                self.buf = self.buf[remainder:]
            else:
                data.write(self.buf)
                self.buf = None
        return data.getvalue()

    # The intended use of this class is to allow fileinput.FileInput objects,
    # which only offer line iteration, to serve as input to Pandas parsers, some
    # of which require buffer-like objects that implement byte reads (e.g. the
    # 'c' text file parsing engine). As for the parsers that still require line
    # iteration, the Pandas internals seem to be a bit inconsistent at the
    # moment as to which methods are actually called. In the interest of maximum
    # compatibility, this class implements more methods around line iteration
    # than might appear to be necessary at first glance.

    def __iter__(self):
        return self.line_iter

    def next(self):
        return next(self.line_iter)

    def readline(self):
        return next(self.line_iter)
