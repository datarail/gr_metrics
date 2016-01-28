from cStringIO import StringIO

import q

class LineReader(object):
    """Adapt a line iterator into something that supports read.

    Note that seek is not supported! The point of this is to use the least
    amount of memory, and seek would require keeping the entire content in
    memory so there's no point in being clever."""

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
