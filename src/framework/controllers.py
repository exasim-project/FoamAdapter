class _CtrlState:
    def __init__(self, outer_iter, correct_iter):
        self.outer_iter = outer_iter
        self.correct_iter = correct_iter

class PIMPLE:
    def __init__(self, max_outer, max_correct):
        self.max_outer = max_outer
        self.max_correct = max_correct
    def __call__(self):
        for outer in range(self.max_outer):
            for correct in range(self.max_correct):
                yield _CtrlState(outer, correct)

class FixedOuter:
    def __init__(self, iterations):
        self.iterations = iterations
    def __call__(self):
        for outer in range(self.iterations):
            yield _CtrlState(outer, 0)
