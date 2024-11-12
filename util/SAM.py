#!/usr/bin/env python3
import copy


class sstate(object):
    def __init__(self):
        self.len = 0
        self.parent = -1
        self.next = dict()
        self.pos = 0


class SAM(object):
    def __init__(self, str=''):
        self.size = 0
        self.states = list()
        self.last = 0
        state0 = sstate()
        self.states.append(state0)
        if str:
            for c in str:
                self.extend(c)

    def extend(self, c):
        self.size += 1
        cur = self.size
        state = sstate()
        state.len = self.states[self.last].len + 1
        state.pos = state.len - 1
        self.states.append(state)
        up = self.last

        while up != -1 and (c not in self.states[up].next):
            self.states[up].next[c] = cur
            up = self.states[up].parent

        if up == -1:
            self.states[cur].parent = 0
        else:
            q = self.states[up].next[c]
            if self.states[up].len + 1 == self.states[q].len:
                self.states[cur].parent = q
            else:
                self.size += 1
                clone = self.size
                clone_state = sstate()
                clone_state.len = self.states[up].len + 1
                clone_state.next = copy.deepcopy(self.states[q].next)
                clone_state.parent = self.states[q].parent
                clone_state.pos = self.states[q].pos
                self.states.append(clone_state)

                while up != -1 and self.states[up].next[c] == q:
                    self.states[up].next[c] = clone
                    up = self.states[up].parent

                self.states[q].parent = clone
                self.states[cur].parent = clone

        self.last = cur

    def traverse(self, str):
        if not str:
            return False, -1

        cur = 0
        for c in str:
            if c in self.states[cur].next:
                cur = self.states[cur].next[c]
            else:
                return False, -1
        return True, self.states[cur].pos - len(str) + 1
